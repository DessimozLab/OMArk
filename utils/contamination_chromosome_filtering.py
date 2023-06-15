import gffutils
from Bio import SeqIO
import os
import sys
import argparse

def build_arg_parser():
    """Handle the parameter sent when executing the script from
 the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="This script use the OMArk contaminant list to identify corresponding areas in a genome and filter these proteins out. It is designed to reduce false positive and false negative by exploiting genomic coordinate. The main idea is that contaminant genes would be from similar contigs or part of the chromsome, missassembled to host Chromosomes, thus grouped tohgether.")
    parser.add_argument('-i', '--omark_folder', type=str, help="OMArk output folder for target proteome.", required=True)
    parser.add_argument('-g', '--gff', type=str, help="Path to the GFF file corresponding to target genome, with location of genes and CDS. Compatible with Ensembl and NCBI format. Please raise an issue to https://github.com/DessimozLab/OMArk/ if you encounter difficulties with your own GFF", required=True)
    parser.add_argument('-f', '--fasta', type=str, help="Path to the original FASTA file from target proteome.",  required=True)
    parser.add_argument('-o', '--output', type=str, help="Prefix path to output files. This script will generate a report file with selected chromsomal stretches and list of genes and proteins; as well as FASTA file with removed contaminants",required=True)

    parser.add_argument('-t', '--threshold', type=float,  help="The threshold (float between 0 and 1) that the proportion of genes in a 'contaminant stretch' must pass to be considered as valid", default=0.5)
    parser.add_argument('-n', '--minimal_gene_nr', type=int,  help="The minimum number of genes in a contaminant stretch to be considered as valid. If the number of genes in a chromosome is lower than this, all of the genes in it must be contaminant to be considered as stretch.", default=5)
    parser.add_argument('-x', '--only_extrem', help="Add if only full-contig stretches or stretches at extremities of contig should be considered as valid",  action='store_true')

    return parser


#Read an OMArk .ump file and extract the contaminant proteins from there
#omark_categ : path to the OMArk .ump file
# return contaminant_ids : a list of protein identifier
def get_contaminants(omark_folder):
    omark_categ = None
    for file in os.listdir(omark_folder):
        if file[-4:]=='.ump':
            omark_categ = os.path.join(omark_folder, file)

    if not omark_categ:
        sys.exit('An .ump file was not found in your folder')
    conta = False
    contaminant_ids = []
    with open(omark_categ) as f:
        for line in f.readlines():
            line = line.replace('\n','')

            if line[0]=='>':
                if line[1:].split('_')[0]=='Contamination':
                    conta=True
                else:
                    conta=False
            elif conta:
                contaminant_ids.append(line)
    return contaminant_ids

#Read a GFF file and find the genomic location of all contaminant proteins in the file
#Inputs:
#contaminant_list : a list of protein id, describing the contaminant proteins
#gff: path to the input GFF file. This script has been tested on and is compatible with Ensembl and NCBI GFF3.
#Outputs: a tuple of three values
#contaminant_position: a list of the positions of contaminant genes in the genome - An item is a tuple of (chromosome_id, start position, end position, gene identifier)
#all_positions: a list of the positions of all the genes in the genome - An item is a tuple of (chromosome_id, start position, end position, gene identifier)
#gene_to_prot : a dictionary linking all gene identifiers (key) to all of their protein identifiers (values as a list)
def get_position_conta(contaminant_list, gff):
    gene_to_prot = {}
    contaminant_positions = []
    all_positions = []
    contaminant_list = [x.split('.')[0] for x in contaminant_list]
    db = gffutils.create_db(gff, ':memory:', merge_strategy="create_unique", keep_order=True)
    for t in db.features_of_type('gene', order_by='start'):
        gene_id = t.id
        gene_id = gene_id.replace('gene:', '')
        gene_is_cont = False
        ordered_child = list(db.children(t, featuretype='CDS', order_by='start'))
        for child in ordered_child:
            type_attribute = ['protein_id', 'Name']
            for att_type in type_attribute:
                    protein = child.attributes.get(att_type, [None])[0]
                    if protein:
                        break
            if protein and protein not in  gene_to_prot.get(gene_id,[]):
                protein_prefix =  protein.split('.')[0]
                prot_in_gene = gene_to_prot.get(gene_id,[])
                prot_in_gene.append(protein)
                gene_to_prot[gene_id] = prot_in_gene
                if protein_prefix in contaminant_list:
                    gene_is_cont = True
        if gene_is_cont:
            contaminant_positions.append((t.chrom, t.start, t.end, gene_id))
        all_positions.append((t.chrom, t.start, t.end, gene_id))
    return contaminant_positions, all_positions, gene_to_prot

    #Function to find stretch of the genomes with high representation of contaminants and output them as a list
#Parameters:
#cont_pos :  list of the positions of contaminant genes in the genome - An item is a tuple of (chromosome_id, start position, end position, gene identifier)
#all_pos: a list of the positions of all the genes in the genome - An item is a tuple of (chromosome_id, start position, end position, gene identifier)
#Thresh (float): the threshold  (between 0 and 1) that the proportion of genes in a "contaminant stretch" must pass to be considered as valid
#min_number_in_stretch (int): The minimum number of genes in a contaminant stretch to be considered as valid. If the number of genes in a chrosome is lower than this, all of the genes in it must be contaminant to be considered as stretch
#force_extremities: a boolean that indicates whether stretch can be only at extremities or also in middle of chromosomes. When true, force the stretches to be at the start or end of chromosomes
#Output:
#stretch_list: a list of contaminant chromosomal stretch. Each entry is a tuple of (chromosome, start position, end position, number of contaminant genes in the stretch, number of genes in the stretch, fraction of contaminant genes in a the stretch)
def infer_contaminant_genome_stretches(cont_pos,all_pos, thresh=0.5, min_number_in_stretch=5, force_extremities=True):

    sorted_possible = sorted(all_pos, key=lambda element: (element[0], element[1]))
    #List of contaminant strech in th genome
    stretch_list  = []
    #Threshold of minimal proportion of gene in a stretch
    thresh = thresh
    #Minimal number of contaminant protein in a stretch
    min_number_in_stretch = 5
    #Strech of contamination start and end
    cur_start = None
    cur_end = None
    #Number of contaminant in stretch
    pos_in_stretch = 0
    #Total number of gene in a stretch
    gene_in_stretch = 0
    prev_chrom = None
    in_strech = False
    #Number of contaminant genes seen from the last begining of a stretcj
    pos_from_last =0
    #Number of gene seen from the start of the chromosome
    from_start = 9999
    #Number of genes seen from the start of the stretch
    from_last = 9999
    for elem in sorted_possible:
        chrom = elem[0]
        start = elem[1]
        end = elem[2]
        if chrom != prev_chrom:

            in_stretch = False
            if (pos_in_stretch>= min_number_in_stretch) or (pos_from_last/from_start ==1):
                if pos_from_last/from_last >= thresh:
                    fraction_gene = pos_in_stretch/gene_in_stretch
                    stretch_list.append((prev_chrom, cur_start, -1, pos_in_stretch ,gene_in_stretch, fraction_gene))
                else:
                    fraction_gene = pos_in_stretch/gene_in_stretch
                    if not force_extremities or (cur_end==-1 or cur_start==1):
                        stretch_list.append((prev_chrom, cur_start, cur_end, pos_in_stretch ,gene_in_stretch, fraction_gene))

            prev_chrom = chrom
            start_chrom = start
            from_start = 0
            from_last = 0
            without_pos = 0
            pos_from_last =0 
            pos_in_stretch = 0
            gene_in_stretch = 0
            cur_start = start
            cur_end = end
        from_start += 1
        from_last += 1
        
        if elem in cont_pos:

            pos_from_last +=1 
            #Differentiate from start to start of a new stretch


            if pos_from_last/from_start >= thresh:
                if not in_stretch:
                    cur_start = 1
                    cur_end = end
                    in_stretch =True
                    gene_in_stretch=from_start
                else:
                    cur_end = end
                    gene_in_stretch+=without_pos+1

                pos_in_stretch += 1            
            elif pos_from_last/from_last >= thresh:

                cur_end = end
                pos_in_stretch+= 1
                gene_in_stretch+=without_pos+1
            else:

                if in_stretch:
                    if pos_in_stretch>= min_number_in_stretch:
                        if not force_extremities or (cur_end==-1 or cur_start==1):
                            fraction_gene = pos_in_stretch/gene_in_stretch
                            stretch_list.append((chrom, cur_start, cur_end, pos_in_stretch, gene_in_stretch,fraction_gene ))
                pos_in_stretch =1
                cur_start = start
                cur_end = end
                pos_from_last = 1
                from_last = 1
                gene_in_stretch = 1
            without_pos = 0
        else:
            without_pos+=1
    return stretch_list

#Get proteins and genes that are located within a given chromosomal stretch
#Parameter:
#contaminant_stretches: a list of tuple representing chromosomal stretches.  Each entry is a tuple of (chromosome identifer, start position, end position, number of contaminant genes in the stretch, number of genes in the stretch, fraction of contaminant genes in a the stretch)
#all_pos: a list of the positions of all the genes in the genome - An item is a tuple of (chromosome_id, start position, end position, gene identifier)
#gene_to_prot : a dictionary linking all gene identifiers (key) to all of their protein identifiers (values as a list)
#Outputs:
#prot_list: a list of protein identifiers corresponding to the chromosomal stretches
#gene_list: a list of gene identifiers and chromsomal position corresponding to the chromosomal stretches as a tuple with (gene identifier, chromosome_id, start position, end position)
def get_genes_in_cont_stretches(contaminant_stretches, all_pos,gene_to_prot):

    prot_list = []
    gene_list = []
    sorted_possible = sorted(all_pos, key=lambda element: (element[0], element[1]))
    for i, elem in enumerate(contaminant_stretches):
            for j, elem2 in enumerate(sorted_possible):

                if elem2[0]<elem[0]:
                    continue
                elif elem2[0]>elem[0]:
                    break
                else:
                    if elem2[0]==elem[0]:
                        if elem2[1]>=elem[1] and (elem2[2]<=elem[2] or elem[2]==-1):
                            gene_list.append(elem2)
                            prot_list +=gene_to_prot[elem2[3]]
    return prot_list, gene_list

#Read a FASTA and copy a version of it with some proteines removed
#Parameters:
#fasta: Path to the original FASTA file
#newfa: Path to the new FASTA file, to be written
#cont_proteins: a list of protein identifier to remove
def filter_proteins(fasta, newfa, cont_proteins):
    keep_record = []
    cont_list = cont_proteins
    cont_list = [x.split('.')[0] for x in cont_list]
        
    for record in SeqIO.parse(fasta, "fasta"):
        rec_id = record.id.split('.')[0]

        if rec_id not in cont_list:
            keep_record.append(record)
    SeqIO.write(keep_record, newfa, "fasta")

#Write a text file indicating which chromosomal stretches were deemed as contaminant, the corresponding genes and the corresponding proteins
#Parameters:
#contaminant_stretches: a list of tuple representing chromosomal stretches.  Each entry is a tuple of (chromosome identifer, start position, end position, number of contaminant genes in the stretch, number of genes in the stretch, fraction of contaminant genes in a the stretch)
#gene_to_remove: a list of gene identifiers and chromsomal position corresponding to the chromosomal stretches as a tuple with (gene identifier, chromosome_id, start position, end position)
#protein_to_remove: a list of protein identifier to be removed
#gene_to_prot : a dictionary linking all gene identifiers (key) to all of their protein identifiers (values as a list)
#outfile: path to the output report file
def write_report(contaminant_stretches, gene_to_remove, protein_to_remove, gene_to_prot, outfile):
    all_str = []
    all_str.append('-- Contaminant contig segments --')
    all_str.append('\t'.join(['Chromosome', 'Start', 'End', 'Number of "contaminants" genes', 'Number of genes']))
    for cont in contaminant_stretches:
        all_str.append('\t'.join([str(x) for x in cont]))
    all_str.append('-- Contaminant genes --')
    all_str.append('\t'.join(['Gene','Chromosome', 'Start', 'End']))
    for gene in gene_to_remove:
        all_str.append('\t'.join([gene[3]]+[str(x) for x in gene[0:3]]))
    all_str.append('-- Contaminant proteins -- ')
    all_str.append('\t'.join(['Protein', 'Gene']))
    prot_to_gene = {}
    for k, v in gene_to_prot.items():
        for prot in v:
             prot_to_gene[prot] = k
    for prot in protein_to_remove:
        all_str.append('\t'.join([prot, prot_to_gene[prot]]))
    with open(outfile, 'w') as of:
        of.write('\n'.join(all_str))



if __name__=='__main__':


    parser = build_arg_parser()  
    arg = parser.parse_args()
    omark_folder = arg.omark_folder
    gff_file  = arg.gff
    original_fasta = arg.fasta
    output = arg.output

    threshold = arg.threshold
    min_num = arg.minimal_gene_nr
    only_extremities=  arg.only_extrem

    outfasta = arg.output+".fa"
    outtext = arg.output+".txt"


    #Get contaminant proteins in OMARK_FOLDER
    contaminant = get_contaminants(omark_folder)
    #Get the positions of the contaminants and genes from the GFF file
    cont_pos, all_pos, gene_to_prot = get_position_conta(contaminant, gff_file)
    #Define contaminant_stretches, with selected parameter
    contaminant_stretches = infer_contaminant_genome_stretches(cont_pos, all_pos,threshold,min_num,only_extremities)
    #List of protein and genes, from contaminant_stretches
    protein_to_remove, gene_to_remove =  get_genes_in_cont_stretches(contaminant_stretches, all_pos, gene_to_prot)
    #Filter protein from a FASTA and create a filtered copy with contaminant removed
    filter_proteins(original_fasta, outfasta, protein_to_remove)
    #Write a textual report noting removed genes
    write_report(contaminant_stretches,gene_to_remove, protein_to_remove, gene_to_prot, outfile =outtext)