'''
            OMArk - Quality assesment of coding-gene repertoire annotation
            (C) 2022 Yannis Nevers <yannis.nevers@unil.ch>
            This file is part of OMArk.
            OMArk is free software: you can redistribute it and/or modify
            it under the terms of the GNU Lesser General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            OMArk is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
            GNU Lesser General Public License for more details.
            You should have received a copy of the GNU Lesser General Public License
            along with OMArk. If not, see <http://www.gnu.org/licenses/>.
        '''
import Bio
import re
import jinja2   
import os
from omark.utils import LOG

#This function is called to make sure the input file correspond to OMArk assumption.
#It return a boolean consisting of whether an errror was detected and print the reason for
#errors
def check_omamerfile(omamerfile):
    expected_headers = ['qseqid','hogid','overlap', 'family-score','subfamily-score','qseqlen','subfamily-medianseqlen']
    try:
        with open(omamerfile,'r') as f:
            firstline = f.readline()
            cat = firstline.strip('\n').split('\t')
            if cat!=expected_headers:
                LOG.error('The input OMAmer file is not well formated. Check that your OMAmer version is >=2.2.')
                return False
    except FileNotFoundError:
        LOG.error('The path to the OMAMer file is not valid.')
        return False
    return True

def check_FASTA(fasta_file):
    try:
        with open(fasta_file, "r") as handle:
            fasta = Bio.SeqIO.parse(handle, "fasta")
            if not any(fasta):  # False when `fasta` is empty, i.e. wasn't a FASTA file
                LOG.error("The FASTA was not in the correct format or was empty.")
                return False
    except FileNotFoundError:
            LOG.error('The path to the FASTA file is not valid.')
            return False
    return True

def check_isoform_file(isoform_file):
    try:
        with open(isoform_file,'r') as f:
            pass
    except FileNotFoundError:
        LOG.error('The path to the isoform file is no valid.')
        return False
    return True

def check_and_create_output_folder(stordir):

    if os.path.isdir(stordir):
        return True
    else:
        try:
            os.mkdir(stordir)
        except FileNotFoundError:
            LOG.error('The path to the output directory is not valid (Its parent directory does not exist).')
            return False
        except PermissionError:
            LOG.error('No permission to write to the output directory path. Please check permissions.')
            return False
    return True

#This function read an OMAmer file (input_) and output two variables:
#alldata -> A list of all OMAmer placement, containing a dictionary correspoding to all of the OMAmer data results
#not_mapped -> A listt of proteins that do not map to any HOGs (no homologs)
def parseOmamer(file):
    alldata = list()
    not_mapped = list()

    with open(file) as f:
        
        firstline = f.readline()
        cat = firstline.strip('\n').split('\t')
        for line in f.readlines():
            data = dict()   
            col = line.strip('\n').split('\t')
            for i in range(len(cat)) :
                data[cat[i]] = col[i]
            if data['hogid']=='na':
                not_mapped.append(col[0])
                continue
            alldata.append(data)
    return alldata, not_mapped 

def parse_isoform_file(file):
    isoform_by_gene = list()
    with open(file) as handle:
        for line in handle.readlines():
            line = line.strip('\n')
            splice = line.split(";")
            isoform_by_gene.append(splice)
    return isoform_by_gene


def select_isoform(isoform_data, alldata):
    indexed_data = {x['qseqid'] : x for x in alldata}
    best_scoring_isoforms = list()
    selected_isoforms = list()
    not_mapped_gene = list()
    for gene in isoform_data:
        main_variant = None
        best_score = 0

        for isoform in gene:
            omamer_res = indexed_data.get(isoform)
            if omamer_res:

                score = float(omamer_res['family-score'])*min(int(omamer_res['qseqlen']),int(omamer_res['subfamily-medianseqlen']))
                if score > best_score:
                    best_score = score
                    main_variant = omamer_res
        if main_variant != None:
            best_scoring_isoforms.append(main_variant)
            selected_isoforms.append(isoform)
        else:
            #If there is no main variant, we have no way to select the best isoform. We select one by default, here the latest appearing in the FASTA file
            not_mapped_gene.append(isoform)
            selected_isoforms.append(isoform)
    return best_scoring_isoforms, not_mapped_gene, selected_isoforms


#This function remove sequences that only match partially to HOGs of interests,
#in order to consider only the most robust placement.
#Will remove sequence for which at least 20% of the sequence share no sequence with the HOG
#as well as proteins that are less than half the median size of the sequences in the HOGs
def filter_partial_matches(omamdata, overlap_threshold=0.8, fragment_threshold = 0.5):
    filter_dat = list()
    partials = list()
    fragments = list()
    for data in omamdata:
        if float(data['overlap'])< overlap_threshold:
            partials.append(data['qseqid'])
        elif float(data['qseqlen'])< fragment_threshold*float(data['subfamily-medianseqlen']):
            fragments.append(data['qseqid'])
        else:
            filter_dat.append(data)
    return filter_dat, partials, fragments



def store_results(storfile, results):
    with open(storfile, 'w') as storage:
        for categ, hoglist in results.items():
            storage.write('>'+categ+'\n')
            for elem in hoglist:
                storage.write(elem+'\n')

def store_list(storfile, data):

    with open(storfile, 'w') as storage:
        for elem in data:
            storage.write(elem+'\n')




def write_templated_report(template_file, storfile, results, results_proteomes, selected_lineage, species_report):

    env = jinja2.Environment(
    loader= jinja2.PackageLoader("omark", "assets"),
    autoescape=jinja2.select_autoescape()
    )
    template =  env.get_template(template_file)

    all_stats = organize_results(results, results_proteomes, selected_lineage, species_report)

    with open(storfile, 'w') as outfile:
        outfile.write(template.render(all_stats))


#This function create a detailed dictionnary of all OMArk stats, by post-processing the results of the different analysis.
def organize_results(results, results_proteomes, selected_lineage, species_report):

    ancestral_lineage = selected_lineage.decode()


    single_nr = len(results['Single'])+len(results['Overspecific_S'])+len(results['Underspecific'])
    dup_nr = len(results['Duplicated'])+ len(results['Overspecific_D'])
    dup_exp_nr = len(results['Overspecific_D'])
    dup_unexp_nr = len(results['Duplicated'])
    missing_nr = len(results['Lost'])
    cons_hog_nr = single_nr+dup_nr+ missing_nr
    single_percent = 100*single_nr/cons_hog_nr
    dup_percent = 100*dup_nr/cons_hog_nr
    dup_exp_percent = 100*dup_exp_nr/cons_hog_nr
    dup_unexp_percent = 100*dup_unexp_nr/cons_hog_nr
    missing_percent = 100*missing_nr/cons_hog_nr


    consistent_nr =  len(results_proteomes['Correct'])
    consistent_partial_nr = len(results_proteomes['Correct_Partial'])
    consistent_fragment_nr = len(results_proteomes['Correct_Fragment'])
    inconsistent_nr =  len(results_proteomes['Erroneous'])
    inconsistent_partial_nr = len(results_proteomes['Erroneous_Partial'])
    inconsistent_fragment_nr = len(results_proteomes['Erroneous_Fragment'])
    contamination_nr = len(results_proteomes['Contamination'])
    contamination_partial_nr = len(results_proteomes['Contamination_Partial'])
    contamination_fragment_nr = len(results_proteomes['Contamination_Fragment'])
    no_map_nr = len(results_proteomes['Not_Placed'])
    protein_nr = consistent_nr + inconsistent_nr + contamination_nr + no_map_nr
    
    consistent_percent = 100*consistent_nr/protein_nr
    consistent_partial_percent = 100*consistent_partial_nr/protein_nr
    consistent_fragment_percent = 100*consistent_fragment_nr/protein_nr
    inconsistent_percent = 100*inconsistent_nr/protein_nr
    inconsistent_partial_percent = 100*inconsistent_partial_nr/protein_nr
    inconsistent_fragment_percent = 100*inconsistent_fragment_nr/protein_nr
    contamination_percent = 100*contamination_nr/protein_nr
    contamination_partial_percent = 100*contamination_partial_nr/protein_nr
    contamination_fragment_percent = 100*contamination_fragment_nr/protein_nr
    no_map_percent = 100*no_map_nr/protein_nr

    contaminants = list()
    main = True
    for values in species_report:     
        species_info = {"name" : values[0], "protein_nr" : values[2], "protein_percent" : 100*values[2]/protein_nr, 'taxid': values[3]}
        if main:
            main_clade = species_info
            main = False
        else:
            contaminants.append(species_info)
    all_stats = {  "ancestral_lineage" : ancestral_lineage,
                        "cons_hog_nr" : cons_hog_nr,
                        "single_nr" : single_nr, 
                        "dup_nr" : dup_nr,
                        "dup_exp_nr" : dup_exp_nr,
                        "dup_unexp_nr" : dup_unexp_nr,
                        "missing_nr" : missing_nr,
                        "single_percent" : single_percent,
                        "dup_percent": dup_percent,
                        "dup_exp_percent" : dup_exp_percent,
                        "dup_unexp_percent" : dup_unexp_percent,
                        "missing_percent" : missing_percent,
                        "protein_nr" : protein_nr,
                        "consistent_nr" : consistent_nr,
                        "consistent_partial_nr" : consistent_partial_nr,
                        "consistent_fragment_nr" : consistent_fragment_nr,
                        "inconsistent_nr" : inconsistent_nr,
                        "inconsistent_partial_nr" : inconsistent_partial_nr,
                        "inconsistent_fragment_nr" : inconsistent_fragment_nr,
                        "contamination_nr" : contamination_nr,
                        "contamination_partial_nr" : contamination_partial_nr,
                        "contamination_fragment_nr" : contamination_fragment_nr,
                        "no_map_nr" : no_map_nr,
                        "consistent_percent" : consistent_percent,
                        "consistent_partial_percent" : consistent_partial_percent,
                        "consistent_fragment_percent" : consistent_fragment_percent,
                        "inconsistent_percent" : inconsistent_percent,
                        "inconsistent_partial_percent" : inconsistent_partial_percent,
                        "inconsistent_fragment_percent" : inconsistent_fragment_percent,
                        "contamination_percent" : contamination_percent,
                        "contamination_partial_percent"  : contamination_partial_percent,
                        "contamination_fragment_percent" : contamination_fragment_percent,
                        "no_map_percent" : no_map_percent,
                        "main_clade" : main_clade,
                        "contaminants" : contaminants }
    return all_stats




#Function to write the OMArk .sum file. Deprecated and replaced by the more gene 'write_templated_report'
def store_summary(storfile, results, results_proteomes, selected_lineage, contaminant = False, prot_clade = False):
    with open(storfile,'w') as storage:
        storage.write('#The selected clade was '+selected_lineage.decode()+"\n")
        total = len(results['Single'])+len(results['Duplicated'])+ len(results['Overspecific_S']) + len(results['Overspecific_D'])+ len(results['Underspecific']) + len(results['Lost'])
        storage.write('#Number of conserved HOGs is: '+str(total)+'\n')
        tot_genes = len(results_proteomes['Not_Placed'])+len(results_proteomes['Contamination'])+len(results_proteomes['Correct'])+len(results_proteomes['Erroneous'])
        storage.write('#Results on conserved HOGs is:\n')
        storage.write('#S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing\n')
        storage.write(f'S:{len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"])},D:{len(results["Duplicated"])+len(results["Overspecific_D"])}[U:{len(results["Duplicated"])},E:{len(results["Overspecific_D"])}],M:{len(results["Lost"])}\n') 
        storage.write(f'S:{100*(len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"]))/total:4.2f}%,D:{100*(len(results["Duplicated"])+len(results["Overspecific_D"]))/total:4.2f}%[U:{100*len(results["Duplicated"])/total:4.2f}%,E:{100*len(results["Overspecific_D"])/total:4.2f}%],M:{100*len(results["Lost"])/total:4.2f}%\n') 
        storage.write('#On the whole proteome, there is '+str(tot_genes)+' proteins\n')
        storage.write('#Of which:\n')
        storage.write('#A:Placements in accurate lineage[P:Partial hits,F:Fragmented], E: Erroneous placements[P:Partial hits,F:Fragmented], C: Likely contamination[P:Partial hits,F:Fragmented], N: no mapping \n')
        storage.write(f'A:{len(results_proteomes["Correct"])}[P:{len(results_proteomes["Correct_Partial"])},F:{len(results_proteomes["Correct_Fragment"])}],E:{len(results_proteomes["Erroneous"])}[P:{len(results_proteomes["Erroneous_Partial"])},F:{len(results_proteomes["Erroneous_Fragment"])}],C:{len(results_proteomes["Contamination"])}[P:{len(results_proteomes["Contamination_Partial"])},F:{len(results_proteomes["Contamination_Fragment"])}],N:{len(results_proteomes["Not_Placed"])}\n')
        storage.write(f'A:{100*len(results_proteomes["Correct"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Correct_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Correct_Fragment"])/tot_genes:4.2f}%],E:{100*len(results_proteomes["Erroneous"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Erroneous_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Erroneous_Fragment"])/tot_genes:4.2f}%],C:{100*len(results_proteomes["Contamination"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Contamination_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Contamination_Fragment"])/tot_genes:4.2f}%],N:{100*len(results_proteomes["Not_Placed"])/tot_genes:4.2f}%\n')
        #storage.write(f'C:{100*len(found_cons)/tot_genes:4.2f}%,L:{100*(len(nicons)-len(nic))/tot_genes:4.2f}%,O:{100*len(nic)/tot_genes:4.2f}%,U:{100*len(unmap)/tot_genes:4.2f}%\n')
        if contaminant:
            storage.write('#From HOG placement, the detected species are:\n')
            storage.write("#Clade\tPercentage of clade's HOGs\tNumber of associated proteins\n")

            count=0
            for values in contaminant:
                if count==1:
                    storage.write('#Including possible contaminant:\n')
                storage.write('\t'.join([str(x) for x in values])+'\n')
                count+=1
                #storage.write('\t'+str(len(prot_clade[values[0]][0][2]))+'\n')


def store_contaminant_FASTA(stordir, basefile_name, prot_clade, original_FASTA_file):
    seqs_by_id = dict()
    with open(original_FASTA_file) as handle:
        for record in Bio.SeqIO.parse(handle, "fasta"):
            seqs_by_id[record.id] = record
    for key, value in prot_clade.items():
        seqs_from_cont = list()
        for level_data in value:
            level = level_data[0]
            clade = level_data[1]
            for prot_data in level_data[2]:
                    seq = seqs_by_id[prot_data[1]]
                    seq.description = seq.description +" Level="+str(level)+" ["+clade+"]"
                    seqs_from_cont.append(seq)
            with open(stordir+"/"+basefile_name+"_"+re.sub("[^0-9a-zA-Z]+", "_",key)+".fasta", "w") as out_handle:
                    Bio.SeqIO.write(seqs_from_cont, out_handle, 'fasta')           

def store_incorrect_map_FASTA(stordir, basefile_name, not_mapped, incorrect_plac, original_FASTA_file):
    seqs_by_id = dict()
    with open(original_FASTA_file) as handle:
        for record in Bio.SeqIO.parse(handle, "fasta"):
            seqs_by_id[record.id] = record

    seqs_not_map = list()
    seqs_mapped = list()
    seqs_misplaced = list()
    for seqid, seq in seqs_by_id.items():
        if seqid in not_mapped:
            seqs_not_map.append(seq)
        elif seqid in incorrect_plac:
            seqs_misplaced.append(seq)
        else: 
            seqs_mapped.append(seq)
  
    with open(stordir+"/"+basefile_name+"_mapped.fasta", "w") as out_handle:
        Bio.SeqIO.write(seqs_mapped, out_handle, 'fasta')    
    with open(stordir+"/"+basefile_name+"_no_hits.fasta", "w") as out_handle:
        Bio.SeqIO.write(seqs_not_map, out_handle, 'fasta')
    with open(stordir+"/"+basefile_name+"_misplaced.fasta", "w") as out_handle:
                Bio.SeqIO.write(seqs_misplaced, out_handle, 'fasta')   


def store_close_level(storfile, data):
        with open(storfile ,'w') as castor:
                castor.write('>Sampled\n')
                castor.write(data['Sampled']+"\n")
                castor.write('>Closest\n')
                castor.write(data['Closest']+"\n")
                if 'All' in data:
                      castor.write('>All'+'\n')
                      for taxid, num in data['All'].items():
                            castor.write(str(taxid)+'\t'+str(num)+'\n')
