#!/usr/bin/env python
'''
            OMArk - Quality assesment of coding-gene repertoire annotation
            (C) 2023 Yannis Nevers <yannis.nevers@unil.ch>
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

'''
This scripts is meant to help with data exploration of OMArk results.
It takes as input the path to OMArk results and OMAmer result file.

It is meant to help with annotation quality exploration and offer possibility for correction of the annotation, if possible.
There are multiple modes available for this: a fragment correction mode, a missing genes correction mode and a assembly completeness assesment modes.
Note that this script can help generate file to help with this task but by no mean does it output corrected annotations and will require an expert user to perform this.
A Notebook is available at https://github.com/DessimozLab/OMArk/utils that implement this in an interactive and documented way.
'''

import glob
import pandas as pd
from omadb import Client
import omadb.OMARestAPI
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
import gffutils
from pyfaidx import Fasta
import argparse


def build_arg_parser():
    """Handle the parameter sent when executing the script from
 the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Assists user in extract OMArk data to expore assembly completeness, or identify and correct fragmented or missing genes from OMArk.")
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "assembly completeness" command
    parser_assembly = subparsers.add_parser('assembly', help='Assembly completeness mode - To use with the conserved_hogs.txt from omark -c command. Will download example sequences from these HOGs to obtain fast annotation for completeness assesment of the assembly. Phase 1 corresponds to sequence extraction from OMA, with two modes. Phase 2 to extracting coding sequence from a miniprot output to use as input for OMArk.')
    parser_assembly.set_defaults(func=assembly_mode)
    
    parser_assembly.add_argument('-p', '--phase', help='Phase of the assembly completeness mode. Phase 1 corresponds to sequence extraction from OMA, with two modes. Phase 2 to extracting coding sequence from a miniprot output to use as input for OMArk.', choices=['1','2'], default='1')
    parser_assembly.add_argument('-gf', '--gene_finder', help='Phase 1. Gene finder mode. If activated, the output FASTA file will contain multiple representative by HOG.', action='store_true')
    parser_assembly.add_argument('-o', '--omark_folder', help='Phase 1. The folder output by the omark -c sub-command, containng a conserved_hogs.txt file.')
    parser_assembly.add_argument('-f', '--output_prefix', help='Phase 1 and 2. The path prefix for the output file(s)', default='output_file')
    parser_assembly.add_argument('-m', '--miniprot_output', help='Phase 2. GFF output from Miniprot to use as input of phase 2.')
    parser_assembly.add_argument('-g', '--genomic_fasta', help='Phase 2. The FASTA file of the genome assembly.')

    parser_fragment = subparsers.add_parser('fragment', help='"Fragment mode" of the script. Output sequences for HOGs with fragmented sequences according to OMArk, to identify the full sequence on the assembly if it exists. Output a FASTA file to use with tool such as Miniprot and map to the assembly')
    parser_fragment.set_defaults(func=fragment_mode)
    
    parser_fragment.add_argument('-m', '--omamer_mapping', help='Path to the OMAmer file for the annotation.', type=str, required=True)
    parser_fragment.add_argument('-o', '--omark_folder', help='Path to the OMArk output folder for the annotation.',type=str, required=True )
    parser_fragment.add_argument('-f', '--output_file', help='Path to the output file.', default='fragment_output.fa')

    parser_missing = subparsers.add_parser('missing', help='"Missing mode" of the script. Output sequences for HOGs with missing sequences according to OMArk, to identify the full sequence on the assembly if it exists. Output a FASTA file to use with tool such as Miniprot and map to the assembly. Another option can be used to output the expected genomic context of these genes to look for it on the assembly.')
    parser_missing.set_defaults(func=missing_mode)
    parser_missing.add_argument('-s', '--sequence_or_synteny', help='Different sunmode of the missing mode, as an it. 1.Only sequence 2. Only synteny 3. Sequence and synteny.', choices=['1','2','3'], default='3' )
    parser_missing.add_argument('-m', '--omamer_mapping', help='Path to the OMAmer file for the annotation.', required=True)
    parser_missing.add_argument('-o', '--omark_folder', help='Path to the OMArk output folder for the annotation.', required=True)
    parser_missing.add_argument('-f', '--output_prefix', help='Prefix path to the output file(s).', default='missing_output')

    
    return parser

def extract_completeness_HOGs(omq_file):
    """
        Extracts conserved HOGs data and their respective completeness status in OMArk output

        Args:
            omq_file (str): The path to the input OMArk .omq file

        Returns:
            pandas.DataFrame: A DataFrame with two columns: 'HOG' containing HOgs data and 'Completeness_Category'
            containing completeness status.

        Example:
            omq_file = 'data.txt'
            df = extract_completeness_HOGs(omq_file)
            print(df)

        Output:
                HOG                   Completeness_Category
            0   gene_data_1           section1
            1   gene_data_2           section1
            2   gene_data_3           section2
            ...
        """
    gene_data = []
    section_data = []
    section = None
    with open(omq_file, 'r') as inf:
        for line in inf.readlines():
            line = line.strip('\n')
            if line.startswith('>'):
                section = line[1:]
            else:
                gene_data.append(line)
                section_data.append(section)
    compHOG_df = pd.DataFrame({"HOG": gene_data, 'Completeness_Category' : section_data })
    return compHOG_df

def extract_consistency_genes(ump_file):
    """
    Extracts gene and consistenct data from an input OMArk ump file and returns a pandas DataFrame.

    Args:
        ump_file (str): The path to the input OMArk ump file containing consistency information.

    Returns:
        pandas.DataFrame: A DataFrame with three columns: 'gene' containing gene data, 'Consistency_Category'
        containing taxonomical conssistency status, and 'structure' containing structural consistency status.

    Example:
        ump_file = 'data.txt'
        df = extract_consistency_genes(ump_file)
        print(df)

    Output:
               gene          Consistency_Category    structure
        0   gene_data_1           section1            None
        1   gene_data_2           section1            structure1
        2   gene_data_3           section2            structure2
        ...
    """
    gene_data = []
    section_data = []
    structure_data = []
    section = None
    structure = None
    with open(ump_file, 'r') as inf:
        for line in inf.readlines():
            line = line.strip('\n')
            if line.startswith('>'):
                sect_data = line[1:].split('_')
                if len(sect_data)==2:
                    section = sect_data[0]
                    structure = sect_data[1]
                else:
                    section = line[1:]
                    structure = None
            else:
                gene_data.append(line)
                section_data.append(section)
                structure_data.append(structure)
    comp_gene_df = pd.DataFrame({"gene": gene_data, 'Consistency_Category' : section_data, 'structure': structure_data })
    return comp_gene_df  

def extract_omamer_results(omamer_file):
    """
    Extracts all data from an OMAmer output tsv file.

    Args:
        omamer_file (str): The path to the input omamer.

    Returns:
        pandas.DataFrame: A DataFrame with columns containing the same categories as an omamer file.

    Example:
        omamer_file = 'data.txt'
        df = extract_omamer_results(omamer_file)
        print(df)

    Output:
        ...
    """
    all_data = {}
    gene_data = []
    hog_data = []
    len_data = []
    hog_medlen_data = []
    with open(omamer_file, 'r') as inf:
        firstline = inf.readline()
        while firstline[0]=='!':
            firstline = inf.readline()
        cat_header = firstline.strip('\n').split('\t')
        for line in inf.readlines():
            cat = line.strip('\n').split('\t')
            for index, val in enumerate(cat):
                full_list = all_data.get(cat_header[index], [])
                full_list.append(val)
                all_data[cat_header[index]] = full_list

    gene_hog_df = pd.DataFrame(all_data)
    gene_hog_df= gene_hog_df.rename(columns={'qseqid':'gene', 'hogid':"HOG"})
    gene_hog_df['qseqlen'] = pd.to_numeric(gene_hog_df['qseqlen'], 'coerce')
    if 'subfamily-medianseqlen' in cat_header:
        gene_hog_df['subfamily_medianseqlen'] = pd.to_numeric(gene_hog_df['subfamily-medianseqlen'],'coerce')
    else:
        gene_hog_df['subfamily_medianseqlen'] = pd.to_numeric(gene_hog_df['subfamily_medianseqlen'],'coerce')

    return gene_hog_df

def get_data_total(omark_folder, omamer_file):
    """
    Retrieves and merges data from OMArk and OMAmer results to create a comprehensive DataFrame.

    Args:
        omark_folder (str): The path to the OMARk output folder.
        omamer_file (str): The path to the omamer file.

    Returns:
        tuple: A tuple containing the following DataFrames:
            - full_data: The merged DataFrame containing data from completeness, consistency, and omamer results, merged into a comprehensive dataframe.
            - completeness_data: The DataFrame containing completeness data extracted from OMArk ump file.
            - consistency_data: The DataFrame containing consistency data extracted from OMArk omq file.
            - omamer_data: The DataFrame containing sequence mapping data data extracted from the omamer file.
    """

    hog_file = glob.glob(omark_folder+'/*.omq')[0]
    gene_file = glob.glob(omark_folder+'/*.ump')[0]
    completeness_data = extract_completeness_HOGs(hog_file)
    consistency_data = extract_consistency_genes(gene_file)
    
    omamer_data = extract_omamer_results(omamer_file)
    const_gene_data = consistency_data.merge(omamer_data, on='gene')
    full_data = const_gene_data.merge(completeness_data, on='HOG', how='outer')
    return full_data, completeness_data, consistency_data, omamer_data

def get_level(omark_folder):
    """
    Retrieves the ancestral lineage taxonomic level from a tax file in the specified OMArk folder.

    Args:
        omark_folder (str): The path to the OMArk folder.

    Returns:
        str: The ancestral lineage taxonomic level.

    Example:
        omark_folder = 'data_folder'
        level = get_level(omark_folder)
        print(level)

    Output:
        'Galloanserae'
    """
    tax_file = glob.glob(omark_folder+'/*.tax')[0]
    with open(tax_file) as taxf:
        taxf.readline()
        level = taxf.readline().strip('\n')
    return level

def get_sequences_hog(hog_ids, nseq = 1, medseqlen = None, level=None):
    """
    Retrieves sequences for a list of HOG IDs from the OMA Browser API.

    Args:
        hog_ids (list): A list of HOG IDs for which sequences need to be retrieved.
        nseq (int, optional): The number of sequences to retrieve per HOG ID. Defaults to 1.
        medseqlen (dict, optional): A dictionary mapping HOG IDs to their median sequence lengths. Defaults to None.
        level (str, optional): The taxonomic level at which to retrieve the sequences. Defaults to None.

    Returns:
        dict: A dictionary mapping HOG IDs to their corresponding sequences.

    Example:
        hog_ids = ['HOG:D0621448.1c.5a', 'HOG:D0686017.2a.2b.1b']
        nseq = 2
        medseqlen = {'HOG:D0621448.1c.5a': 100, 'HOG:D0686017.2a.2b.1b': 150}
        level = 'species'
        sequences = get_sequences_hog(hog_ids, nseq, medseqlen, level)
        print(sequences)

    Output:
        {'HOG:D0621448.1c.5a': ['MVEKILDR...', 'MVEKILDK...'], 'HOG:D0686017.2a.2b.1b': ['MADEFRGG...', 'MAEEFRGG...']}
    """
    c = Client()

    seq_hogs= {}
    for hog_id in tqdm(hog_ids):
        try:
            hog = c.hogs[hog_id]
            apir_members = hog['members_url']
            apir_members.uri =  '/'.join(apir_members.uri.split('/')[0:-1])
            if level:
                apir_members.uri += f'?level={level}'
            hog_members = apir_members()['members']
            hog_seq = []
            for entry in hog_members:
                if len(hog_seq)>=nseq :
                    break
                if not medseqlen or entry.sequence_length >= 0.8*medseqlen[hog_id]:
                    hog_entry  =entry['entry_url']()
                    seq =hog_entry['sequence']
                    hog_seq.append(seq)
        except omadb.OMARestAPI.ClientException as ce:
            print(hog_id)
            print(repr(ce))
            hog_seq = None
        seq_hogs[hog_id] = hog_seq    
    return seq_hogs


def write_FASTA_fragmented_HOGs(sequence_of_hog, hog2genes, output_file):
    """
    Writes sequences of HOGs with detected fragments to a FASTA file.

    Args:
        sequence_of_hog (dict): A dictionary mapping HOG IDs to their corresponding sequences.
        output_file (str): The path to the output file where the sequences will be written.

    Returns:
        None

    Example:
        sequence_of_hog = {
            'HOG:D0621448.1c.5a': ['MVEKILDR...', 'MVEKILDK...'], 
            'HOG:D0686017.2a.2b.1b': ['MADEFRGG...', 'MAEEFRGG...']
        }
        output_file = 'output.fasta'
        write_FASTA_fragmented_HOGs(sequence_of_hog, output_file)
    """
    all_rec = []
    for hog, seqs in sequence_of_hog.items():
        if seqs:
            for index, seq in enumerate(seqs):
                identifier = f'{hog}.{index+1}'
                desc = f"({','.join(hog2genes[hog])})"
                record = SeqRecord(
                    Seq(seq),
                    id=identifier,
                    description=desc
                )
                all_rec.append(record)
    with open(output_file,'w') as out:
        SeqIO.write(all_rec, out, "fasta")
        
def write_FASTA_missing_HOGs(sequence_of_hog,  output_file):
    """
    Writes sequences of missing HOGs to a FASTA file.

    Args:
        sequence_of_hog (dict): A dictionary mapping HOG IDs to their corresponding sequences.
        output_file (str): The path to the output file where the sequences will be written.

    Returns:
        None

    Example:
        sequence_of_hog = {
            'HOG:D0621448.1c.5a': ['MVEKILDR...', 'MVEKILDK...'], 
            'HOG:D0686017.2a.2b.1b': ['MADEFRGG...', 'MAEEFRGG...']
        }
        output_file = 'output.fasta'
        write_FASTA_missing_HOGs(sequence_of_hog, output_file)
    """
    all_rec = []
    for hog, seqs in sequence_of_hog.items():
        if seqs:
            for index, seq in enumerate(seqs):
                identifier = f'{hog}.{index+1}'
                record = SeqRecord(
                    Seq(seq),
                    id=identifier,
                )
                all_rec.append(record)
    with open(output_file,'w') as out:
        SeqIO.write(all_rec, out, "fasta")

def omark_input_from_gff(gff_file, fasta, output_fasta, output_splice):
    """
    Reads a GFF file  (for example outputted from miniprot), extracts mRNA features, and generates input files for OMArK pipeline.

    Args:
        gff_file (str): Path to the GFF file containing mRNA features.
        fasta (str): Path to the FASTA file containing the genomic sequences.
        output_fasta (str): Path to the output FASTA file where translated sequences will be written.
        output_splice (str): Path to the output splice mapping file.

    Returns:
        None

    Example:
        gff_file = 'annotations.gff'
        fasta = 'genome.fasta'
        output_fasta = 'translated_sequences.fasta'
        output_splice = 'splice_mapping.txt'
        omark_input_from_gff(gff_file, fasta, output_fasta, output_splice)
    """
    gffdb = gffutils.create_db(gff_file, ':memory:', merge_strategy="create_unique", keep_order = True, checklines=1,force_gff=True)
    seq_records = []
    splice_mapping = []
    fasta = Fasta(fasta)
    have_splice_variants = False
    prev_chr = None
    transcript_of_gene = {'+' : [], '-' : []}
    # Process mRNA features from the GFF file
    for g in gffdb.features_of_type('mRNA', order_by=('seqid','start')):
        new_gene = True
        gene_id = g.id
        chro = g.seqid

        # Check if the chromosome has changed (useful for splice mapping)
        if chro != prev_chr:
            prev_end = {'+': None, '-' : None}
        start_gene = g.start
        end_gene  = g.end
        strand = g.strand
         # Check if the start position is before the previous end position (Transcript overlapping)
        if prev_end[strand]:
            if start_gene < prev_end[strand]:
                new_gene = False
        prev_end[strand] = end_gene
        prev_chr = chro
        seq_combined = ''
        # Find the child features (CDS) and order them by start position                
        for translated_feature in ['CDS', 'cds']:
            ordered_child = list(gffdb.children(g, featuretype=translated_feature, order_by='start'))
            if len(ordered_child)!=0: 
                break
        # Reverse the order of CDS if the strand is negative
        if strand=='-':
            ordered_child.reverse()
        start=True

        # Concatenate the sequences of CDS chunks
        for i in ordered_child:
            if i.start > i.end:
                new_end = i.start
                i.start = i.end
                i.end = new_end
            #The use_strand option make sure the sequence is read in reverse complement on the reverse strand
            if start:
                if i.frame=='1':
                    seq_combined+='NN'
                elif i.frame=='2':
                    seq_combined+='N'
                start=False
            seq = i.sequence(fasta, use_strand=True) 
            seq_combined += seq
        seq_combined = Seq(seq_combined)
        # Translate the sequence using Biopython
        seq_transl = seq_combined.translate()
        # Create a SeqRecord for the translated sequence
        rec = SeqRecord( seq_transl,  id=gene_id, description= '')
        seq_records.append(rec)

        # Track transcript IDs for splice mapping
        if new_gene and transcript_of_gene[strand]!=[]:
            splice_mapping.append(transcript_of_gene[strand])
            transcript_of_gene[strand] = [] 
        transcript_of_gene[strand].append(gene_id)

    # Append remaining transcript IDs to splice mapping
    for strand in ['+', '-']:
        if transcript_of_gene[strand] != []:
            splice_mapping.append(transcript_of_gene[strand])
    # Write translated sequences to FASTA file
    with open(output_fasta, 'w') as out:
        SeqIO.write(seq_records, out, 'fasta')

    # Write splice mapping to file
    with open(output_splice , 'w') as ous:
        for splice_data in splice_mapping:
            ous.write(';'.join(splice_data)+'\n')

def get_synteny_hog(hog_ids, level):
    """
    Retrieve ancestral synteny adjacencies for a given list of HOG IDs at a specified taxonomic level.

    Parameters:
    - hog_ids (list): A list of HOG IDs for which synteny groups will be retrieved.
    - level (int): The taxonomic level level at which to retrieve the ancestral synteny.

    Returns:
    - synteny_groups (dict): A dictionary mapping each HOGs ID to its corresponding syntenic neighbourhood. 
                            The synteny group is represented as a list of HOGs IDs, ordered based on their order on the gemome.
                            If an error occurs during retrieval, the synteny group for that hog ID will be [HOG ID] itself since it means there are no known ancestral neighbour.

    """
    import networkx as nx

    c = Client()
    synt = c.synteny

    synteny_groups = {}
    for hog in tqdm(hog_ids):
        try:
            g = synt.neighbourhood(hog, level=level)
            adj = {x[0] : x[1] for x in g.adjacency()}
            main = adj[hog]
            hright = None
            hleft= None
            rightid = list(main.keys())[0]
            right_neigh = adj[rightid]
            for neighbour_id , neighbour in right_neigh.items():
                if neighbour_id != hog:
                    hright = neighbour_id
            if len(main)>1:
                leftid = list(main.keys())[1]
                left_neigh = adj[leftid]
                for neighbour_id , neighbour in left_neigh.items():
                    if neighbour_id != hog:
                        hleft = neighbour_id
                reordered_adj = [leftid,hog,rightid]
                if hleft:
                    reordered_adj = [hleft] + reordered_adj
                if hright:
                    reordered_adj.append(hright)
            else:
                reordered_adj = [hog,rightid]
                if hright:
                    reordered_adj.append(hright)
            synteny_groups[hog] = reordered_adj
        except omadb.OMARestAPI.ClientException as ce:
            print(hog)
            print(repr(ce))
            synteny_groups[hog] = [hog]
    return synteny_groups

def translate_to_genomic_context(synteny_groups, omamer_map):
    """
    Translate synteny groups to their genomic context in the genome of interest using an OMAmer mapping.

    Parameters:
    - synteny_groups (dict): A dictionary mapping HOG IDs to their corresponding ancestral synteny contexts.
                             The synteny groups should be represented as lists of HOG IDs.
    - omamer_map (dict): A dictionary mapping HOG IDs to genes in the target genome.

    Returns:
    - expected_neighbourhood (list): A sorted list of tuples, where each tuple contains the HOG ID, 
                                     the number of genes retrieved in the context, and the genomic context of the missing gene.
                                     The genomic context is represented as a list of protein IDs, where the target HOG is labeled as 'Target'
                                     and other missing HOGs in the context are labeled as 'Missing'.

    """
    expected_neighbourhood = {}
    for hog, synteny in synteny_groups.items():
        local_context = [x for x in synteny]
        context_score = 0
        for index, syn_hog in enumerate(synteny):
            if syn_hog == hog:
                local_context[index] = 'Target'
            else:
                local_context[index] = omamer_map.get(syn_hog, 'Missing')
                if local_context[index] != 'Missing':
                    context_score += 1
        expected_neighbourhood[hog] = (context_score, local_context)
    expected_neighbourhood = sorted(expected_neighbourhood.items(), key=lambda x:x[1][0], reverse=True)
    return expected_neighbourhood

def read_conserved_hogs(hogfile):
    """
    Read conserved HOG IDs and extract the ancestral lineage level from a conserved_hogs.txt file output by OMArk -c mode.

    Parameters:
    - hogfile (str): The path to the conserved_hogs.txt file output by OMArk -c mode.

    Returns:
    - hog_list (list): A list of conserved HOG IDs extracted from the file.
    - level (str): The ancestral lineage taxonomic level extracted from the file.
    """
    hog_list = []
    with open(hogfile, 'r') as hf:
        for line in hf.readlines():
            line= line.strip('\n')    
            if line.startswith('#'):
                if line[1:].startswith('Ancestral lineage'):
                    level = ':'.join(line.split(':')[1:]).strip(' ')
            else:
                hog_list.append(line)
    return hog_list, level

def write_synteny_file(context, outfile):
    """
    Write synteny context information to a file.

    Parameters:
    - context (list): A list of tuples containing synteny context information.
                      Each tuple should contain the HOG ID, the number of present neighbors, and the expected neighborhood.
                      The expected neighborhood is represented as a hyphen-separated string of hog IDs.
    - outfile (str): The path to the output file.

    Returns:
    - None
    """
    with open(outfile, 'w') as out:
        out.write('\t'.join(['Missing HOG', 'Nr of present neighbour', 'Expected neighbourhood'])+'\n')
        for hog, data in context:
            out.write('\t'.join([hog]+[str(data[0])]+ ['--'.join(data[1])])+'\n')

def assembly_mode(args):
    
    phase = int(args.phase)
    output_file = args.output_prefix
    if phase ==1:
        omark_folder = args.omark_folder
        if omark_folder == None:
            print("Please specify the OMArk folder to read from (-0 [OMArk_folder].")
            exit()
        gfinder = args.gene_finder
        hog_list_file = omark_folder + '/conserved_HOGs.txt'
        hog_list, level = read_conserved_hogs(hog_list_file)
        nseq=1
        if gfinder:
            nseq=20
        sequence_of_hog = get_sequences_hog(hog_list, nseq=nseq, level=level)
        write_FASTA_missing_HOGs(sequence_of_hog, output_file+'.fa')

    elif phase==2:
        miniprot_gff = args.miniprot_output
        genomic_fasta = args.genomic_fasta
        if miniprot_gff==None:
            print('Please specify a Miniprot GFF file (-m [miniprot_gff])')
            exit()
        if genomic_fasta==None:
            print('Please specify a genomic FASTA file (-g [FASTA_file])')
            exit()
        omark_input_from_gff(miniprot_gff, genomic_fasta,output_file+'.fa' , output_file+'.splice')


def fragment_mode(args):

    
    omamer = args.omamer_mapping
    omark_folder = args.omark_folder
    output_file = args.output_file 
    full_df, completeness_df, consistency_df, omamer_df = get_data_total(omark_folder, omamer)
    kept_info = ['gene', 'Consistency_Category', 'structure', 'HOG', 'qseqlen', 'subfamily_medianseqlen', 'Completeness_Category']
    fragment_df =  full_df[full_df['structure']=='Fragment'][kept_info]
    hog_with_fragment = list(fragment_df['HOG'].unique())
    possible_fragments = full_df[full_df['HOG'].isin(hog_with_fragment)][kept_info].sort_values(by='HOG')
    possible_fragments = possible_fragments[possible_fragments['qseqlen']<0.8*possible_fragments['subfamily_medianseqlen']]
    #Extract uniq HOGs
    uniq_HOGs = list(possible_fragments['HOG'].unique())
    hog_to_medseqlen  = {k: v for k, v in zip(possible_fragments['HOG'], possible_fragments['subfamily_medianseqlen'])}
    hog_genes = {}
    for hog, seq in zip(possible_fragments['HOG'],possible_fragments['gene']):
        glist = hog_genes.get(hog, [])
        glist.append(seq)
        hog_genes[hog] =glist
        
    sequence_of_hog = get_sequences_hog(uniq_HOGs, medseqlen=hog_to_medseqlen)
    write_FASTA_fragmented_HOGs(sequence_of_hog, hog_genes, output_file ) 

def missing_mode(args):
    seqsynt_mode = int(args.sequence_or_synteny)
    omamer = args.omamer_mapping
    omark_folder = args.omark_folder
    output_file = args.output_prefix
    full_df, completeness_df, consistency_df, omamer_df = get_data_total(omark_folder, omamer )
    level = get_level(omark_folder)
    missing_df =  full_df[full_df['Completeness_Category']=='Lost']
    uniq_HOGs = list(missing_df['HOG'].unique())
    
    if seqsynt_mode in [1,3]:
        sequence_of_hog = get_sequences_hog(uniq_HOGs, level=level)
        write_FASTA_missing_HOGs(sequence_of_hog, output_file+'.fa')
    if seqsynt_mode in [2,3]:
        synteny_groups = get_synteny_hog(uniq_HOGs, level)
        omamer_map  = {x : y for x,y in zip(list(omamer_df['HOG']),list(omamer_df['gene']) )}
        expected_neighbourhood = translate_to_genomic_context(synteny_groups, omamer_map)
        write_synteny_file(expected_neighbourhood, output_file+'synteny.tsv')


if __name__ == "__main__":
    parser = build_arg_parser()  
    args = parser.parse_args()
    if hasattr(args, "func"):
    
        args.func(args)
    else:
        parser.print_usage()


