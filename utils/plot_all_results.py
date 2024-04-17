#!/usr/bin/env python
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

'''
This scripts is meant to help with data exploration and creation of plot with OMArk results from multiple species.
As input, it needs a folder containing all of OMArk output folders.
It will identify the summarised result file and read (.sum).
By default, the filename prefix will be used as species name in the plot.
You can override this behaviour by providing a mapping file where the file prefix is associated with a file name.
Optionally, you can provide a TaxId and order the data according to their taxonomic classification.
This script is provided with a companion Jupyter Notebook that allows more plotting parameters and data visualisation as a pandas Dataframe.
'''

import pandas as pd
import os
import re
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import argparse 

def build_arg_parser():
    """Handle the parameter sent when executing the script from
 the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Generate a multiple species plot for all OMArk folder in a given path.")
    parser.add_argument('-i', '--input', type=str, help="Folder with OMArk results.", required=True)
    parser.add_argument('-o', '--output', type=str, help="Path to the figure to be outputted (PNG format).", default="./omark_multiple.png")
    parser.add_argument('-m', '--mapping', type=str, help="Path to a mapping file. Mapping file shoud be in the tsv format and have a header that lists 'Filename', 'Species name' and 'TaxId'. Filename should indicate the prefix of the sum file. Refer to the template for example.", default=None)
    parser.add_argument('-t', '--taxonomy', help="Flag whether the data should be ranked by taxonomy order. Only works if the TaxId is provided in the mapping file.", action='store_true')

    return parser

#Parse a summary file from OMArk and return a tuple representing the data. First element of the tuple is a Dictionnary containing all OMArk metrics about the proteome of interest. 
#The second element is a list of dictionnary, with one entry per contaminant. The list is empty if there is no detected contaminant
def parse_sum_file(sumfile):

    filebase = re.sub(r'\.sum$', '', os.path.basename(sumfile))

    main_data = dict()
    with open(sumfile) as omaqsum:
            in_cont = False
            detected_species = list()
            for line in omaqsum.readlines():
                if not in_cont:
                    protein_nr_line = re.search(r"On the whole proteome, there \w+ ([0-9]+) proteins",line)
                    if protein_nr_line:
                        main_data['Protein_number'] =int(protein_nr_line.group(1))
                    #S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing
                    #resultline =  re.search("F:([-0-9.]+)%,D:([-0-9.]+)%,O:([-0-9.]+)%,U:([-0-9.]+)%,L:([-0-9.]+)", line)
                    resultline =  re.search(r"S:([-0-9.]+)%,D:([-0-9.]+)%\[U:([-0-9.]+)%,E:([-0-9.]+)%\],M:([-0-9.]+)", line)

                    if resultline:
                        results = resultline
                        main_data['Filename'] = filebase
                        main_data['Species name'] = filebase

                        main_data['Complete'] = float(resultline.group(1))+float(resultline.group(2))
                        main_data['Single'] = float(resultline.group(1))
                        main_data['Duplicated'] = float(resultline.group(2))
                        main_data['Expected_Duplicated'] = float(resultline.group(3)) 
                        main_data['Unexpected_Duplicated'] = float(resultline.group(4))              
                        main_data['Missing'] = float(resultline.group(5))

                    #C:Placements in correct lineage [P:Partial hits, F:Fragmented], E: Erroneous placement [P:Partial hits, F:Fragmented], N: no mapping 
                    #conservline =  re.search("C:([-0-9.]+)%,L:([-0-9.]+)%,O:([-0-9.]+)%,U:([-0-9.]+)%", line)
                    conservline =  re.search(r"A:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],I:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],C:([-0-9.]+)%\[P:([-0-9.]+)%,F:([-0-9.]+)%\],U:([-0-9.]+)%", line)
                    if conservline:
                        main_data['Consistent'] = float(conservline.group(1))
                        main_data['Consistent_Partially_Mapping'] = float(conservline.group(2))
                        main_data['Consistent_Fragments']  = float(conservline.group(3))
                        main_data['Consistent_Structurally_Consistent'] =  float(conservline.group(1))-float(conservline.group(2))-float(conservline.group(3))
                        main_data['Inconsistent'] = float(conservline.group(4))
                        main_data['Inconsistent_Partially_Mapping'] = float(conservline.group(5))
                        main_data['Inconsistent_Fragments'] = float(conservline.group(6))
                        main_data['Inconsistent_Structurally_Consistent'] = float(conservline.group(4))-float(conservline.group(5))-float(conservline.group(6))
                        main_data['Contaminant']= float(conservline.group(7))
                        main_data['Contaminant_Partially_Mapping']= float(conservline.group(8))
                        main_data['Contaminant_Fragments']= float(conservline.group(9))
                        main_data['Contaminant_Structurally_Consistent'] = float(conservline.group(7))-float(conservline.group(8))-float(conservline.group(9))
                        main_data['Unknown'] = float(conservline.group(10))

                    if line == '#From HOG placement, the detected species are:\n':
                        in_cont = True
                else:                        
                    if line[0]=='#' or line=='\n':
                        continue
                    detected_species.append(line.strip('\n').split('\t'))
            main_data['Contaminant'] = [x[0] for x in detected_species[1:]]
            main_data['Detected_Main_Species'] = detected_species[0][0]

            contaminant_data = list()
            if len(detected_species)>1:
                for sup_detect_species in detected_species[1:]:
                    contaminant_data.append({'Filename': filebase, 'Species_name': filebase, 'Main_Taxon' : detected_species[0][0], 'Contaminant' : sup_detect_species[0],
                                        'Contaminant_Taxid': sup_detect_species[1], 'Number_of_Proteins' : sup_detect_species[2]  })
    return main_data, contaminant_data


#Create a dataframe from a list of OMArk folder. It will look for the summarized result file and output it as two Dataframes: one containing all results and another containing only information about contamination events detected in the dataset
def create_df_from_results(file_list):
    
    all_main_data = list()
    all_cont_data = list()
    for path in file_list:
        file_found = False
        for file in os.listdir(path):
            if file[-3:] == 'sum':
                file = os.path.join(path, file)
                file_found = True
                main_data, cont_data = parse_sum_file(file)
                all_main_data.append(main_data)
                all_cont_data += cont_data
        if not file_found:
            #Warning: no file found for path
            pass
    omark_df = pd.DataFrame(all_main_data, columns=['Protein_number','Filename','Species name', 'Complete', 'Single', 'Duplicated','Expected_Duplicated', 
                                         'Unexpected_Duplicated', 'Missing', 'Consistent', 'Consistent_Partially_Mapping', 'Consistent_Fragments',
                                         'Consistent_Structurally_Consistent', 'Inconsistent', 'Inconsistent_Partially_Mapping', 'Inconsistent_Fragments',
                                         'Inconsistent_Structurally_Consistent', 'Contaminant', 'Contaminant_Partially_Mapping', 'Contaminant_Fragments', 'Contaminant_Structurally_Consistent',
                                         'Unknown'])
    cont_df = pd.DataFrame(all_cont_data, columns=['Filename', 'Species_name', 'Main_Taxon', 'Contaminant', 'Contaminant_Taxid', 'Number_of_Proteins'])
    omark_df.sort_values(by='Filename', inplace=True)

    return omark_df, cont_df

#Update the dataframe with species name and taxonmomic informations provided in the mapping file, and optionally order the data taxonomically
def integrate_external_data(omark_df, cont_df, mapping_file, taxonomy_order=False):
    
    mapping = read_mapping_file(mapping_file)
    omark_df['Species name'] = [mapping.get(x,{}).get('Species name',x) for x in omark_df['Filename']]
    omark_df['Taxid']= [mapping.get(x,{}).get('Taxid',None) for x in omark_df['Filename']]
    omark_df['Taxonomy'] = [mapping.get(x,{}).get('Taxonomy',None) for x in omark_df['Filename']]
    if taxonomy_order:
        omark_df = omark_df.sort_values(by='Taxonomy')
    cont_df['Species name'] = [mapping.get(x,{}).get('Species name',x) for x in cont_df['Filename']]
    
    return omark_df,cont_df

#Extract data from an user provided mapping file.
def read_mapping_file(mapfile):
    mapping = dict()
    with open(mapfile) as f:
        datalabel = f.readline().strip('\n').split('\t')
        for elem in datalabel:
            if elem.lower()=='taxid':
                #Load the taxonomy database only if the taxonomy option is used
                from ete3 import NCBITaxa
                ncbi =  NCBITaxa()
        for line in f.readlines():
            filename = None
            taxid = None
            taxonomy = None
            species_name = None
            cat = line.strip('\n').split('\t')
            for i, elem in enumerate(cat):

                if datalabel[i].lower()=='filename':
                    filename = elem
                elif datalabel[i].lower()=='species name':
                    species_name = elem
                elif datalabel[i].lower()=='taxid':
                    taxid = elem
                    try:
                        taxonomy = ncbi.translate_to_names(ncbi.get_lineage(taxid))
                    except ValueError:
                        taxonomy= None
            mapping[filename] = dict()
            if species_name:
                mapping[filename]["Species name"] = species_name
            if taxid:
                mapping[filename]["Taxid"] = species_name
                mapping[filename]["Taxonomy"] = taxonomy
    return mapping

#Plot OMArk results from a Dataframe of results, using the Species name as label. Width, height, fontsize and presence of xtick labels can be set as argument of the function
def plot_omark_df(omark_df, savefile=None, width=None, height=None, no_labels=False, fontsize=12):

    if width==None:
        width = max(len(omark_df)*0.6,6)
    if height==None:
        height = 10
    fig, axes = plt.subplots(3,1,figsize = (width,height))
    omark_df.plot.bar(y='Protein_number', width=1.0, ax=axes[0], lw=0.2, edgecolor='whitesmoke',)

    omark_df.plot.bar(y=['Single', 'Duplicated','Missing'], label=['Single', 'Duplicated', 'Missing'], color = ['#46bea1ff', '#cfee8eff','#ed1c5aff'],  ylim=(0,100), width=1.0, stacked=True, ax=axes[1], rasterized=True)
    xticks = axes[0].get_xaxis().set_visible(False)

    xticks = axes[1].get_xaxis().set_visible(False)
    omark_df.plot.bar(y=['Consistent_Structurally_Consistent','Consistent_Partially_Mapping', 'Consistent_Fragments', 
                               'Contaminant_Structurally_Consistent', 'Contaminant_Partially_Mapping', 'Contaminant_Fragments',
                               'Inconsistent_Structurally_Consistent', 'Inconsistent_Partially_Mapping', 'Inconsistent_Fragments',
                               'Unknown'], x='Species name',  label=['Consistent', '','','Contamination','','','Inconsistent',
                                '','','Unknown'],color =['#3db7e9ff','#3db7e9ff','#3db7e9ff','#e69f00ff','#e69f00ff','#e69f00ff','#8a5df1ff','#8a5df1ff','#8a5df1ff','#000000ff' ],
                               edgecolor='white', lw=0.2, ylim=(0,100), width=1.0, stacked=True, ax=axes[2], rasterized=True)
    xticks = axes[2].xaxis.get_major_ticks()
    bars = axes[2].patches
    hatches= ['','///','///','','///','///','','///','///','']
    col = ['white', 'darkslategray', 'whitesmoke','white', 'darkslategray', 'whitesmoke','white', 'darkslategray', 'whitesmoke','white']
    for index, bar  in enumerate(bars):
        bar.set_hatch(hatches[index//len(omark_df)])
        bar.set_edgecolor(col[index//len(omark_df)])

    partial_patch = mpatches.Patch(facecolor='white',lw=0.1, hatch='///', edgecolor='darkslategray', label='Partial mapping', alpha=.99)
    fragment_patch = mpatches.Patch(facecolor='darkslategray',lw=0.1, hatch='///', edgecolor='white', label='Fragments', alpha=.99)
    custom_legend = [partial_patch,fragment_patch]
    axes[2].legend(title='Consistency', handles=axes[2].get_legend_handles_labels()[0]+custom_legend,loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize)
    axes[1].legend(title='Completeness',loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize)
    axes[0].legend().remove()

    xticks = axes[2].xaxis.get_major_ticks()

    axes[2].set_xlabel('Species', fontsize=fontsize)
    axes[0].set_ylabel('Protein number', fontsize=fontsize)
    axes[1].set_ylabel('Percentage of conserved HOGs', fontsize=fontsize)
    axes[2].set_ylabel('Percentage of proteome', fontsize=fontsize)
    axes[2].set_ylim(axes[1].get_ylim()[::-1])
    axes[2].tick_params(axis='x',labelsize=fontsize)
    if no_labels:
        ticks = axes[2].get_xaxis().set_visible(False)

    plt.subplots_adjust(wspace=0, hspace=0.0)
    if savefile:
        plt.savefig(savefile, bbox_inches='tight')

if __name__=='__main__':
    parser = build_arg_parser()  
    arg = parser.parse_args()
    input_folder = arg.input
    output_figure = arg.output
    mapping_file = arg.mapping


    folders = [os.path.join(input_folder, x) for x in os.listdir(input_folder) if os.path.isdir(os.path.join(input_folder, x))]
    main_df, cont_df =  create_df_from_results(folders)
    if mapping_file:
        main_df, cont_df = integrate_external_data(main_df, cont_df, mapping_file, taxonomy_order=arg.taxonomy)
    if len(main_df)>0:
        plot_omark_df(main_df, savefile=output_figure)
    else:
        print('No data was found in this folder.')
