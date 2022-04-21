import argparse
import os
import numpy as np
import omamer
import omamer.database
from omamer.hierarchy import get_descendants, get_leaves, get_root_leaf_offsets , get_children
import omamer_species_placement as osp
import files as io
import species_determination as spd
import omamer_utils as utils



def build_arg_parser():
    """Handle the parameter sent when executing the script from the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Compute an OMA quality score from the OMAmer file of a proteome.")   
    parser.add_argument('-f', '--file', help="The OMAmer file to read." )	
    parser.add_argument('-d', '--database', help="The OMAmer database.")
    parser.add_argument('-o', '--outputFolder', help="The folder containing output data the script wilp generate.")
    parser.add_argument('-t', '--taxid', help='Taxonomic identifier', default=None)
    parser.add_argument('-of', '--og_fasta', help='Original FASTA file', default=None)
    return parser

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


#Reorder informations from different treatment into a single cohesive results for stats about the set of protein coding genes
def score_whole_proteome(found_clade, not_in_clade, partials, fragments, not_mapped, contaminants):

    proteome_res = dict()

    proteome_res['Correct'] = found_clade
    confirmed_contaminants = set(not_in_clade).intersection(set(contaminants))
    proteome_res['Contamination'] = confirmed_contaminants
    misplaced = set(not_in_clade) - set(confirmed_contaminants)
    proteome_res['Erroneous'] = misplaced
    proteome_res['Correct_Partial'] = set(found_clade).intersection(set(partials))
    proteome_res['Correct_Fragment'] =set(found_clade).intersection(set(fragments))
    proteome_res['Erroneous_Partial'] = set(misplaced).intersection(set(partials))
    proteome_res['Erroneous_Fragment'] =set(misplaced).intersection(set(fragments))
    proteome_res['Contamination_Partial'] = set(confirmed_contaminants).intersection(set(partials))
    proteome_res['Contamination_Fragment'] =set(confirmed_contaminants).intersection(set(fragments))
    proteome_res['Not_Placed'] = not_mapped
    return proteome_res

#Deprecated
def print_results(res):
    print(len(res['Found']))
    print(len(res['Lost']))
    print(len(res['Overspecific']))
    print(len(res['Underspecific']))
    print(len(res['Duplicated']))

#Useful if we allow paralogs
#root_hog = list(filter(lambda x : x['ParentOff'] == -1, tabi))
#print(root_hog)

#Mutliple level of a same HOG can be counted as is   
def get_conserved_hogs(clade, hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,   cprot_buff, chog_buff, tax_buff, hogtax_buff,  duplicate, threshold=0.9 ) :
    found_hog = list()
    poss_hog = list()
    seen_hog = list()
    other_cl_hog = list()
    lineage = utils.get_full_lineage_omamer(clade.encode('ascii'), tax_tab, tax_buff, True)
    sp_target = utils.get_species_from_taxon(clade, tax_tab, sp_tab, tax_buff)

    for f in fam_tab:

        hog_off = f['HOGoff']
        hog_num = f['HOGnum']
        hogs = hog_tab[hog_off : hog_off+hog_num]
        for x in hogs:
        #for x in range(hog_off, hog_off+hog_num):
                #tax_ind = hogtax_tab[x: x+2]
                tax_ind = x['HOGtaxaOff']
                tax_num = x['HOGtaxaNum']
                tax_off = hogtax_buff[tax_ind: tax_ind+tax_num]
                tax_name = tax_tab[tax_off]['ID']
                if clade.encode('ascii') in tax_name :
                       poss_hog.append(x)
                       seen_hog.append(x['ID'])
    
    for t in poss_hog:
        # Maybe useful if we allow paralogs
        #all_desc = get_descendant_HOGs(t, tabi, chog_buff)
        sp_hog=list()        
        clade_name = tax_tab[t['TaxOff']]['ID']
        if clade_name not in lineage :
                continue
        prot_num = t["ChildrenProtNum"]
        prot_off = t["ChildrenProtOff"]
	
        if duplicate:
                all_desc = utils.get_descendant_HOGs(t, hog_tab, chog_buff)
                for desc in all_desc:
                        desc_tax_name = tax_tab[desc['TaxOff']]["ID"]
                        if desc_tax_name not in lineage :
                               continue
                        sp_hog += [x[0].decode() for x in utils.get_species_from_omamer(desc,prot_tab, sp_tab, cprot_buff)]
        sp_hog += [x[0].decode() for x in utils.get_species_from_omamer(t,prot_tab, sp_tab, cprot_buff)]
        inter = set(sp_hog).intersection(set(sp_target))
        #print(len(inter))
        #print(len(sp_target))
        if len(inter)>=threshold*len(sp_target):
            found_hog.append(t)
        #else:
        #    print(t)
        #omamer.hierarchy.get_descendant_species_taxoffs(hog_off, tabi, chog_buff, cprot_buff, prot2speoff, speoff2taxoff
    return found_hog, poss_hog



def get_root_HOGs_descendants(lineage, tax_tab, hog_tab, fam_tab, tax_buff):

    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    descendants = tax_tab[omamer.hierarchy.get_descendants(tax2tax_off[lineage], tax_tab, tax_buff)]['ID']
    desc_root_HOGs = list()
    for f in fam_tab:
        off = f['TaxOff']
        sp = tax_tab[off]["ID"]
        if sp in descendants:
            hog = hog_tab[f['HOGoff']]
            desc_root_HOGs.append(hog)
    return desc_root_HOGs
    


def found_with_omamer(omamer_data, conserved_hogs, hog_tab, chog_buff):
    all_subf = list()
    all_prot = list()    
    found = list()
    seen_hog_id = list()
    results = { 'Single':[] , 'Lost' : [], 'Duplicated': [], 'Underspecific':[], 'Overspecific_S': [], 'Overspecific_D': []}
    for data in omamer_data:
        all_prot.append(data['qseqid'])
        all_subf.append(data['hogid'])

    #Checking all HOGs from the conserved HOG list
    for hog in conserved_hogs :
        done = False
        identifier = hog['OmaID'].decode()
        if identifier in all_subf:
            nb_found = all_subf.count(identifier)
            st_ind = 0
            for i in  range(nb_found):
                ind = all_subf.index(identifier, st_ind)
                found.append(all_prot[ind])
                st_ind = ind+1
            if nb_found > 1:
                results['Duplicated'].append(identifier)

            else :
                results['Single'].append(identifier)
            done = True

        count_os = 0
        for subhog in [x['OmaID'].decode() for x in utils.get_descendant_HOGs(hog, hog_tab, chog_buff)]:
		
            if subhog in all_subf:
                if subhog not in seen_hog_id:
                    seen_hog_id.append(subhog)
                    nbf = all_subf.count(subhog)
                    st_ind = 0
                    for i in range(nbf):
                        ind = all_subf.index(subhog, st_ind)
                        found.append(all_prot[ind])
                        st_ind = ind+1
                count_os+=1
        if not done and count_os>0:
            if count_os> 1:
                results['Overspecific_D'].append(identifier)
            else:
                results['Overspecific_S'].append(identifier)
            done = True
            
            
        for superhog in [x['OmaID'].decode() for x in utils.get_ancestral_HOGs(hog, hog_tab, chog_buff)]:
            if superhog in all_subf:
                if superhog not in seen_hog_id:
                    seen_hog_id.append(superhog)
                    nbf = all_subf.count(superhog)
                    st_ind = 0
                    for i in range(nbf):
                        ind = all_subf.index(superhog, st_ind)
                        found.append(all_prot[ind])
                        st_ind = ind+1
                if not done:
                    results['Underspecific'].append(identifier)
                    done = True
                    #break
        if not done:
            results['Lost'].append(identifier)

    not_in_clade = list(set(all_prot).difference(set(found)))    
    return results, found, not_in_clade

def get_omamer_qscore(omamerfile, dbpath, stordir, taxid=None, contamination= True, original_FASTA_file = None, force = True):

    db = omamer.database.Database(dbpath)
    #Variables
    hog_tab = db._hog_tab[:]
    prot_tab = db._prot_tab
    sp_tab = db._sp_tab
    tax_tab = db._tax_tab[:]
    fam_tab = db._fam_tab
    cprot_buff = db._cprot_arr
    tax_buff = db._ctax_arr
    chog_buff = db._chog_arr
    hogtax_buff = db._hog_taxa_buff
	 
    allres = dict()
    #Store the temporary results in a file to avoid recomputing and make it computationally feasible
    if os.path.isfile(omamerfile):

        basefile = '.'.join(omamerfile.split('/')[-1].split('.')[:-1])
        if force or not os.path.isfile(stordir+'/'+basefile+".omq"): 
            #Extract OMAmer data
            omamdata, not_mapped  = io.parseOmamer(omamerfile)
            #Get only full match for placement
            full_match_data, partials, fragments = filter_partial_matches(omamdata)
            #Determine species and contamination
            placements = spd.get_present_lineages(full_match_data, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff)      
            #Get the proteins placed in species correspoding to each placement in a dictionary
            prot_clade = spd.get_prot_by_clades(placements, omamdata, hog_tab, tax_tab, tax_buff, chog_buff)
            contaminant_prots = spd.get_contaminant_proteins(placements, prot_clade)
            #Reorganize the placements to consider the species with most proteins to be the main one. (Needed in edge cases where the proteins of the contaminant
            #is an exhaustive set)
            placements = spd.reorganized_placement(placements, prot_clade)

            #Procedure when the user do not give taxonomu information. Will use the main species from the placement
            if taxid==None:
   
                likely_clade =  placements[0][0].encode()
            #Otherwise find the closest clade in the OMAmer database from the given taxid
            else :
                lin = spd.get_lineage_ncbi(taxid)
                likely_clade = spd.find_taxa_from_ncbi(lin, tax_tab, sp_tab,tax_buff)
            #Get the first parent of the chosen clade with at least 5 species
            closest_corr = spd.get_sampled_taxa(likely_clade, 5 , tax_tab, sp_tab, tax_buff)
            #Conshog : HOG with 80% representative of the target lineage
            #Cladehog : HOG with at least 1 representative of the target lineage present in the common ancestir
            conshog, cladehog = get_conserved_hogs(closest_corr.decode(), hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,  cprot_buff,chog_buff, tax_buff, hogtax_buff, True, threshold=0.8)
            lineage_rhog = get_root_HOGs_descendants(closest_corr, tax_tab, hog_tab, fam_tab,tax_buff)

            cladehog = cladehog + lineage_rhog

            #Two modes? : Normal and listing unexpected protein mapping?

            wholeres, found_clade, nic = found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)

            #wholeres, whfound, nic = found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)
           
            res_completeness, found_cons, nicons = found_with_omamer(omamdata ,conshog, hog_tab, chog_buff)

            res_proteomes = score_whole_proteome(found_clade, nic, partials, fragments, not_mapped, contaminant_prots)


            #Store the taxonomic choices
            io.store_close_level(stordir+'/'+basefile+".tax", {'Sampled': str(closest_corr.decode()),
                                                                                'Closest' : str(likely_clade.decode())})

            #Write optionnal taxa files
            if original_FASTA_file:
                if contamination :
                    io.store_contaminant_FASTA(stordir, basefile, prot_clade, original_FASTA_file)
                io.store_incorrect_map_FASTA(stordir, basefile, not_mapped, nic, original_FASTA_file)
            #Write results files
            io.store_results(stordir+'/'+basefile+".ump", {'Unmapped' : not_mapped, 'UnClade' : nic})
            io.store_results(stordir+'/'+basefile+".omq", res_completeness) 
            io.store_summary(stordir+'/'+basefile+".sum",
                            res_completeness, res_proteomes, closest_corr, placements, prot_clade)



if __name__=='__main__':
	
    print('Setting up')
    parser = build_arg_parser()  
    arg = parser.parse_args()
    omamerfile = arg.file
    print(omamerfile)
    dbpath = arg.database
    outdir = arg.outputFolder
    print(outdir)
    taxid = arg.taxid
    print(taxid)
    original_fasta = arg.og_fasta
    get_omamer_qscore(omamerfile, dbpath, outdir, taxid, original_FASTA_file = original_fasta)
    print('Done')

