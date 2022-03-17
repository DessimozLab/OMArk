import pyoma.browser.db
import argparse
import os
import sys
import time
#import pyoma
#import pyham
import inspect
import urllib 
import sys
import gzip
import Bio
import ete3
import re
import numpy as np
from tqdm import tqdm
import pickle
import pandas
if sys.version_info[0] == 3:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
import omamer
import omamer.database
from omamer.hierarchy import get_descendants, get_leaves, get_root_leaf_offsets , get_children
import omamer_species_placement as osp
#import taxonomic_placement as tp


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


#def getCloseTaxaOMAm():
#alltaxa = dict()
#for omamapping in omamerdats:

# get taxonomic levels
def get_hog_taxa(hog_off, sp_tab, prot_tab, hog_tab, cprot_buff, tax_tab, chog_buff):
    '''
    Compute all HOG taxonomic level induced by child HOGs or member proteins.
    '''
    taxa = set()
    
    # add taxa induced by member proteins
    cprot_taxa = np.unique(sp_tab[prot_tab[get_children(hog_off, hog_tab, cprot_buff)]['SpeOff']]['TaxOff'])
    for tax_off in cprot_taxa:
        taxa.update(get_root_leaf_offsets(tax_off, tax_tab['ParentOff']))
    
    # add taxa induced by child HOGs (thus exluding their own taxon)
    chog_taxa = np.unique(hog_tab[get_children(hog_off, hog_tab, chog_buff)]['TaxOff'])
    for tax_off in chog_taxa:
        taxa.update(get_root_leaf_offsets(tax_off, tax_tab['ParentOff'])[:-1])
    
    # remove taxa older than the HOG root-taxon
    hog_tax_off = hog_tab[hog_off]['TaxOff']
    taxa = taxa.difference(get_root_leaf_offsets(hog_tax_off, tax_tab['ParentOff'])[:-1])
    
    return taxa


def get_hog2taxa(hog_tab, sp_tab, prot_tab, cprot_buff, tax_tab, chog_buff):
    '''
    Precompute compact hog2taxa.
    '''
    buff_off = 0
    hog_taxa_idx = [buff_off]
    hog_taxa_buff = []
    for hog_off in tqdm(range(hog_tab.size)):
        taxa = get_hog_taxa(hog_off, sp_tab, prot_tab, hog_tab, cprot_buff, tax_tab, chog_buff)
        buff_off += len(taxa)
        hog_taxa_idx.append(buff_off)
        hog_taxa_buff.extend(taxa)
    return np.array(hog_taxa_idx, dtype=np.int64), np.array(hog_taxa_buff, dtype=np.int16)

def getCloseTaxa(omamerdata, dbpath, tax=None):
    dbObj = pyoma.browser.db.Database(dbpath)
    alltaxa = dict()
    j=0
    descendant = None
    for omamapping in omamerdata:
        print('Family')
        print('----')
        print(omamapping['subfamily'])
        j+=1

        suby = dbObj.get_subhogs(omamapping['subfamily'])
        #The database do not handle too many request, need to reset the connection
        if(j%1==0):
            dbObj = pyoma.browser.db.Database(dbpath)
        if tax :
            curr_taxlist = list()
            for elem in suby:
                descendant = get_children(elem.level, tax)
                if len(descendant)> len(curr_taxlist):
                    curr_taxlist=descendant
            descendant = curr_taxlist
        seen = list()
        for i in suby:
            if i.hog_id==omamapping['subfamily']:
                print(i.level)
                if i.level not in seen:
                    seen.append(i.level)
                if i.level in alltaxa:
                    alltaxa[i.level]+=1
                else:
                    alltaxa[i.level]=1
                if tax:
                    if i.level not in descendant:
                        print([x.level for x in suby])
                        continue
                    descendant.remove(i.level)
        if tax:
            for loss in descendant:
                if loss in alltaxa:
                    alltaxa[loss] -=1
                else :
                    alltaxa[loss]=-1
    alltaxa = {k: v for k, v in reversed(sorted(alltaxa.items(), key=lambda item: item[1]))}
    return alltaxa


def get_full_lineage_omamer(taxname, tax_tab, tax_buff = False,  descendant = False):
    lineage = list()
    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    reached = False
    #print(tax2tax_off) 
    current_tax = tax_tab[tax2tax_off[taxname]]
    while not reached: 
        lineage.append(current_tax['ID'])
        ancestor_tax = current_tax['ParentOff']
        if ancestor_tax!=-1:
                current_tax  = tax_tab[ancestor_tax]
        else:
                reached = True
    if descendant :
        print(lineage)
        print(tax_tab[get_descendants(tax2tax_off[taxname], tax_tab, tax_buff)]['ID'])
        lineage += tax_tab[get_descendants(tax2tax_off[taxname], tax_tab, tax_buff)]['ID'].tolist()
    return lineage
    
	


def get_close_taxa_omamer(omamerdata, hog_tab, tax_tab, ctax_buff, chog_buff, allow_hog_redun =True):
    
    alltaxa = dict()
    j=0
    descendant = None
    seen_hogs = list()

    hog_off2subf = hog_tab['OmaID'] 
    subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
    tax_off2tax = tax_tab['ID'] 
    
    for omamapping in omamerdata:
        j+=1
        if omamapping['hogid'] == 'na':
            continue	
        hog_off = subf2hog_off[omamapping['hogid'].encode('ascii')]       
        taxa = get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
        if not allow_hog_redun:
            if omamapping['hogid'] in seen_hogs:
                continue
            else:
                seen_hogs.append(omamapping['hogid'])
        for taxon in taxa:
            taxname = tax_off2tax[taxon]
            if taxname in alltaxa :
                alltaxa[taxname]+=1
            else: 
                alltaxa[taxname]=1
    #print(len(alltaxa))
    #print(alltaxa)
    
    alltaxa = {k: v for k, v in sorted(alltaxa.items(), key=lambda item: item[1], reverse=True)}
    #alltaxa = { k:v for k in sorted(alltaxa.iteritems(), key=itemgetter(1), reverse=True)}
    return alltaxa

def get_HOGs_taxa_omamer(omamerdata, hog_tab, tax_tab, ctax_buff, chog_buff, allow_hog_redun =True):
    tax_HOGs = dict()
    alltaxa = dict()
    j=0
    descendant = None
    seen_hogs = list()
 
    hog_off2subf = hog_tab['OmaID'] 
    subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
    tax_off2tax = tax_tab['ID']
    for omamapping in omamerdata:
        j+=1
        if omamapping['hogid'] == 'na':
            continue    
        hog_off = subf2hog_off[omamapping['hogid'].encode('ascii')]       
        taxa = get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
        if not allow_hog_redun:
            if omamapping['hogid'] in seen_hogs:
                continue
            else:
                seen_hogs.append(omamapping['hogid'])
  
        for taxon in taxa:
            taxname = tax_off2tax[taxon]
            if taxname in alltaxa :
                alltaxa[taxname]+=1
                tax_HOGs[taxname].append((omamapping['hogid'],omamapping['qseqid']))
            else: 
                alltaxa[taxname]=1
                tax_HOGs[taxname] = list()
                tax_HOGs[taxname].append((omamapping['hogid'],omamapping['qseqid']))

    #print(len(alltaxa))
    #print(alltaxa)
    
    alltaxa = {k: v for k, v in sorted(alltaxa.items(), key=lambda item: item[1], reverse=True)}

    return alltaxa, tax_HOGs

def get_lineage_comp(alltaxa, clade, tax_tab, tax_buff):
    lineage = get_full_lineage_omamer(clade.encode('ascii'), tax_tab, tax_buff, True)
    compatible = 0
    non_comp = 0
    for k, v in alltaxa.items():
        if k in lineage:
            compatible+=v
        else:
            non_comp+=v
    return compatible, non_comp


def get_lower_noncontradicting(alltaxa, tax_tab):
    current_lower_name = None
    current_lower_lineage = None
    for name, count  in alltaxa.items():
        lineage = get_full_lineage_omamer(name, tax_tab)
        if not current_lower_name:
               current_lower_name = name
               current_lineage = lineage
        elif name in current_lineage:
               pass
        else:
               if current_lower_name in lineage:
                       current_lower_name = name
                       current_lineage = lineage
               else:
                       return current_lower_name
    return current_lower_name

def get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff):
    '''
    implied because include taxa having lost their copy
    '''
    #Get the tax-off of the target HOG
    tax_off = hog_tab[hog_off]['TaxOff']
    #Get all taxa that descend from the taxon of target HOG
    hog_taxa = set()    
    #hog_taxa = set(get_descendant_taxa(tax_off, tax_tab, ctax_buff))
    
    hog_taxa.add(tax_off)
    chogs_taxa = set()
    #Substract the taxa that are in a subhog of this family
    #for chog_off in _children_hog(hog_off, hog_tab, chog_buff):
    #    ctax_off = hog_tab[chog_off]['TaxOff']
    #    chogs_taxa.add(ctax_off)
    #    chogs_taxa.update(get_descendant_taxa(ctax_off, tax_tab, ctax_buff))
    return hog_taxa.difference(chogs_taxa)

def get_lineage_ncbi(taxid):
        lineage = list()
        ncbi = ete3.NCBITaxa()
        linid = ncbi.get_lineage(taxid)
        linmap = ncbi.get_taxid_translator(linid)
        lineage = [linmap[x] for x in linid]
        return lineage

def find_taxa_from_ncbi(lineage, tax_tab, sp_tab, tax_buff):
        spec  = []
        for tax in reversed(lineage):
                try:
                    spec = get_species_from_taxon(tax, tax_tab, sp_tab,tax_buff)
                except KeyError:
               	    continue
                if len(spec)>=1:
                        return tax.encode('ascii')
        return None

def getLineage(taxid, tax_tab, sp_tabm, tax_buff):
	ncbi = ete3.NCBITaxa()
	name = ncbi.get_taxid_translator(taxid)
	sp_tax = get_species_from_taxon(name, tax_tab, sp_tab, tax_buff)
	



def get_species_from_taxid(taxid, tax_tab, sp_tab, tax_bugg):
	i=-1
	for taxi in tax_tab:
		print(tax_tab)

def get_species_from_taxon(taxname, tax_tab, sp_tab, tax_buff):
    tax_off2tax = tax_tab['ID'] 
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    tax_off = tax2tax_off[taxname.encode('ascii')]
    sp_off_in_tax = omamer.hierarchy.get_leaves(tax_off, tax_tab, tax_buff)
    sp_tax =[ tax_tab[x][0].decode() for x in sp_off_in_tax]
    return sp_tax

#Could be more efficient. For now fetch multiple time children of Tax. See if need improvement
#Return a dictionary of set of species by clade.
def get_spec_by_tax(tax_tab, sp_tab, tax_buff):
    i=0
    spec_by_tax = dict()
    for tax_entry in tax_tab:
        taxname = tax_entry['ID']
        sp_off_in_tax = omamer.hierarchy.get_leaves(i, tax_tab, tax_buff)
        sp_tax =[ tax_tab[x][0].decode() for x in sp_off_in_tax]
        spec_by_tax[taxname] = set(sp_tax)
        i+=1
    return spec_by_tax

def get_species_from_omamer(hog, prot_tab, spe_tab, cprot_buff) :
    sp_list = list()
    chog_off = hog["ChildrenProtOff"]
    prots = cprot_buff[chog_off : chog_off + hog["ChildrenProtNum"]]

    for p in prots:        
        spe_off = prot_tab[p][1]
        sp_list.append(spe_tab[spe_off])

    return sp_list

def get_ancestral_HOGs(hog, hog_tab, chog_buff):
    all_hogs = list()
    hog_off = hog['ParentOff'] 
    if hog_off != -1:
        anc_hog = hog_tab[hog_off]
        all_hogs.append(anc_hog)
        all_hogs += get_ancestral_HOGs(anc_hog, hog_tab, chog_buff)
    return all_hogs

def get_descendant_HOGs(hog, hog_tab, chog_buff):
    all_hogs = list()
    hog_off = hog['ChildrenOff'] 
    hog_num = hog['ChildrenNum']
    desc_hog = chog_buff[hog_off : hog_off+hog_num]
    subhogs = [ hog_tab[x] for x in desc_hog]
    all_hogs += subhogs
    for subhog in subhogs :
        all_hogs += get_descendant_HOGs(subhog, hog_tab, chog_buff)
    return all_hogs


def get_nb_hogs_by_clade(hog_tab, tax_tab):
    hog_by_tax = dict()
    for hog_off in range(hog_tab.size):
        taxoff = hog_tab[hog_off]['TaxOff']
        tax = tax_tab[taxoff]["ID"]
        if tax not in hog_by_tax:
            hog_by_tax[tax]=0
        hog_by_tax[tax]+=1
    return hog_by_tax

#Species placements

#Get all lineages from which proteins come in the analyzed proteomes considering the HOGs where the placement was done.
def get_present_lineages(omamdata, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff):
    cutoff_percentage = 0.01
    #Get taxa in which placement were made with the number of placement, couting individual HOGs only once
    all_tax = get_close_taxa_omamer(omamdata, hog_tab, tax_tab, tax_buff, chog_buff,  allow_hog_redun =False)
    
    #Condider only taxa with more than 2 hits
    filter_all_tax = {key: value for (key, value) in all_tax.items() if value > 1 }

    #Consider only taxa in which at least one percent of the registered HOGs has a hit
    all_taxa_perc = dict()
    hog_by_tax = get_nb_hogs_by_clade(hog_tab, tax_tab)
    for k, v in filter_all_tax.items():
        all_taxa_perc[k] = float(v)/float(hog_by_tax[k])
    all_taxa_perc = {k: v for k, v in sorted(all_taxa_perc.items(), key=lambda item: item[1], reverse=True)}
    filter_all_taxa_perc = {key: value for (key, value) in all_taxa_perc.items() if value > cutoff_percentage }

    #Create a tree with all target lineages and uses it to find likely taxa mixture
    t =tree_from_taxlist(filter_all_taxa_perc, tax_tab )
    tax_to_spec = get_spec_by_tax(tax_tab, sp_tab, tax_buff)
    all_plac  = get_likely_spec(t,filter_all_taxa_perc,tax_to_spec)

    return all_plac


def tree_from_taxlist(all_taxa, tax_tab):

    t = ete3.Tree(name='LUCA')
    existing_node = ['LUCA']
    curr_node = t
    #Creating the tree only using lineage present in OMAmer, need to get the full lineage for this
    for name, count  in all_taxa.items():
        lineage = get_full_lineage_omamer(name, tax_tab)
        for clade in reversed(lineage):
            clade = clade.decode()
            if clade not in existing_node:
                curr_node = curr_node.add_child(name=clade)
                existing_node.append(clade)
            else:
                curr_node = t&clade

    return t

#Rather than checking depth: check significant linearity to main branch + depth higher node. Maybe try with a new version
#Return a list of tuple: (best ranking clade [allow for selection of lesser clade], maximum score, depth of the selected clade, continuity, depth of best_clade]?
def get_likely_spec(t, score, tax_to_spec):
    cur_score = score.get(t.name.encode(),0)
    cur_sp = tax_to_spec[t.name.encode()]
    if t.is_leaf():
        return [(t.name, cur_score,0)]
    all_child = list()
    for child in t.get_children():
        all_child += get_likely_spec(child,score, tax_to_spec)
    max_score = 0
    best_ranking = None
    qualified = list()
    contaminants =  list()
    new_main = None
    low_depth = True

    for child in all_child:
        if child[1] > max_score:

            best_ranking = child
            max_score = child[1]
        child_sp = tax_to_spec[child[0].encode()]
        if child[1]>cur_score*(len(child_sp)/len(cur_sp)) :

            qualified.append(child)
            if child[2]>1:
                low_depth = False

    
    if len(qualified)==1 or (len(qualified)>1 and not low_depth):
        for qual in qualified:
            name = qual[0]
            score = qual[1]
            depth = qual[2]
            if score==max_score:
                new_main = (name, score,depth+1)
            else:
                contaminants.append((name, score,depth+1))
    if not new_main:
        #There is a case where there is no selected branch, but with some qualified clade. It means the branch
        #with the most representation did not pass the threshold. In these case, we chose the current node: most general.
        new_main = (t.name, cur_score, 0)
        contaminants = list()
    return [new_main] + contaminants

def get_prot_by_clades(all_plac, omamdata, hog_tab, tax_tab, tax_buff, chog_buff):
    all_tax, prot_by_tax = get_HOGs_taxa_omamer(omamdata, hog_tab, tax_tab, tax_buff, chog_buff,  allow_hog_redun = True)
    t_complete = tree_from_taxlist(all_tax, tax_tab)
    prot_clade = compute_protein_breakdowm(all_plac, t_complete, prot_by_tax)
    return prot_clade



def compute_protein_breakdowm(all_plac,t, prot_by_taxa):
    prots_by_clade = dict()
    clade_list = [ x[0] for x in all_plac]
    node_list = [t&x for x in clade_list]
    ancestor_list = list()
    for node in node_list:
        cur_closest_ancestor = None
        cur_depth = 0
        for snd_node in node_list:
            if node!=snd_node:
                    ca = t.get_common_ancestor(node, snd_node)
                    depth = t.get_distance(ca)
                    if not cur_closest_ancestor or depth>cur_depth:
                        cur_closest_ancestor = ca
                        cur_depth = depth
        ancestor_list.append(cur_closest_ancestor)
    for i in range(len(node_list)):
        clade = clade_list[i]
        cur_node = node_list[i]
        ancestor = ancestor_list[i]
        seen_nodes = list()
        level = 0
        all_prots = list()

        while cur_node != ancestor:
            all_prots.append((level, cur_node.name, prot_by_taxa.get(cur_node.name.encode(),[])))
            seen_nodes.append(cur_node.name)
            all_children = cur_node.get_descendants()
            for node_child in all_children:
                if node_child.name not in seen_nodes:
                    all_prots.append((level, node_child.name, prot_by_taxa.get(node_child.name.encode(),[])))
                    seen_nodes.append(node_child.name)

            cur_node = cur_node.up
            level += 1 
        
        prots_by_clade[clade] = all_prots
    return prots_by_clade

def reorganized_placement(placements, prot_by_clade):
    new_placements = list()
    for clade in placements:
        count = 0
        clname = clade[0]
        for levels in prot_by_clade[clname]:
            count += len(levels[2])
        #Drop depth because it does not matter here
        new_placements.append((clname,clade[1],count))
    new_placements.sort(key = lambda x: x[2], reverse=True)
    return new_placements

#Reorder informations from different treatment into a single cohesive results for stats about the set of protein coding genes
def score_whole_proteome(found_clade, not_in_clade, partials, fragments, not_mapped):

    proteome_res = dict()

    proteome_res['Correct'] = found_clade
    proteome_res['Erroneous'] = not_in_clade
    proteome_res['Correct_Partial'] = set(found_clade).intersection(set(partials))
    proteome_res['Correct_Fragment'] =set(found_clade).intersection(set(fragments))
    proteome_res['Erroneous_Partial'] = set(not_in_clade).intersection(set(partials))
    proteome_res['Erroneous_Fragment'] =set(not_in_clade).intersection(set(fragments))
    proteome_res['Not_Placed'] = not_mapped
    return proteome_res

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
    lineage = get_full_lineage_omamer(clade.encode('ascii'), tax_tab, tax_buff, True)
    sp_target = get_species_from_taxon(clade, tax_tab, sp_tab, tax_buff)

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
                all_desc = get_descendant_HOGs(t, hog_tab, chog_buff)
                for desc in all_desc:
                        desc_tax_name = tax_tab[desc['TaxOff']]["ID"]
                        if desc_tax_name not in lineage :
                               continue
                        sp_hog += [x[0].decode() for x in get_species_from_omamer(desc,prot_tab, sp_tab, cprot_buff)]
        sp_hog += [x[0].decode() for x in get_species_from_omamer(t,prot_tab, sp_tab, cprot_buff)]
        inter = set(sp_hog).intersection(set(sp_target))
        #print(len(inter))
        #print(len(sp_target))
        if len(inter)>=threshold*len(sp_target):
            found_hog.append(t)
        #else:
        #    print(t)
        #omamer.hierarchy.get_descendant_species_taxoffs(hog_off, tabi, chog_buff, cprot_buff, prot2speoff, speoff2taxoff
    return found_hog, poss_hog


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
        for subhog in [x['OmaID'].decode() for x in get_descendant_HOGs(hog, hog_tab, chog_buff)]:
		
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
            
            
        for superhog in [x['OmaID'].decode() for x in get_ancestral_HOGs(hog, hog_tab, chog_buff)]:
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

def get_omamer_qscore(omamerfile, dbpath, stordir, taxid=None, unmapped=True, contamination= True, original_FASTA_file = None, force = True):

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
            omamdata, not_mapped  = parseOmamer(omamerfile)
            #Get only full match for placement
            full_match_data, partials, fragments = filter_partial_matches(omamdata)
            #Determine species and contamination
            placements = get_present_lineages(full_match_data, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff)      
            #Get the proteins placed in species correspoding to each placement in a dictionary
            prot_clade = get_prot_by_clades(placements, omamdata, hog_tab, tax_tab, tax_buff, chog_buff)
            #Reorganize the placements to consider the species with most proteins to be the main one. (Needed in edge cases where the proteins of the contaminant
            #is an exhaustive set)
            placements = reorganized_placement(placements, prot_clade)

            #Procedure when the user do not give taxonomu information. Will use the main species from the placement
            if taxid==None:
   
                likely_clade =  placements[0][0].encode()
            #Otherwise find the closest clade in the OMAmer database from the given taxid
            else :
                lin = get_lineage_ncbi(taxid)
                likely_clade = find_taxa_from_ncbi(lin, tax_tab, sp_tab,tax_buff)
            #Get the first parent of the chosen clade with at least 5 species
            closest_corr = osp.get_sampled_taxa(likely_clade, 5 , tax_tab, sp_tab, tax_buff)
            #Store the taxonomic choices
            store_close_level(stordir+'/'+basefile+".tax", {'Sampled': str(closest_corr.decode()),
                                                                                                        'Closest' : str(likely_clade.decode())})
            #Conshog : HOG with 80% representative of the target lineage
            #Cladehog : HOG with at least 1 representative of the target libeage
            conshog, cladehog = get_conserved_hogs(closest_corr.decode(), hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,  cprot_buff,chog_buff, tax_buff, hogtax_buff, True, threshold=0.8)
            #Two modes? : Normal and listing unexpected protein mapping?
            if unmapped :
                print('HOGs')
                print(len(cladehog))
                wholeres, found_clade, nic = found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)
                print('Unmapped')
                print(len(not_mapped))
                print('Not mapped to clade')				
                print(len(nic))
                store_results(stordir+'/'+basefile+".ump", {'Unmapped' : not_mapped, 'UnClade' : nic})
            #wholeres, whfound, nic = found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)
            if original_FASTA_file:
                if contamination :
                    store_contaminant_FASTA(stordir, basefile, prot_clade, original_FASTA_file)
                store_incorrect_map_FASTA(stordir, basefile, not_mapped, nic, original_FASTA_file)
            res_completeness, found_cons, nicons = found_with_omamer(omamdata ,conshog, hog_tab, chog_buff)

            res_proteomes = score_whole_proteome(found_clade, nic, partials, fragments, not_mapped)

            store_results(stordir+'/'+basefile+".omq", res_completeness) 
            store_summary(stordir+'/'+basefile+".sum",
                            res_completeness, res_proteomes, placements, prot_clade)


def store_results(storfile, results):
	with open(storfile, 'w') as storage:
		for categ, hoglist in results.items():
			storage.write('>'+categ+'\n')
			for elem in hoglist:
				storage.write(elem+'\n')

def store_summary(storfile, results, results_proteomes, contaminant = False, prot_clade = False):
    with open(storfile,'w') as storage:
        storage.write('#The clade used was ' +"\n")
        total = len(results['Single'])+len(results['Duplicated'])+ len(results['Overspecific_S']) + len(results['Overspecific_D'])+ len(results['Underspecific']) + len(results['Lost'])
        storage.write('#Number of conserved HOGs is: '+str(total)+'\n')
        tot_genes = len(results_proteomes['Not_Placed'])+len(results_proteomes['Correct'])+len(results_proteomes['Erroneous'])
        storage.write('#Results on conserved HOGs is:\n')
        storage.write('#S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing\n')
        storage.write(f'S:{len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"])},D:{len(results["Duplicated"])+len(results["Overspecific_D"])}[U:{len(results["Duplicated"])},E:{len(results["Overspecific_D"])}],M:{len(results["Lost"])}\n') 
        storage.write(f'S:{100*(len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"]))/total:4.2f}%,D:{100*(len(results["Duplicated"])+len(results["Overspecific_D"]))/total:4.2f}%[U:{100*len(results["Duplicated"])/total:4.2f}%,E:{100*len(results["Overspecific_D"])/total:4.2f}%],M:{100*len(results["Lost"])/total:4.2f}%\n') 
        storage.write('#On the whole proteome, there is '+str(tot_genes)+' proteins\n')
        storage.write('#Of which:\n')
        storage.write('#C:Placements in correct lineage[P:Partial hits,F:Fragmented],E: Erroneous placement[P:Partial hits,F:Fragmented],N: no mapping \n')
        storage.write(f'C:{len(results_proteomes["Correct"])}[P:{len(results_proteomes["Correct_Partial"])},F:{len(results_proteomes["Correct_Fragment"])}],E:{len(results_proteomes["Erroneous"])}[P:{len(results_proteomes["Erroneous_Partial"])},F:{len(results_proteomes["Erroneous_Fragment"])}],N:{len(results_proteomes["Not_Placed"])}\n') 
        storage.write(f'C:{100*len(results_proteomes["Correct"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Correct_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Correct_Fragment"])/tot_genes:4.2f}%],E:{100*len(results_proteomes["Erroneous"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Erroneous_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Erroneous_Fragment"])/tot_genes:4.2f}%],N:{100*len(results_proteomes["Not_Placed"])/tot_genes:4.2f}%\n') 

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


def store_close(storfile, close):
	with open(storfile, 'w') as castor:
		for taxid, num in close.items():
			castor.write(str(taxid)+'\t'+str(num)+'\n')

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

