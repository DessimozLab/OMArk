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

from omamer.hierarchy import get_descendants
import ete3
import omark.omamer_utils as outils
from omark.utils import LOG

# Check whether the taxid seem to exist in the NCBI database (Using the get_lineage function)
def check_taxid(taxid):
        try:
            ncbi = ete3.NCBITaxa()
            linid = ncbi.get_lineage(taxid)
        except ValueError:
            LOG.error('The provided taxid is not found in the NCBI database')
            return False
        return True

#Get all lineages from which proteins come in the analyzed proteomes considering the HOGs where the placement was done.
def get_present_lineages(omamdata, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff):
    #This cutoff is made to avoid some false positives. Count lineage only if more than 0.001 of its HOGs are represented
    cutoff_percentage = 0.001
    #Condider only taxa with more than 2 hits
    cutoff_nb_prot = 2

    #Get taxa in which placement were made with the number of placement, couting individual HOGs only once
    all_tax = get_close_taxa_omamer(omamdata, hog_tab, tax_tab, tax_buff, chog_buff,  allow_hog_redun =False)
    
    filter_all_tax = {key: value for (key, value) in all_tax.items() if value > cutoff_nb_prot }

    #Consider only taxa in which at least one percent of the registered HOGs has a hit
    all_taxa_perc = dict()
    hog_by_tax = outils.get_nb_hogs_by_clade(hog_tab, tax_tab)
    proportion_hog_dup = outils.get_prop_duplicated(hog_tab, tax_tab, chog_buff)
    for k, v in filter_all_tax.items():
        all_taxa_perc[k] = float(v)/float(hog_by_tax[k])
    all_taxa_perc = {k: v for k, v in sorted(all_taxa_perc.items(), key=lambda item: item[1], reverse=True)}
    filter_all_taxa_perc = {key: value for (key, value) in all_taxa_perc.items() if value >  cutoff_percentage}

    #Create a tree with all target lineages and uses it to find likely taxa mixture
    t =tree_from_taxlist(filter_all_taxa_perc, tax_tab )
    tax_to_spec = outils.get_spec_by_tax(tax_tab, sp_tab, tax_buff)
    all_plac  = get_likely_spec(t,filter_all_taxa_perc,filter_all_tax, tax_to_spec, proportion_hog_dup)

    return all_plac


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
        taxa = outils.get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
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


def tree_from_taxlist(all_taxa, tax_tab):

    t = ete3.Tree(name='LUCA')
    existing_node = ['LUCA']
    curr_node = t
    #Creating the tree only using lineage present in OMAmer, need to get the full lineage for this
    all_names = all_taxa.keys()
    name_to_lineage = outils.get_full_lineage_omamer(all_names, tax_tab)

    for name, count in all_taxa.items():
        lineage = name_to_lineage[name]
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
def get_likely_spec(t, score, numb, tax_to_spec, prop_inherited):
    cur_score = score.get(t.name.encode(),0)
    cur_numb = numb.get(t.name.encode(),0)
    cur_sp = tax_to_spec[t.name.encode()]
    if t.is_leaf():
        return [(t.name, cur_score, cur_numb,0)]
    all_child = list()
    for child in t.get_children():
        all_child += get_likely_spec(child, score, numb, tax_to_spec, prop_inherited)
    max_score = 0
    best_ranking = None
    qualified = list()
    contaminants =  list()
    new_main = None
    low_depth = 0

    for child in all_child:
        if child[1] >max_score:

            best_ranking = child
            max_score = child[1]
        child_sp = tax_to_spec[child[0].encode()]
        if child[2]>cur_numb*(prop_inherited.get(t.name.encode(), dict()).get(child[0].encode(),0))  and child[1]>cur_score*(len(child_sp)/len(cur_sp))  :
            qualified.append(child)
            if child[3]<=1:
                low_depth += 1

    #If there is only one possibility, and there are no multiple possibilities at a low level directly below
    if len(qualified)==1 or (len(qualified)>1 and low_depth <=1):
        
        for qual in qualified:
            qname = qual[0]
            qscore = qual[1]
            qnumb = qual[2]
            
            qdepth = qual[3]
            if qscore==max_score:
                new_main = (qname, qscore, qnumb,qdepth+1)
            else:
                contaminants.append((qname, qscore, qnumb,qdepth+1))
    if not new_main:
        #There is a case where there is no selected branch, but with some qualified clade. It means the branch
        #with the most representation did not pass the threshold. In these case, we chose the current node: most general.
        new_main = (t.name, cur_score, cur_numb, 0)
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
        taxa = outils.get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
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

    alltaxa = {k: v for k, v in sorted(alltaxa.items(), key=lambda item: item[1], reverse=True)}

    return alltaxa, tax_HOGs


def get_contaminant_proteins(placements, prot_by_clade):
    all_contaminants = list()
    for x in placements[1:]:
        spec = x[0]
        for level  in prot_by_clade[spec]:
            proteins_id = [x[1] for x in level[2]]
            all_contaminants += proteins_id
    return all_contaminants


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

def add_taxid(placements, tax_tab):
    all_sp_name = [x[0] for x in placements]
    #Get the correspondance species name/taxid
    name_to_taxid = outils.get_name_to_taxid(all_sp_name, tax_tab)
    #Add the taxid at the end of each clade descriptor list
    new_placements = [(x[0],x[1], x[2], name_to_taxid[x[0]]) for x in placements]
    return new_placements



#Return the closest ancestor of a clade with more than a threshold of species in omamer
def get_sampled_taxa(clade, threshold_species, tax_tab, sp_tab, tax_buff):
    name_to_lineage = outils.get_full_lineage_omamer([clade], tax_tab)
    lineage = name_to_lineage[clade]
    for tax in lineage:
        species = outils.get_species_from_taxon(tax, tax_tab, sp_tab, tax_buff)
        if len(species)>=threshold_species:
              return tax
    return None


# USER DETERMINED PLACEMENT
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
                    spec = outils.get_species_from_taxon(tax, tax_tab, sp_tab,tax_buff)
                except KeyError:
                    continue
                if len(spec)>=1:
                        return tax.encode('ascii')
        return None
