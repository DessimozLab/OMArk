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
import os

ETE_TAXA_PATH=None

#Change the global variable setting the path to the NCBI Taxonomy database. To use when users want to use a local copy
def set_ete_taxa_path(path):
    global ETE_TAXA_PATH
    ETE_TAXA_PATH = path

#Initialize the ete3 database, using a local path if it was defined before.
def get_ete_ncbi_db():
    if ETE_TAXA_PATH is None:
        ncbi = ete3.NCBITaxa()
    else:
        ncbi = ete3.NCBITaxa(dbfile=ETE_TAXA_PATH)
    return ncbi

# Check whether the taxid seem to exist in the NCBI database (Using the get_lineage function)
def check_taxid(taxid):
        try:
            ncbi = get_ete_ncbi_db()
            linid = ncbi.get_lineage(taxid)
        except ValueError:
            LOG.error('The provided taxid is not found in the NCBI database')
            return False
        return True

def check_rank(taxonomic_rank):
    #Rank options above genus - copied from https://github.com/ropensci/taxize/issues/835
    rank_options = [ 'domain', 'superkingdom', 'kingdom', 'subkingdom', 'infrakingdom', 'superphylum', 'phylum' ,'division', 'subphylum','subdivision', 'infradivision', 'superclass',
                    'class', 'subclass', 'infraclass', 'subterclass', 'parvclass', 'megacohort', 'supercohort', 'cohort', 'subcohort', 'infracohort' , 'superorder',
                    'order', 'suborder', 'infraorder', 'parvorder' , 'superfamily', 'family', 'subfamily', 'supertribe', 'tribe', 'subtribe']
    if taxonomic_rank not in rank_options:
        available_options = ','.join(rank_options)
        LOG.error(f'The provided taxonomic rank is not a valid options. Valid options are {available_options}')
        return False
    return True

#Get all lineages from which proteins come in the analyzed proteomes considering the HOGs where the placement was done.
def get_present_lineages(omamdata, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff, family_score_filter=70, cutoff_percentage=0.000):
    #This cutoff is made to avoid some false positives. Count lineage only if more than 0.001 of its HOGs are represented
    cutoff_percentage = cutoff_percentage
    #Condider only taxa with more than 2 hits
    cutoff_nb_prot = 2
    if family_score_filter:
        omamdata= [x for x in omamdata if float(x['family-score'])>family_score_filter]

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
    if len(filter_all_taxa_perc)==0:
        LOG.warning('Not enough omamer placements for automatic lineage detection. Please check your file is complete.')
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

    if 'hoglevel' not in omamerdata[0]:
        hog_off2subf = hog_tab['OmaID']
        subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
        tax_off2tax = tax_tab['ID']
    
    for omamapping in omamerdata:
        j+=1
        if omamapping['hogid'] == 'na':
            continue
        if not allow_hog_redun:
            if omamapping['hogid'] in seen_hogs:
                continue
            else:
                seen_hogs.append(omamapping['hogid'])
        if 'hoglevel' not in omamapping:
            hog_off = subf2hog_off[omamapping['hogid'].encode('ascii')]       
            taxon = outils.get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
            taxname = tax_off2tax[taxon].decode()
        else:
            taxname = omamapping['hoglevel']
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

    #Creating the tree only using lineage present in OMAmer, need to get the full lineage for this
    all_names = all_taxa.keys()
    if len(all_names)>0:
        name_to_lineage = outils.get_full_lineage_omamer(all_names, tax_tab)
        oldest_ancestor = name_to_lineage[list(all_names)[0]][-1]
    else:
        oldest_ancestor = outils.get_root_clade(tax_tab)
    t = ete3.Tree(name=oldest_ancestor)
    existing_node = [oldest_ancestor]
    curr_node = t
    for name, count in all_taxa.items():
        lineage = name_to_lineage[name]
        for clade in reversed(lineage):
            clade = clade
            if clade not in existing_node:
                curr_node = curr_node.add_child(name=clade)
                existing_node.append(clade)
            else:
                curr_node = t&clade

    return t


#Rather than checking depth: check significant linearity to main branch + depth higher node. Maybe try with a new version
#Return a list of tuple: (best ranking clade [allow for selection of lesser clade], maximum score, depth of the selected clade, continuity, depth of best_clade]?
def get_likely_spec(t, score, numb, tax_to_spec, prop_inherited):
    cur_score = score.get(t.name,0)
    cur_numb = numb.get(t.name,0)
    cur_sp = tax_to_spec[t.name]
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
        child_sp = tax_to_spec[child[0]]
        if child[2]>cur_numb*(prop_inherited.get(t.name, dict()).get(child[0],0))  and child[1]>cur_score*(len(child_sp)/len(cur_sp))  :
            qualified.append(child)
            if child[3]<=1:
                low_depth += 1


    #If there is only one possibility, and there are no multiple possibilities at a low level directly below
    if len(qualified)==1 or (len(qualified)>1 and ( (not ( low_depth>0 and len(cur_sp)<20)) and (not (low_depth==len(qualified))))):
        pass
    #If we have multiple possibilities, will still keep those with an higher score
    else:
        saved_qualified = list()
        for q in qualified:
            if q[1]>cur_score:
                saved_qualified.append(q)
        qualified = saved_qualified
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
    prot_clade = compute_protein_breakdown(all_plac, t_complete, prot_by_tax,True)
    return prot_clade



def compute_protein_breakdown(all_plac,t, prot_by_taxa, include_uncertain=True):
    #Attribute proteins to the species detected from placemebt
    prots_by_clade = dict()
    #Lists of clades detected from placement and position in the tree
    clade_list = [ x[0] for x in all_plac]
    node_list = [t&x for x in clade_list]
    ancestor_list = list()
    for node in node_list:
        all_closest_ancestor = []
        for snd_node in node_list:
            ca = t.get_common_ancestor(node, snd_node)
            depth = t.get_distance(ca)
            all_closest_ancestor.append(ca)
        #Get a list of closest common ancestor to the closest detected species
        ancestor_list.append(all_closest_ancestor)


    seen_nodes = list()
    for i in range(len(node_list)):
        clade = clade_list[i]
        cur_node = node_list[i]
        f_node = node_list[i]
        ancestor = ancestor_list[i]
        level = 0
        all_prots = list()
        included_species = [clade]
        seen_all = False
        while not seen_all:
            unexplored_contaminant = False
            if cur_node in ancestor and cur_node!=f_node:
                for index, elem in enumerate(ancestor):
                    if elem==cur_node:
                        included_species.append(clade_list[index])
                        #If all other descendant have not yet be seen, going to the next species
                        if clade_list[index] not in seen_nodes:
                            unexplored_contaminant = True
            if unexplored_contaminant:
                break
            target_clade = included_species[0] if len(included_species)==1 else tuple(sorted(included_species))

            if (type(target_clade)==tuple and not include_uncertain) or unexplored_contaminant:
                break
            all_prots = prots_by_clade.get(target_clade, [])
            if cur_node.name not in seen_nodes:
                all_prots.append((level, cur_node.name, prot_by_taxa.get(cur_node.name,[])))
                seen_nodes.append(cur_node.name)
            all_children = cur_node.get_descendants()
            for node_child in all_children :
                if node_child.name not in seen_nodes:
                    all_prots.append((level, node_child.name, prot_by_taxa.get(node_child.name,[])))
                    seen_nodes.append(node_child.name)

            level += 1
            prots_by_clade[target_clade] = all_prots
            if cur_node.is_root():
                seen_all = True
            else:
                cur_node = cur_node.up
    return prots_by_clade

def get_HOGs_taxa_omamer(omamerdata, hog_tab, tax_tab, ctax_buff, chog_buff, allow_hog_redun =True):
    tax_HOGs = dict()
    alltaxa = dict()
    j=0
    descendant = None
    seen_hogs = list()
    if not 'hoglevel' in omamerdata[0]:
        hog_off2subf = hog_tab['OmaID']
        subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
        tax_off2tax = tax_tab['ID']
    for omamapping in omamerdata:
        j+=1
        if omamapping['hogid'] == 'na':
            continue    
        if not allow_hog_redun:
                if omamapping['hogid'] in seen_hogs:
                    continue
                else:
                    seen_hogs.append(omamapping['hogid'])
        if not 'hoglevel' in omamapping:
            hog_off = subf2hog_off[omamapping['hogid'].encode('ascii')]
            taxon = outils.get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
            taxname = tax_off2tax[taxon].decode()
        else:
            taxname = omamapping['hoglevel']
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
    main_spec = placements[0][0]
    contaminant_spec = [x[0] for x in placements[1:]]
    for spec_list, levels in prot_by_clade.items():
        if not (spec_list==main_spec or (type(spec_list)==tuple and main_spec in spec_list )):
            for level in levels:
                proteins_id = [x[1] for x in level[2]]                
                all_contaminants += proteins_id
    return all_contaminants

def add_uncertain_contaminants(placements, prot_by_clade):
    all_contaminants = list()
    main_spec = placements[0][0]
    contaminant_spec = [x[0] for x in placements[1:]]
    count_proteins = 0
    new_prot_by_clade = {}
    for spec_list, levels in prot_by_clade.items():
        if type(spec_list)==tuple and main_spec not in spec_list:
            for level in levels:
                proteins_id = [x[1] for x in level[2]]           
                count_proteins += len(proteins_id)
            cur_levels = new_prot_by_clade.get('Ambiguous_contaminant',[])
            new_prot_by_clade['Ambiguous_contaminant'] = cur_levels+levels
        elif type(spec_list)!=tuple:
            new_prot_by_clade[spec_list] = levels
    if count_proteins!=0:
        placements.append(('Ambiguous contaminant', 0, count_proteins, -1))
    
    return placements, new_prot_by_clade

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

def reorganize_placements_from_taxid(placements, likely_clade,tax_tab, tax_buff):
    max_overlap = 0
    best_fit = list()
    all_clades = [likely_clade] + [x[0] for x in placements]
    lineages = outils.get_full_lineage_omamer(all_clades, tax_tab, tax_buff, False)

    main_lineage = set(lineages[likely_clade])
    for candidate in placements:
        candidate_lineage = set(lineages[candidate[0]])
        overlap = len(main_lineage.intersection(candidate_lineage))
        if overlap>max_overlap:
            max_overlap = overlap
            best_fit = candidate
    placements = [best_fit] + [x for x in placements if x!=best_fit]
    return placements

def add_taxid(placements, tax_tab):
    all_sp_name = [x[0] for x in placements]
    #Get the correspondance species name/taxid
    name_to_taxid = outils.get_name_to_taxid(all_sp_name, tax_tab)
    #Add the taxid at the end of each clade descriptor list
    new_placements = [(x[0],x[1], x[2], name_to_taxid[x[0]]) for x in placements]
    return new_placements



#Return the closest ancestor of a clade with more than a threshold of species in omamer
def get_sampled_taxa(clade, threshold_species, tax_tab, sp_tab, tax_buff,taxonomic_rank =None):
    ncbi = get_ete_ncbi_db()
    name_to_lineage = outils.get_full_lineage_omamer([clade], tax_tab)
    lineage = name_to_lineage[clade]

    is_taxonomic_rank = False
    selected_tax = None
    selected_phylum_or_higher = True
    selected_rank = None
    chosen_rank = False
    previous_seen = None

    name_to_taxid = outils.get_name_to_taxid([x for x in lineage], tax_tab)
    tax_clade = name_to_taxid[clade]
    # Getting the full lineage from NCBI
    lineage_ncbi = ncbi.get_lineage(tax_clade)
    ranks = ncbi.get_rank(lineage_ncbi)
    lineage_rank = [ranks.get(taxid, '') for taxid in lineage_ncbi]
    index_rank = -1
    if taxonomic_rank:
        if taxonomic_rank in lineage_rank:
            index_rank = lineage_rank.index(taxonomic_rank)

    for tax in lineage:
        species = outils.get_species_from_taxon(tax, tax_tab, sp_tab, tax_buff)
        taxid = name_to_taxid[tax]
        try:
            #Assessing the position of the given taxid in the lineage list, to compare with the rank that may be requested
            sp_rank_ncbi = lineage_ncbi.index(taxid)
        except ValueError:
            #If the taxid is not in index, setting it to a value that will make it avoid selecting the taxon
            sp_rank_ncbi = len(lineage_ncbi)
        rank = ncbi.get_rank([taxid]).get(taxid,'')
        if rank == taxonomic_rank:
            is_taxonomic_rank = True
        if rank == 'phylum' and selected_tax:
                #we assume all sub-phylum taxa are a subset of a phylum taxa, any clades selected below phylum will be marked as such
                selected_phylum_or_higher = False
        if len(species)>=threshold_species:
            if not selected_tax or rank==taxonomic_rank:
                selected_tax =  tax
                selected_rank = rank

                if rank==taxonomic_rank:
                    chosen_rank = True
                    is_taxonomic_rank = True
            elif not chosen_rank and sp_rank_ncbi<index_rank and previous_seen is not None:
                chosen_rank = True
                selected_tax = previous_seen[0]
                selected_rank = previous_seen[1]
        previous_seen = (tax, rank)

    if taxonomic_rank:
        if is_taxonomic_rank:
            LOG.info(f'The ancestral lineage is selected at provided taxonomic rank: {taxonomic_rank}.')
        elif chosen_rank:
            LOG.info(f"The chosen rank {taxonomic_rank} is not explicitely stored in the OMAmer database. Instead chosing the closest clade that share the same species set: {selected_tax}")
        else:
            LOG.info(f'The provided taxonomic rank {taxonomic_rank} was not an option (too narrow or absent from our lineage option). Default ancestral lineage will be used.')

    if selected_phylum_or_higher:
        LOG.warning("The selected ancestral lineage is from the phylum rank or higher which means the target species' taxonomic division is not well sampled in our database. The results may lack accuracy.")
    if selected_rank in ['genus', 'subgenus', 'section', 'subsection', 'species group', 'species subgroup', 'species']:
        LOG.warning(f'The selected ancestral lineage is from the {selected_rank} taxonomic rank. Consider trying a broader rank to validate your results.')
    return selected_tax


# USER DETERMINED PLACEMENT
def get_lineage_ncbi(taxid):
        lineage = list()
        ncbi = get_ete_ncbi_db()
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
                        return tax
        return None
