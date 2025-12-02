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


import omamer
import omamer.database
from omamer.hierarchy import get_leaves
import argparse
import numpy as np

#Get the closest taxa based on omamerdata (A list of subfamily ID). For each Hog, get their level and Return an ordered  dict with 
#taxonomic level as keys and their count as value
def get_close_taxa(omamerdata, hog_tab, tax_tab, ctax_buff, chog_buff):

    alltaxa = dict()
    descendant = None


    hog_off2subf = hog_tab['OmaID']
    subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))    
    tax_off2tax = tax_tab['ID']
    for omamapping in omamerdata:
        hog_off = subf2hog_off[omamapping['subfamily'].encode('ascii')]
        taxon = hog_tab[hog_off]['TaxOff']
        
        taxname = tax_off2tax[taxon]
        if taxname in alltaxa :
            alltaxa[taxname]+=1
        else:
            alltaxa[taxname]=1
    #Ordering
    alltaxa = {k: v for k, v in sorted(alltaxa.items(), key=lambda item: item[1], reverse=True)}
    return alltaxa

def get_queriesp_from_output(omamerdata, hog_tab):
    q2hog_off = np.full(len(omamerdata, -1, dtype=np.int64))
    hog_off2subf = hog_tab['OmaID']
    subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
    tax_off2tax = tax_tab['ID']
    for i in range(len(omamerdata)):
        hog_off = subf2hog_off[omamapping[i]['subfamily'].encode('ascii')]
        q2hog_off[i] = np.int(hog_off)	
    return q2hog_off


def _place_queries(query_offsets, q2fam_off, q2fam_score, fam_tab, q2hog_bestpath, q2hog_scores, fst, sst):
    q2hog_off = np.full(q2fam_off.size, -1, dtype=np.int64)
    for i in numba.prange(query_offsets.size):
        q = query_offsets[i]

        # if below family threshold or equal to 0, skip 
        if (q2fam_score[q] < fst) or (q2fam_score[q] == 0):
            continue

        # find best scoring subfamily
        fam_hog_nr = fam_tab[q2fam_off[q]]['HOGnum']
        fam_bestpath = q2hog_bestpath[q, :fam_hog_nr]
        fam_hog_scores = q2hog_scores[q, :fam_hog_nr]
        best_j = None
        best_s = np.inf
        for j in np.argwhere(fam_bestpath).flatten():
            s = fam_hog_scores[j]
            if s < best_s and s >= sst:
                best_j = j
                best_s = s
        if best_j is not None:
            root_hog_off = fam_tab[q2fam_off[q]]['HOGoff']
            q2hog_off[i] = np.int(best_j + root_hog_off)
    return q2hog_off

#Based on an ordered dict of taxa, get the lowest taxonomic level that is not in conflict with higher level
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

#Return the full lineage of a clade based on OMAmer taxonomic data, as an array of taxa (ordered from initial node to root)
def get_full_lineage_omamer(taxname, tax_tab):
    lineage = list()
    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    reached = False

    current_tax = tax_tab[tax2tax_off[taxname]]
    while not reached:
        lineage.append(current_tax['ID'])
        ancestor_tax = current_tax['ParentOff']
        if ancestor_tax!=-1:
                current_tax  = tax_tab[ancestor_tax]
        else:
                reached = True
    return lineage

#Return the closest ancestor of a clade with more than a threshold of species in omamer
def get_sampled_taxa(clade, threshold_species, tax_tab, sp_tab, tax_buff):
    lineage = get_full_lineage_omamer(clade, tax_tab)
    for tax in lineage:
        species = get_species_from_taxon(tax, tax_tab, sp_tab, tax_buff)
        if len(species)>=threshold_species:
              return tax
    return None

#Return a list of species name, corresponding to a given taxonomic clade (given as name)
def get_species_from_taxon(taxname, tax_tab, sp_tab, tax_buff):

    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    tax_off = tax2tax_off[taxname]
    sp_off_in_tax = get_leaves(tax_off, tax_tab, tax_buff)
    sp_tax = [ tax_tab[x][0] for x in sp_off_in_tax]

    return sp_tax


#---------------------






