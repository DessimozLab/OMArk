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
import omark.omamer_utils as outils
from omamer.hierarchy import get_descendants



#Mutliple level of a same HOG can be counted as is   
def get_conserved_hogs(clade, hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,   cprot_buff, chog_buff, tax_buff, hogtax_buff,  duplicate, threshold=0.9 ) :
    found_hog = list()
    poss_hog = list()
    seen_hog = list()
    other_cl_hog = list()
    clade_to_lineage = outils.get_full_lineage_omamer([clade.encode('ascii')], tax_tab, tax_buff, True)
    lineage = clade_to_lineage[clade.encode('ascii')]
    sp_target = outils.get_species_from_taxon(clade, tax_tab, sp_tab, tax_buff)

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
                all_desc = outils.get_descendant_HOGs(t, hog_tab, chog_buff)
                for desc in all_desc:
                        desc_tax_name = tax_tab[desc['TaxOff']]["ID"]
                        if desc_tax_name not in lineage :
                               continue
                        sp_hog += [x[0].decode() for x in outils.get_species_from_omamer(desc,prot_tab, sp_tab, cprot_buff)]
        sp_hog += [x[0].decode() for x in outils.get_species_from_omamer(t,prot_tab, sp_tab, cprot_buff)]
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
    descendants = tax_tab[get_descendants(tax2tax_off[lineage], tax_tab, tax_buff)]['ID']
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
        for subhog in [x['OmaID'].decode() for x in outils.get_descendant_HOGs(hog, hog_tab, chog_buff)]:
        
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
            
            
        for superhog in [x['OmaID'].decode() for x in outils.get_ancestral_HOGs(hog, hog_tab, chog_buff)]:
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
