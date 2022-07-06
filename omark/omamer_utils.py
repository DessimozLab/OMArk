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

from omamer.hierarchy import get_descendants, get_leaves, get_root_leaf_offsets, get_children
from tables.exceptions import HDF5ExtError
import omamer.database
from omark.utils import LOG

def check_database(dbpath):
    #Check errors in the database parameter.
    valid = True
    try:
        db = omamer.database.Database(dbpath)
        hog_tab = db._hog_tab
        prot_tab = db._prot_tab
        sp_tab = db._sp_tab
        tax_tab = db._tax_tab
        fam_tab = db._fam_tab
        cprot_buff = db._cprot_arr
        tax_buff = db._ctax_arr
        chog_buff = db._chog_arr
        hogtax_buff = db._hog_taxa_buff
        db.close()

    except OSError:
        LOG.error('Path to the OMAmer database is not valid.')
        valid = False
    except HDF5ExtError:
        LOG.error('The OMAmer database is not a valid HDF5 file.')
        valid = False
    except AttributeError:
        LOG.error('The provided HDF5 database is not a correct OMAmer database.')
        valid = False
        db.close()
    return valid

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

    return hog_taxa


#Could be more efficient. For now fetch multiple time children of Tax. See if need improvement
#Return a dictionary of set of species by clade.
def get_spec_by_tax(tax_tab, sp_tab, tax_buff):
    i=0
    spec_by_tax = dict()
    for tax_entry in tax_tab:
        taxname = tax_entry['ID']
        sp_off_in_tax = get_leaves(i, tax_tab, tax_buff)
        sp_tax =[ tax_tab[x][0].decode() for x in sp_off_in_tax]
        spec_by_tax[taxname] = set(sp_tax)
        i+=1
    return spec_by_tax


def get_nb_hogs_by_clade(hog_tab, tax_tab):
    hog_by_tax = dict()
    for hog_off in range(hog_tab.size):
        taxoff = hog_tab[hog_off]['TaxOff']
        tax = tax_tab[taxoff]["ID"]
        if tax not in hog_by_tax:
            hog_by_tax[tax]=0
        hog_by_tax[tax]+=1
    return hog_by_tax

def get_prop_duplicated(hog_tab, tax_tab, chog_buff):
    transition = dict()
    hog_by_tax = dict()
    for hog_off in range(hog_tab.size):
        t = hog_tab[hog_off]
        taxoff = hog_tab[hog_off]['TaxOff']
        tax = tax_tab[taxoff]["ID"]
        if tax not in hog_by_tax:
            hog_by_tax[tax]=0
        descendants = get_descendant_HOGs(t, hog_tab, chog_buff)
        if len(descendants)>0:
            if tax not in transition:
                transition[tax] = dict()
            seen_sp = list()        
            for desc in descendants:
                
                taxd_off = desc['TaxOff']
                tax_d = tax_tab[taxd_off]['ID']
                if tax_d not in transition[tax]:
                    transition[tax][tax_d] = 0
                if not tax_d in seen_sp:
                    transition[tax][tax_d] += 1
                    seen_sp.append(tax_d)
        hog_by_tax[tax]+=1

    prop_duplicated = dict()
    for x in transition:
        tot = hog_by_tax[x]
        if x not in prop_duplicated:
            prop_duplicated[x] = dict()
        for desc in transition[x]:
            prop_duplicated[x][desc] = transition[x][desc]/hog_by_tax[x]
    return prop_duplicated


def get_species_from_taxon(taxname, tax_tab, sp_tab, tax_buff):
    tax_off2tax = tax_tab['ID'] 
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    if type(taxname)==str:
        taxname = taxname.encode('ascii')
    tax_off = tax2tax_off[taxname]
    sp_off_in_tax = get_leaves(tax_off, tax_tab, tax_buff)
    sp_tax =[ tax_tab[x][0].decode() for x in sp_off_in_tax]
    return sp_tax


def get_full_lineage_omamer(taxnames, tax_tab, tax_buff = False,  descendant = False):
    name_to_lineage = dict()
    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    for taxname in taxnames:
        lineage = list()
        reached = False
        current_tax = tax_tab[tax2tax_off[taxname]]
        while not reached: 
            lineage.append(current_tax['ID'])
            ancestor_tax = current_tax['ParentOff']
            if ancestor_tax!=-1:
                    current_tax  = tax_tab[ancestor_tax]
            else:
                    reached = True
        if descendant :

            lineage += tax_tab[get_descendants(tax2tax_off[taxname], tax_tab, tax_buff)]['ID'].tolist()
        name_to_lineage[taxname] = lineage
        
    return name_to_lineage

def get_name_to_taxid(taxnames, tax_tab):
    name_to_taxid = dict()
    tax_off2tax = tax_tab['ID']
    tax2tax_off = dict(zip(tax_off2tax, range(tax_off2tax.size)))
    for name in taxnames:
        tax = tax_tab[tax2tax_off[name.encode()]]
        taxid = tax['TaxID']
        name_to_taxid[name] = taxid
    return name_to_taxid

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

def get_species_from_omamer(hog, prot_tab, spe_tab, cprot_buff) :
    sp_list = list()
    chog_off = hog["ChildrenProtOff"]
    prots = cprot_buff[chog_off : chog_off + hog["ChildrenProtNum"]]

    for p in prots:        
        spe_off = prot_tab[p][1]
        sp_list.append(spe_tab[spe_off])

    return sp_list
