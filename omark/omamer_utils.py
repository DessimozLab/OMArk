from omamer.hierarchy import get_descendants, get_leaves, get_root_leaf_offsets, get_children



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


def get_full_lineage_omamer(taxname, tax_tab, tax_buff = False,  descendant = False):
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
    if descendant :

        lineage += tax_tab[get_descendants(tax2tax_off[taxname], tax_tab, tax_buff)]['ID'].tolist()
    return lineage

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
