import pyoma.browser.db
import argparse
import os
import sys
import time
import pyoma
import pyham
import inspect
import urllib 
import sys
import gzip
import Bio
import ete3
import re

import pandas
if sys.version_info[0] == 3:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
import omamer
import omamer.database
from omamer.hierarchy import get_descendant_taxa, _children_hog
import omamer_species_placement as osp
import taxonomic_placement as tp


def build_arg_parser():
	"""Handle the parameter sent when executing the script from the terminal

	Returns
	-----------
	A parser object with the chosen option and parameters"""

	parser = argparse.ArgumentParser(description="Compute an OMA quality score from the OMAmer file of a proteome.")   
	parser.add_argument('-f', '--file', help="The OMAmer file to read." )	
	parser.add_argument('-d', '--database', help="The OMAmer database.")
	parser.add_argument('-o', '--outputFolder', help="The folder containing output data the script wilp generate.")
	parser.add_argument('-m', '--oma', help="Path to the OMA Database")

	return parser



def parseOmamer(file):
    alldata = list()
    with open(file) as f:
        
        firstline = f.readline()
        cat = firstline.strip('\n').split('\t')
        for line in f.readlines():
            data = dict()
            col = line.strip('\n').split('\t')
            for i in range(len(cat)) :
                data[cat[i]] = col[i]
            alldata.append(data)
    return alldata


#def getCloseTaxaOMAm():
#alltaxa = dict()
#for omamapping in omamerdats:



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
    
	


def get_close_taxa_omamer(omamerdata, hog_tab, tax_tab, ctax_buff, chog_buff):
    
    alltaxa = dict()
    j=0
    descendant = None
    

    hog_off2subf = hog_tab['OmaID'] 
    subf2hog_off = dict(zip(hog_off2subf, range(hog_off2subf.size)))
    tax_off2tax = tax_tab['ID'] 
    
    for omamapping in omamerdata:
        j+=1
	
        hog_off = subf2hog_off[omamapping['subfamily'].encode('ascii')]       
        taxa = get_hog_implied_taxa(hog_off, hog_tab, tax_tab, ctax_buff, chog_buff)
        tax_off2tax = tax_tab['ID'] 
        for taxon in taxa:
            taxname = tax_off2tax[taxon]
            if taxname in alltaxa :
                alltaxa[taxname]+=1
            else: 
                alltaxa[taxname]=1
    print(len(alltaxa))
    print(alltaxa)
    
    alltaxa = {k: v for k, v in sorted(alltaxa.items(), key=lambda item: item[1], reverse=True)}
    #alltaxa = { k:v for k in sorted(alltaxa.iteritems(), key=itemgetter(1), reverse=True)}
    return alltaxa

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
    sp_off_in_tax = omamer.hierarchy.get_descendant_species(tax_off, tax_tab, tax_buff)
    sp_tax =[ sp_tab[x][0].decode() for x in sp_off_in_tax]
    return sp_tax

def get_species_from_omamer(hog, prot_tab, spe_tab, cprot_buff) :
    sp_list = list()
    chog_off = hog["ChildrenProtOff"]
    prots = cprot_buff[chog_off : chog_off + hog["ChildrenProtNum"]]

    for p in prots:        
        spe_off = prot_tab[p][3]
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
    hog_off = hog['ChildrenHOGoff'] 
    hog_num = hog['ChildrenHOGnum']
    desc_hog = chog_buff[hog_off : hog_off+hog_num]
    subhogs = [ hog_tab[x] for x in desc_hog]
    all_hogs += subhogs
    for subhog in subhogs :
        all_hogs += get_descendant_HOGs(subhog, hog_tab, chog_buff)
    return all_hogs



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
def get_conserved_hogs(clade, hog_tab, prot_tab, sp_tab, tax_tab, cprot_buff, chog_buff, tax_buff, duplicate ) :
    found_hog = list()
    lineage = get_full_lineage_omamer(clade.encode('ascii'), tax_tab)
    sp_target = get_species_from_taxon(clade, tax_tab, sp_tab, tax_buff)

    for t in hog_tab:
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
        if len(inter)>0.9*len(sp_target):
            found_hog.append(t)
        #omamer.hierarchy.get_descendant_species_taxoffs(hog_off, tabi, chog_buff, cprot_buff, prot2speoff, speoff2taxoff
    return found_hog


def found_with_omamer(omamer_data, conserved_hogs, hog_tab, chog_buff):
    all_subf = list()
    results = { 'Found':[] , 'Lost' : [], 'Duplicated': [], 'Underspecific':[], 'Overspecific': []}
    for data in omamer_data:
        
        all_subf.append(data['subfamily'])

    for hog in conserved_hogs :
        identifier = hog['OmaID'].decode()
        if identifier in all_subf:
            nb_found = all_subf.count(identifier)
            if nb_found > 1:
                results['Duplicated'].append(identifier)

            else :
                results['Found'].append(identifier)

        else:
            in_subhog = False
            
            for subhog in [x['OmaID'].decode() for x in get_descendant_HOGs(hog, hog_tab, chog_buff)]:


                if subhog in all_subf:
                    results['Overspecific'].append(identifier)
                    in_subhog = True
            if not in_subhog:                   
                in_superhog = False
                for superhog in [x['OmaID'].decode() for x in get_ancestral_HOGs(hog, hog_tab, chog_buff)]:
                    if superhog in all_subf:

                        results['Underspecific'].append(identifier)
                        in_superhog = True
                        break

                if not in_superhog:
                    results['Lost'].append(identifier)
    return results

def get_omamer_qscore(omamerfile,db, omadbpath,  stordir):
	
	#Variables
	hog_tab = db._hog_tab[:]
	prot_tab = db._prot_tab
	sp_tab = db._sp_tab
	tax_tab = db._tax_tab[:]
	cprot_buff = db._cprot_arr
	tax_buff = db._ctax_arr
	chog_buff = db._chog_arr
	
	allres = dict()
	#Store the temporary results in a file to avoid recomputing and make it computationally feasible
	if os.path.isfile(omamerfile):
		taxid = omamerfile.split('/')[-1].split('_')[1].strip('.fasta')
		#print('Working with '+omamerfile)
		#print(taxid)
		#print(stordir+omamerfile.split('/')[-1].strip('.fasta')+".omq")
		#print("Before working")
		if not os.path.isfile(stordir+omamerfile.split('/')[-1].strip('.fasta')+".omq"): 

			#print('Parse OMAmer')
			omamdata = parseOmamer(omamerfile)
			#print(len(omamdata))
			#print('get Close')
			close = get_close_taxa_omamer(omamdata, hog_tab, tax_tab, tax_buff, chog_buff)
			#close = getCloseTaxa(omamdata, omadbpath)
			#print('Close taxa found')
			closest =  get_lower_noncontradicting(close, tax_tab)
			closest_corr = osp.get_sampled_taxa(closest, 2 , tax_tab, sp_tab, tax_buff)
			conshog = get_conserved_hogs(closest_corr.decode(), hog_tab, prot_tab, sp_tab, tax_tab, cprot_buff,chog_buff, tax_buff, True)
			
			res = found_with_omamer(omamdata ,conshog, hog_tab, chog_buff)
			store_close_level(stordir+omamerfile.split('/')[-1].strip('.fasta')+".tax", {'Sampled': str(closest_corr.decode()), 
													'Closest' : str(closest.decode()),
													'All' : close })
			store_results(stordir+omamerfile.split('/')[-1].strip('.fasta')+".omq", res)   

def store_results(storfile, results):
	with open(storfile, 'w') as storage:
		for categ, hoglist in results.items():
			storage.write('>'+categ+'\n')
			for elem in hoglist:
				storage.write(elem+'\n')

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
	db = omamer.database.Database(dbpath)
	print(db)
	outdir = arg.outputFolder
	print(outdir)
	omadb = arg.oma
	print(omadb)
	get_omamer_qscore(omamerfile, db, omadb ,  outdir)
	print('Done')

