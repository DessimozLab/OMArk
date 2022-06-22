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

import argparse
import os
import omamer.database
import omark.files as io
import omark.species_determination as spd
import omark.omamer_utils as utils
import omark.scoring as sc
import omark.graphics as graph

def get_omamer_qscore(omamerfile, dbpath, stordir, taxid=None, contamination= True, original_FASTA_file = None, force = True):

    db = omamer.database.Database(dbpath)
    #Variables
    hog_tab = db._hog_tab[:]
    prot_tab = db._prot_tab
    sp_tab = db._sp_tab
    tax_tab = db._tax_tab[:]
    fam_tab = db._fam_tab[:]
    cprot_buff = db._cprot_arr
    tax_buff = db._ctax_arr[:]
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
            full_match_data, partials, fragments = io.filter_partial_matches(omamdata)


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
            conshog, cladehog = sc.get_conserved_hogs(closest_corr.decode(), hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,  cprot_buff,chog_buff, tax_buff, hogtax_buff, True, threshold=0.8)
            lineage_rhog = sc.get_root_HOGs_descendants(closest_corr, tax_tab, hog_tab, fam_tab,tax_buff)

            cladehog = cladehog + lineage_rhog

            #Two modes? : Normal and listing unexpected protein mapping?

            wholeres, found_clade, nic = sc.found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)

            #wholeres, whfound, nic = found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff)
           
            res_completeness, found_cons, nicons = sc.found_with_omamer(omamdata ,conshog, hog_tab, chog_buff)

            res_proteomes = sc.score_whole_proteome(found_clade, nic, partials, fragments, not_mapped, contaminant_prots)


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
            #Write graphical representation
            graph.plot_omark_results(stordir+'/'+basefile+".pdf", res_completeness, res_proteomes)

def launcher(args):
    omamerfile = args.file
    dbpath = args.database
    outdir = args.outputFolder
    taxid = args.taxid
    original_fasta = args.og_fasta
    get_omamer_qscore(omamerfile, dbpath, outdir, taxid, original_FASTA_file = original_fasta)

