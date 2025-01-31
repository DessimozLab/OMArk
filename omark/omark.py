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
import omark.omamer_utils as outils
import omark.scoring as sc
import omark.graphics as graph
from omark.utils import LOG, set_log_level
import sys

def get_omamer_qscore(omamerfile, dbpath, stordir, taxid=None, contamination= True, original_FASTA_file = None, force = True, isoform_file = None, taxonomic_rank= None, extract_conserved_only=True):



	   
    allres = dict()
    
    basefile = '.'.join(omamerfile.split('/')[-1].split('.')[:-1])
    if force or not os.path.isfile(stordir+'/'+basefile+".omq"): 
        #Extract OMAmer data
        LOG.info('Extracting data from input file: '+omamerfile)

        omamdata, not_mapped  = io.parseOmamer(omamerfile)

        #Check omamer version:
        if "hoglevel" in omamdata[0]:
            omamer_version = "2.0.0"
        else:
            omamer_version = "0.2.0"
        db = omamer.database.Database(dbpath)
        if omamer_version == '2.0.0':
            hog_tab = db._db_HOG[:]  # db._hog_tab[:]
            prot_tab = db._db_Protein  # db._prot_tab
            sp_tab = db._db_Species  # db._sp_tab
            tax_tab = db._db_Taxonomy[:]  # db._tax_tab[:]
            fam_tab = db._db_Family[:]  # db._fam_tab[:]
            cprot_buff = db._db_ChildrenProt  # db._cprot_arr
            tax_buff = db._db_ChildrenTax[:]  # db._ctax_arr[:]
            chog_buff = db._db_ChildrenHOG  # db._chog_arr
            hogtax_buff = db._db_HOGtaxa  # db._hog_taxa_buff
            hog_id_buff = db._db_HOGIDBuffer[:]
        else:
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
            hog_id_buff = None

        if isoform_file:
            LOG.info('An isoform_file was provided.')
            LOG.info('Extracting data from isoform file '+isoform_file)
            isoform_data = io.parse_isoform_file(isoform_file)
            omamdata, not_mapped, selected_isoforms = io.select_isoform(isoform_data, omamdata, omamer_version)

        #Get only full match for placements
        full_match_data, partials, fragments = io.filter_partial_matches(omamdata)

        LOG.info('Determinating species composition from HOG placements')

        #Determine species and contamination
        if omamer_version == "2.0.0":
            placements = spd.get_present_lineages(full_match_data, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff)
        else:
            placements = spd.get_present_lineages(full_match_data, hog_tab, tax_tab, tax_buff, sp_tab, chog_buff, family_score_filter=None, cutoff_percentage=0.001)


        #Get the proteins placed in species correspoding to each placement in a dictionary
        prot_clade = spd.get_prot_by_clades(placements, omamdata, hog_tab, tax_tab, tax_buff, chog_buff)
        #Reorganize the placements to consider the species with most proteins to be the main one. (Needed in edge cases where the proteins of the contaminant
        #is an exhaustive set).
        placements = spd.reorganized_placement(placements, prot_clade)
        
        #Procedure when the user do not give taxonomu information. Will use the main species from the placement
        if taxid==None:

            likely_clade =  placements[0][0]
            LOG.info('No taxid was provided. From HOG placements, the query taxon is '+likely_clade)

        #Otherwise find the closest clade in the OMAmer database from the given taxid
        else :
            lin = spd.get_lineage_ncbi(taxid)
            likely_clade = spd.find_taxa_from_ncbi(lin, tax_tab, sp_tab,tax_buff)
            #Used in case the placement main species conflict with the given taxid, will chose the given clade as reference
            placements = spd.reorganize_placements_from_taxid(placements, likely_clade,tax_tab, tax_buff)
            LOG.info('A taxid was provided. The query taxon is '+likely_clade)

        #Get the proteins belonging to contaminants
        contaminant_prots = spd.get_contaminant_proteins(placements, prot_clade)
        #Add the taxid information to the species description list.
        placements = spd.add_taxid(placements, tax_tab)
        placements, prot_clade = spd.add_uncertain_contaminants(placements, prot_clade)
        #Get the first parent of the chosen clade with at least 5 species
        closest_corr = spd.get_sampled_taxa(likely_clade, 5 , tax_tab, sp_tab, tax_buff, taxonomic_rank)

        LOG.info('Ancestral lineage is '+closest_corr)

        #Conshog : HOG with 80% representative of the target lineage
        #Cladehog : HOG with at least 1 representative of the target lineage present in the common ancestir
        LOG.info('Estimating ancestral and conserved HOG content')

        conshog, cladehog = sc.get_conserved_hogs(closest_corr, hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,  cprot_buff,chog_buff, tax_buff, hogtax_buff, True, threshold=0.8)
        lineage_rhog = sc.get_root_HOGs_descendants(closest_corr, tax_tab, hog_tab, fam_tab,tax_buff)

        cladehog = cladehog + lineage_rhog

        LOG.info(str(len(cladehog))+" HOGs are associated to the query's lineage and will be used for consistency assesment")
        LOG.info(str(len(conshog))+' conserved ancestral HOGs will be used from completeness assesment')

        LOG.info('Comparing the query gene repertoire to lineage-associated HOGs')
        wholeres, found_clade, nic = sc.found_with_omamer(omamdata ,cladehog, hog_tab, chog_buff, hog_id_buff)
        res_proteomes = sc.score_whole_proteome(found_clade, nic, partials, fragments, not_mapped, contaminant_prots)

        LOG.info('Comparing the query gene repertoire to conserved ancestral HOGs')
        res_completeness, found_cons, nicons = sc.found_with_omamer(omamdata ,conshog, hog_tab, chog_buff, hog_id_buff)

        db.close()
        LOG.info('Writing OMArk output files')

        #Store the taxonomic choices
        io.store_close_level(stordir+'/'+basefile+".tax", {'Sampled': str(closest_corr),
                                                                            'Closest' : str(likely_clade) })

        #Write optionnal taxa files
        if original_FASTA_file:
            if contamination :
                io.store_contaminant_FASTA(stordir, basefile, prot_clade, original_FASTA_file)
            io.store_incorrect_map_FASTA(stordir, basefile, not_mapped, nic, original_FASTA_file)


        if isoform_file:
            io.store_list(stordir+'/'+basefile+"_selected_isoforms.txt", selected_isoforms)
        #Write results files
        io.store_results(stordir+'/'+basefile+".ump", { key : res_proteomes[key] for key in ['Consistent_Full', 'Consistent_Partial', 'Consistent_Fragment',
                                                                                             'Inconsistent_Full', 'Inconsistent_Partial', 'Inconsistent_Fragment',
                                                                                             'Contamination_Full', 'Contamination_Partial', 'Contamination_Fragment',
                                                                                             'Unknown']})
        io.store_results(stordir+'/'+basefile+".omq", res_completeness) 
        io.write_templated_report('summarized_report.txt', stordir+'/'+basefile+".sum", res_completeness, res_proteomes, closest_corr, placements)

        io.write_templated_report('textual_report.txt', stordir+'/'+basefile+"_detailed_summary.txt", res_completeness, res_proteomes, closest_corr, placements)
        #Write graphical representation
        plot_res = {'pdf': os.path.join(stordir, basefile + ".pdf"),
                    'png': os.path.join(stordir, basefile + ".png")}
        graph.plot_omark_results(plot_res, res_completeness, res_proteomes)



def get_only_conserved_HOGs(dbpath, stordir, taxid, taxonomic_rank=None):
    db = omamer.database.Database(dbpath)

    #Variables
    if hasattr(db, '_hog_tab'):
        hog_tab = db._hog_tab[:]
        prot_tab = db._prot_tab
        sp_tab = db._sp_tab
        tax_tab = db._tax_tab[:]
        fam_tab = db._fam_tab[:]
        cprot_buff = db._cprot_arr
        tax_buff = db._ctax_arr[:]
        chog_buff = db._chog_arr
        hogtax_buff = db._hog_taxa_buff
        hog_id_buff = None
    else:
        hog_tab = db._db_HOG[:]  # db._hog_tab[:]
        prot_tab = db._db_Protein  # db._prot_tab
        sp_tab = db._db_Species  # db._sp_tab
        tax_tab = db._db_Taxonomy[:]  # db._tax_tab[:]
        fam_tab = db._db_Family[:]  # db._fam_tab[:]
        cprot_buff = db._db_ChildrenProt  # db._cprot_arr
        tax_buff = db._db_ChildrenTax[:]  # db._ctax_arr[:]
        chog_buff = db._db_ChildrenHOG  # db._chog_arr
        hogtax_buff = db._db_HOGtaxa  # db._hog_taxa_buff
        hog_id_buff = db._db_HOGIDBuffer[:]
    lin = spd.get_lineage_ncbi(taxid)
    likely_clade = spd.find_taxa_from_ncbi(lin, tax_tab, sp_tab,tax_buff)
    LOG.info('A taxid was provided. The query taxon is '+likely_clade)

    closest_corr = spd.get_sampled_taxa(likely_clade, 5 , tax_tab, sp_tab, tax_buff, taxonomic_rank)
    LOG.info('Ancestral lineage is '+closest_corr)
    LOG.info('Estimating ancestral HOG content')
    conshog, cladehog = sc.get_conserved_hogs(closest_corr, hog_tab, prot_tab, sp_tab, tax_tab, fam_tab,  cprot_buff,chog_buff, tax_buff, hogtax_buff, True, threshold=0.8)
    LOG.info(str(len(conshog))+" HOGs are conserved in the query's lineage.")


    hog_list = [outils.get_hog_id(x, hog_id_buff) for x in conshog]
    io.store_list(stordir+'/conserved_HOGs.txt', hog_list, comment=[f'Ancestral lineage: {closest_corr}'])
    db.close()

def check_parameters(omamerfile, dbpath, stordir, taxid=None, original_FASTA_file = None,  isoform_file = None, taxonomic_rank=None):

    omamerfile_valid = io.check_omamerfile(omamerfile)

    database_valid = outils.check_database(dbpath)

    taxid_valid = spd.check_taxid(taxid)

    output_directory_valid = io.check_and_create_output_folder(stordir)

    if original_FASTA_file:
        fasta_valid = io.check_FASTA(original_FASTA_file, omamerfile)
    else:
        fasta_valid = True

    if isoform_file:
        isoform_valid = io.check_isoform_file(isoform_file, omamerfile)
    else:
        isoform_valid = True

    if taxonomic_rank:
        taxonomic_rank = spd.check_rank(taxonomic_rank)
    else:
        taxonomic_rank = True

    return omamerfile_valid and database_valid and output_directory_valid and taxid_valid and fasta_valid and isoform_valid and taxonomic_rank

def launcher(args):
    omamerfile = args.file
    dbpath = args.database
    outdir = args.outputFolder
    taxid = args.taxid
    original_fasta = args.og_fasta
    isoform_file = args.isoform_file
    taxonomic_rank = args.taxonomic_rank
    verbose = args.verbose
    only_conserved_HOGs = args.output_cHOGs
    ete_ncbi_db = args.ete_ncbi_db
    log_level = 'INFO' if verbose else 'WARNING'
    set_log_level(log_level)
    LOG.info('Starting OMArk')
    if ete_ncbi_db is not None:
        LOG.info(f'A custom path was offered fot the ete3 NCBI database. {ete_ncbi_db} will be used.')
        spd.set_ete_taxa_path(ete_ncbi_db)
    if only_conserved_HOGs: 
        LOG.info('The option to output only conserved_HOGs was selected')
        if not taxid :
            LOG.error('A taxid must be provided for this use-case')
            sys.exit(1)
        if spd.check_taxid(taxid) and io.check_and_create_output_folder(outdir):
            get_only_conserved_HOGs(dbpath, outdir, taxid,taxonomic_rank=taxonomic_rank )
            LOG.info('Done')
        else:
            LOG.error('Exiting because one or more parameters are incorrect')
            sys.exit(1)
    elif check_parameters(omamerfile, dbpath, outdir,taxid,original_fasta,isoform_file,taxonomic_rank):
        LOG.info('Input parameters passed validity check')
        get_omamer_qscore(omamerfile, dbpath, outdir, taxid, original_FASTA_file = original_fasta, isoform_file=isoform_file, taxonomic_rank=taxonomic_rank)
        LOG.info('Done')

    else:
        LOG.error('Exiting because one or more parameters are incorrect')
        sys.exit(1)

