import Bio
import re


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

def store_results(storfile, results):
    with open(storfile, 'w') as storage:
        for categ, hoglist in results.items():
            storage.write('>'+categ+'\n')
            for elem in hoglist:
                storage.write(elem+'\n')

def store_summary(storfile, results, results_proteomes, selected_lineage, contaminant = False, prot_clade = False):
    with open(storfile,'w') as storage:
        storage.write('#The selected clade was '+selected_lineage.decode()+"\n")
        total = len(results['Single'])+len(results['Duplicated'])+ len(results['Overspecific_S']) + len(results['Overspecific_D'])+ len(results['Underspecific']) + len(results['Lost'])
        storage.write('#Number of conserved HOGs is: '+str(total)+'\n')
        tot_genes = len(results_proteomes['Not_Placed'])+len(results_proteomes['Correct'])+len(results_proteomes['Erroneous'])
        storage.write('#Results on conserved HOGs is:\n')
        storage.write('#S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing\n')
        storage.write(f'S:{len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"])},D:{len(results["Duplicated"])+len(results["Overspecific_D"])}[U:{len(results["Duplicated"])},E:{len(results["Overspecific_D"])}],M:{len(results["Lost"])}\n') 
        storage.write(f'S:{100*(len(results["Single"])+len(results["Overspecific_S"])+len(results["Underspecific"]))/total:4.2f}%,D:{100*(len(results["Duplicated"])+len(results["Overspecific_D"]))/total:4.2f}%[U:{100*len(results["Duplicated"])/total:4.2f}%,E:{100*len(results["Overspecific_D"])/total:4.2f}%],M:{100*len(results["Lost"])/total:4.2f}%\n') 
        storage.write('#On the whole proteome, there is '+str(tot_genes)+' proteins\n')
        storage.write('#Of which:\n')
        storage.write('#A:Placements in accurate lineage[P:Partial hits,F:Fragmented], E: Erroneous placements[P:Partial hits,F:Fragmented], C: Likely contamination[P:Partial hits,F:Fragmented], N: no mapping \n')
        storage.write(f'A:{len(results_proteomes["Correct"])}[P:{len(results_proteomes["Correct_Partial"])},F:{len(results_proteomes["Correct_Fragment"])}],E:{len(results_proteomes["Erroneous"])}[P:{len(results_proteomes["Erroneous_Partial"])},F:{len(results_proteomes["Erroneous_Fragment"])}],C:{len(results_proteomes["Contamination"])}[P:{len(results_proteomes["Contamination_Partial"])},F:{len(results_proteomes["Contamination_Fragment"])}],N:{len(results_proteomes["Not_Placed"])}\n')
        storage.write(f'A:{100*len(results_proteomes["Correct"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Correct_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Correct_Fragment"])/tot_genes:4.2f}%],E:{100*len(results_proteomes["Erroneous"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Erroneous_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Erroneous_Fragment"])/tot_genes:4.2f}%],C:{100*len(results_proteomes["Contamination"])/tot_genes:4.2f}%[P:{100*len(results_proteomes["Contamination_Partial"])/tot_genes:4.2f}%,F:{100*len(results_proteomes["Contamination_Fragment"])/tot_genes:4.2f}%],N:{100*len(results_proteomes["Not_Placed"])/tot_genes:4.2f}%\n')
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
