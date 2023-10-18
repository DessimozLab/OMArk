from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
import random
from collections import Counter
from numpy.random import choice


def build_arg_parser():
	"""Handle the parameter sent when executing the script from the terminal

	Returns
	-----------
	A parser object with the chosen option and parameters"""

	parser = argparse.ArgumentParser(description="Generate simulated proteomes from an input directory of proteomes.")   
	parser.add_argument('-i', '--input', help="The directory with input genomes." )	
	parser.add_argument('-o', '--output', help="The directory with output genomes.")
	parser.add_argument('-c', '--contamnants', help="The directory with input proteomes." )
	parser.add_argument('-r', '--reuse', action="store_true", help="Sign to use the minimalized files to generate whole genomes.")
	return parser


def load_proteome(proteome_file):
	return list(SeqIO.parse(proteome_file, "fasta"))

def write_proteome(output, records):
	with open(output,'w') as of:
		SeqIO.write(records, of, 'fasta')


def add_contamination(target_prot, contaminant_prot, genes_number):
	cont_sampling = random.sample(contaminant_prot,genes_number)
	contaminated_prot = target_prot+cont_sampling
	return contaminated_prot, cont_sampling

def draw_partial_proteome(target_prot, proportion):
	number_genes  = int(len(target_prot)*proportion)
	partial_prot = random.sample(target_prot,number_genes)
	return partial_prot

def add_random_sequences(target_prot, proportion, frequency=None):
	number_genes  = int(len(target_prot)*proportion)
	if not frequency:
		added_prot = [generate_random_prot(x) for x in range(number_genes)]
	else:
		added_prot = [generate_random_prot_from_aa_freq(x, frequency) for x in range(number_genes)]

	expanded_prot = target_prot + added_prot
	return expanded_prot, added_prot

def create_fusions(target_prot, proportion):
	removed_seq = []
	fused_prots = []
	number_genes = int(len(target_prot)*proportion)
	number_genes = number_genes if number_genes%2==0 else number_genes-1
	new_rec_list = target_prot.copy()
	subset_seq = random.sample(target_prot,number_genes)
	for i in range(0, len(subset_seq), 2):
		fseq = subset_seq[i]
		sseq = subset_seq[i+1]
		new_seq = SeqRecord(make_fragmented_sequence(fseq.seq, min_size=0.8, max_size=1, side='N')+make_fragmented_sequence(sseq.seq, min_size=0.8, max_size=1, side='C'), id=fseq.id+"_with_"+sseq.id)
		removed_seq.append(fseq.id)
		removed_seq.append(sseq.id)
		new_rec_list.append(new_seq)
		fused_prots.append(new_seq)
	new_rec_list = [x for x in new_rec_list if x.id not in removed_seq ]
	return new_rec_list, fused_prots

def create_fragments(target_prot, proportion):
	removed_seq = []
	fragmented_prot = []
	new_rec_list = target_prot.copy()

	subset_seq = random.sample(target_prot, int(len(target_prot)*proportion))
	for record in subset_seq:
		new_seq = make_fragmented_sequence(record.seq, min_size=0.1, max_size=0.9) 
		removed_seq.append(record.id)
		new_rec_list.append(SeqRecord(new_seq, id=record.id+'_fragments'))
		fragmented_prot.append(SeqRecord(new_seq, id=record.id+'_fragments'))
	new_rec_list = [x for x in new_rec_list if x.id not in removed_seq ]
	return new_rec_list, fragmented_prot

def regenerate_fragmented(target_prot, fragmented_prot):
	fragmented_ids = [x.id.replace('_fragments','') for x in fragmented_prot]
	new_rec_list = [x for x in target_prot if x.id not in fragmented_ids]
	new_proteome = new_rec_list + fragmented_prot
	return new_proteome

def regenerate_fused(target_prot, fused_prots):
	fused_ids = []
	for x in fused_prots:
		fused_identifier = new_seq.split('_with_')
		fident = fused_identifier[0]
		fused_ids.append(fident)
		sident = fused_identifier[1]
		fused_ids.append(sident)

	new_rec_list = [x for x in target_prot if x.id not in fused_ids]
	new_proteome = new_rec_list + fused_prots
	return new_proteome

#Generate a nucleic sequence up to the first stop codon in the right frame
def generate_random_prot( id_num, limit_sequence_size=20):
	sequence = str()
	alphabet = ['A','T','G','C']
	while len(sequence)<limit_sequence_size:
		new_codon = str() 
		nucleic_seq = str()
		while new_codon not in ['TGA', 'TAA', 'TAG']:
		
			new_codon = ''.join(random.choices(alphabet,k=3))
			nucleic_seq=nucleic_seq+new_codon
		nucleic_seq = Seq(nucleic_seq)
		sequence = nucleic_seq.translate(stop_symbol='')

	record = SeqRecord(sequence, id='RandProt_'+str(id_num))
	return record


def generate_random_prot_from_aa_freq(id_num, frequencies,limit_sequence_size=20):
	choices = list(frequencies.keys())
	weight  = list(frequencies.values())
	sequence = 'M'
	while len(sequence)< limit_sequence_size:
		sequence = 'M'
		while sequence[-1]!='*':
			seq_letter = choice(choices, None,
							  p=weight)
			sequence+=seq_letter
	sequence = sequence[:-1]
	record = SeqRecord(Seq(sequence), id='RandProt_'+str(id_num))
	return record


def get_aa_frequencies(records):
	count = Counter()
	totcount = 0
	for rec in records:
		seq = rec.seq+'*'
		totcount+=len(seq)
		count.update(seq)
	return {x:y/totcount for x,y in count.items()}

def make_fragmented_sequence(sequence, min_size=0.1, max_size=0.9, side=None):
	seqlen = len(sequence)-1
	if not side:
		side = random.choice(['N', 'C'])
	if side == 'N':
		cut = random.randrange(int(seqlen*min_size),int(seqlen*max_size)+1)
		fragment_seq = sequence[:cut+1]
	elif side=='C':
		cut = random.randrange(seqlen-int(seqlen*max_size),seqlen-int(seqlen*min_size)+1)

		fragment_seq = sequence[cut:]
	return fragment_seq



def generate_simulations(indir, outdir, contdir, reuse=False):
	contamination_numbers = [10,20,50,100,200,500,1000]
	completeness_proportions = [0.1,.2,.3,.4,.5,.6,.7,.8,.9]
	added_bits_proportions =   [0.1,.2,.3,.4,.5,.6,.7,.8,.9]
	added_mimic_bits_proportions = [0.1,.2,.3,.4,.5,.6,.7,.8,.9]
	fragmented_proportion =  [0.1,.2,.3,.4,.5,.6,.7,.8,.9] 
	fused_proportion =  [0.1,.2,.3,.4,.5,.6,.7,.8,.9] 
	if not os.path.isdir(os.path.join(outdir,'Complete')):
		os.mkdir(os.path.join(outdir,'Complete'))
	if not os.path.isdir(os.path.join(outdir,'Storage')):
		os.mkdir(os.path.join(outdir,'Storage'))
	files = list(os.listdir(indir))
	cont_files = list(os.listdir(contdir))

	for file in files:
		path = os.path.join(indir,file)
		target_proteome = load_proteome(path)

		for p in completeness_proportions:
			outf = 'Comp'+str(p)+'_'+file
			opath = os.path.join(outdir,'Complete',outf)
			spath = os.path.join(outdir,'Storage',outf)
			if reuse and os.path.exists(spath):
				partial_proteome = load_proteome(spath)
			else:
				partial_proteome = draw_partial_proteome(target_proteome, p)
				write_proteome(spath, partial_proteome)
			write_proteome(opath, partial_proteome)

		for p in added_bits_proportions:
			outf = 'Error'+str(p)+'_'+file
			opath = os.path.join(outdir,'Complete',outf)
			spath = os.path.join(outdir,'Storage',outf)
			if reuse and os.path.exists(spath):
				false_seq = load_proteome(spath)
				erroneous_proteome = target_proteome + false_seq
			else:
				
				erroneous_proteome, false_seq = add_random_sequences(target_proteome, p)
				write_proteome(spath, false_seq)

			write_proteome(opath, erroneous_proteome)

		for p in added_mimic_bits_proportions:
			outf = 'RealisticError'+str(p)+'_'+file
			opath = os.path.join(outdir,'Complete',outf)
			spath = os.path.join(outdir,'Storage',outf)
			if reuse and os.path.exists(spath):
				false_seq = load_proteome(spath)
				erroneous_proteome = target_proteome + false_seq
			else:
				frequencies = get_aa_frequencies(target_proteome)
				erroneous_proteome, false_seq = add_random_sequences(target_proteome, p, frequencies)
				write_proteome(spath, false_seq)
			write_proteome(opath, erroneous_proteome)

		for p in fragmented_proportion:
			outf = 'Fragment'+str(p)+'_'+file
			opath = os.path.join(outdir,'Complete',outf)
			spath = os.path.join(outdir,'Storage',outf)
			if reuse and os.path.exists(spath):
				fragments = load_proteome(spath)
				fragmented_proteome = regenerate_fragmented(target_proteome, fragments)
			else:
				fragmented_proteome, fragments = create_fragments(target_proteome, p)
				write_proteome(spath, fragments)
			write_proteome(opath, fragmented_proteome)

		for p in fused_proportion:
			outf = 'Fused'+str(p)+'_'+file
			opath = os.path.join(outdir,'Complete',outf)
			spath = os.path.join(outdir,'Storage',outf)
			if reuse and os.path.exists(spath):
				fused_proteins = load_proteome(spath)
				fused_proteome = regenerate_fusion(target_proteome, fused_proteins)
			else:
				fused_proteome, fused_proteins = create_fusions(target_proteome, p)
				write_proteome(spath, fused_proteins)
			write_proteome(opath, fused_proteome)

		for other_file in cont_files:
			cpath = os.path.join(contdir, other_file)
			contaminant_proteome = load_proteome(cpath)
			for n in contamination_numbers:				
				outf = 'Cont'+str(n)+'_'+file+'_'+other_file
				opath = os.path.join(outdir,'Complete',outf)
				spath = os.path.join(outdir,'Storage',outf)	
				if reuse and os.path.exists(spath):
					contamination = load_proteome(spath)
					contaminated_proteome = target_proteome + contamination
				else:
					contaminated_proteome, contamination = add_contamination(target_proteome, contaminant_proteome,n)
					write_proteome(spath, contamination)

				write_proteome(opath, contaminated_proteome)
if __name__=='__main__':

	parser = build_arg_parser()  
	arg = parser.parse_args()
	input_folder = arg.input
	output_folder = arg.output
	cont_folder = arg.contamnants
	reuse  = arg.reuse
	generate_simulations(input_folder, output_folder, cont_folder, reuse)


