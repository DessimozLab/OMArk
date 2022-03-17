from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
import random


def build_arg_parser():
    """Handle the parameter sent when executing the script from the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Generate simulated proteomes from an input directory of proteomes.")   
    parser.add_argument('-i', '--input', help="The directory with input genomes." )	
    parser.add_argument('-o', '--output', help="The directory with output genomes.")
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

def add_random_sequences(target_prot, proportion):
	number_genes  = int(len(target_prot)*proportion)
	added_prot = [generate_random_prot(x) for x in range(number_genes)]
	expanded_prot = target_prot + added_prot
	return expanded_prot, added_prot

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

def generate_simulations(indir, outdir, reuse=False):
	contamination_numbers = [10,20,50,100,200,500,1000]
	completeness_proportions = [0.1,.2,.3,.4,.5,.6,.7,.8,.9]
	added_bits_proportions =   [0.1,.2,.3,.4,.5,.6,.7,.8,.9]
	if not os.path.isdir(os.path.join(outdir,'Complete')):
		os.mkdir(os.path.join(outdir,'Complete'))
	if not os.path.isdir(os.path.join(outdir,'Storage')):
		os.mkdir(os.path.join(outdir,'Storage'))
	files = list(os.listdir(indir))

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
				erronoeous_proteome = target_proteome + false_seq
			else:
				
				erronoeous_proteome, false_seq = add_random_sequences(target_proteome, p)
				write_proteome(spath, false_seq)

			write_proteome(opath, erronoeous_proteome)
		for other_file in files:
			cpath = os.path.join(indir, other_file)
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
	reuse  = arg.reuse
	generate_simulations(input_folder, output_folder, reuse)


