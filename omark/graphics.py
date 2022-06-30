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


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#Color code for every category:
F_SINGLE = '#9af19cff'
F_DUP = '#359b73ff'
MISSING = '#d55e00ff'

CONSISTENT = '#3db7e9ff'
INCONSISTENT = '#8a5df1ff'
CONTAMINATION = '#e69f00ff'
UNKNOWN = '#000000ff'


def plot_omark_results(format_filename_dict, results, results_proteomes, fragment_info = True):
	fig, axes = plt.subplots(2, 1, figsize = (3,10))
	total_completeness = len(results['Single'])+len(results['Duplicated'])+ len(results['Overspecific_S']) + len(results['Overspecific_D'])+ len(results['Underspecific']) + len(results['Lost'])
	single = 100*(len(results['Single'])+len(results['Overspecific_S'])+len(results['Underspecific']))/total_completeness
	dup = 100*(len(results['Duplicated']) + len(results['Overspecific_D']))/total_completeness
	miss = 100*len(results['Lost'])/total_completeness
	axes[0].set_ylim(0,100)
	axes[0].bar(x=['a'],height=single, label='Single', color= F_SINGLE, alpha=.99)
	axes[0].bar(x=['a'],height=dup, bottom = single,label='Duplicated', color= F_DUP, alpha=.99)
	axes[0].bar(x=['a'],height=miss, bottom = single+dup, label='Missing', color=MISSING, alpha=.99)
	axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=11)

	total_proteomes = len(results_proteomes['Not_Placed'])+len(results_proteomes['Correct'])+len(results_proteomes['Contamination'])+len(results_proteomes['Erroneous'])
	correct = 100*len(results_proteomes['Correct'])/total_proteomes
	cont = 100*len(results_proteomes['Contamination'])/total_proteomes
	incorrect = 100*len(results_proteomes['Erroneous'])/total_proteomes
	nomap = 100*len(results_proteomes['Not_Placed'])/total_proteomes
	axes[1].set_ylim(0,100)
	custom_legend = []
	if fragment_info:

		correct_partial = 100*len(results_proteomes['Correct_Partial'])/total_proteomes
		correct_fragment = 100*len(results_proteomes['Correct_Fragment'])/total_proteomes
		correct_full = correct-(correct_partial+correct_fragment)
		incorrect_partial = 100*len(results_proteomes['Erroneous_Partial'])/total_proteomes
		incorrect_fragment = 100*len(results_proteomes['Erroneous_Fragment'])/total_proteomes
		incorrect_full = incorrect-(incorrect_partial+incorrect_fragment)
		cont_partial = 100*len(results_proteomes['Contamination_Partial'])/total_proteomes
		cont_fragment = 100*len(results_proteomes['Contamination_Fragment'])/total_proteomes
		cont_full = cont-(cont_partial+cont_fragment)

		axes[1].bar(x=['a'],height=correct_full, label='Consistent', color= CONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=correct_partial, bottom = correct_full, hatch='///', color= CONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=correct_fragment, bottom=correct_full+correct_partial, hatch='\\\\\\', color= CONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=cont_full, bottom = correct,label='Contaminant', color= CONTAMINATION, alpha=.99)
		axes[1].bar(x=['a'],height=cont_partial, bottom = correct+cont_full, hatch='///', color= CONTAMINATION, alpha=.99)
		axes[1].bar(x=['a'],height=cont_fragment, bottom = correct+cont_full+cont_partial, hatch='\\\\\\', color= CONTAMINATION, alpha=.99)
		axes[1].bar(x=['a'],height=incorrect_full, bottom =correct+cont,label='Inconsistent', color= INCONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=incorrect_partial, bottom = correct+cont+incorrect_full, hatch='///', color= INCONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=incorrect_fragment, bottom = correct+cont+incorrect_full+incorrect_partial, hatch='\\\\\\', color= INCONSISTENT, alpha=.99)

		partial_patch = mpatches.Patch(facecolor='white', hatch='///', label='Partial mapping', alpha=.99)
		fragment_patch = mpatches.Patch(facecolor='white', hatch='\\\\\\', label='Fragments', alpha=.99)
		custom_legend = [partial_patch,fragment_patch]

	else:
		axes[1].bar(x=['a'],height=correct, label='Consistent', color= CONSISTENT, alpha=.99)
		axes[1].bar(x=['a'],height=cont, bottom = correct,label='Contaminant', color= CONTAMINATION, alpha=.99)
		axes[1].bar(x=['a'],height=incorrect, bottom = correct+cont, label='Inconsistent', color=INCONSISTENT, alpha=.99)
	axes[1].bar(x=['a'],height=nomap, bottom = correct+cont+incorrect, label='Unknown', color=UNKNOWN, alpha=.99)
	axes[1].legend(handles=axes[1].get_legend_handles_labels()[0]+custom_legend, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=11)

	xticks = axes[1].xaxis.get_major_ticks()
	axes[0].set_ylabel('Proportion of conserved HOGs', fontsize=11)
	axes[1].set_ylabel('Proportion of proteomes', fontsize=11)
	axes[0].get_xaxis().set_visible(False,)
	axes[1].get_xaxis().set_visible(False,)

	axes[1].set_ylim(axes[1].get_ylim()[::-1])

	plt.subplots_adjust(wspace=0, hspace=0)

	for fmt, filename in format_filename_dict.items():
		plt.savefig(filename, format=fmt, bbox_inches="tight")
