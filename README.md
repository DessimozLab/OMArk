# OMArk

OMArk is a software for proteome (protein-coding gene repertoire) quality assessment. It provides measures of proteome completeness, characterizes the consistency of all protein coding genes with regard to their homologs, and identifies the presence of contamination from other species.
OMArk relies on the OMA orthology database, from which it exploits orthology relationships, and on the OMAmer software for fast placement of all proteins into gene families.

## Installation

You can use OMArk by installing the package through conda:

``conda install -c bioconda omark``

Alternatively, it can also be installed through Pypi:

``pip install omark``

Or by cloning this repository and installing it manually with your Python installer.  

Example command from the git directory:
``pip install .``


You can then use it on your Python environment by calling it as a command line tool.

#### OMAmer Database


OMArk relies on an OMAmer database to run. You can download one from the latest release of the OMA database on the ["Current release"](https://omabrowser.org/oma/current/) page of the OMA Browser.  
 
For all OMArk features to work correctly, it is recommended that this database covers a wide range of species. Thus we recommend using one constructed from the whole OMA database, often called [**LUCA.h5**](https://omabrowser.org/All/LUCA.h5) .   
Using a database for a more restricted taxonomic range (Metazoa, Viridiplantae, Primates) would limit the ability of OMArk to detect contamination or to identify sequences of species that belong outside this range.  

Alternatively, an OMAmer database is available through: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10034236.svg)](https://doi.org/10.5281/zenodo.10034236) - File :OMAmerDB.gz.  
This is the LUCA.h5 database constructed from the December 2021 release of the OMA database and is the one that was used for the OMArk [preprint](https://www.biorxiv.org/content/10.1101/2022.11.25.517970v1).


## Usage

Required arguments: ``-f (--file)``, ``-d (--database)``

    usage: omark [-h] (-f FILE | -c) -d DATABASE [-o OUTPUTFOLDER] [-t TAXID] [-of OG_FASTA] [-i ISOFORM_FILE] [-r TAXONOMIC_RANK]  [-v]


## Arguments
| Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [``-f`` ``--file``](#markdown-header--file)                  |                 |Path to an OMAmer search output file (Default mode)                                                                                                          |
| [``-c`` ``--output_cHOGs``](#markdown-header--output_cHOGs)  | False           |Switch OMArk mode to only computing a list of conserved HOGs and output it as list. Can be used to obtain a set of genes on which to train models.           |
| [``-d`` ``--db``](#markdown-header--database)                |                 | Path to an OMAmer database                                                                                                                                  |
| [``-o`` ``--outputFolder``](#markdown-header--outputFolder)  | ./omark_output/ | Path to the folder into which OMArk results will be output. OMArk will create it if it does not exist.                                                      |
| [``-t`` ``--taxid``](#markdown-header--taxid)                | None            | NCBI taxid corresponding to the input proteome (Optional).                                                                                                  |
| [``-of`` ``--og_fasta``](#markdown-header--og_fasta)         | None            | The original proteomes file. Provide if you want optional FASTA file to be outputted by OMArk (Sequences by categories, sequences by detected species, etc) |
| [``-i``, ``--isoform_file``](#markdown-header--isoform_file) | None            | A text file, listing all isoforms of each gene as semi-colon separated values, with one gene per line. Use if your input proteome include more than one protein per gene. See the [Splicing isoforms](#splicing-isoforms) section.|
| [``-r`` ``--taxonomic-rank``](#markdown-header--taxonomic-rank)| None           |The taxonomic rank (genus, order, family...) that should be used as ancestral lineage if possible.                                                           |
| [``-e`` ``--ete_ncbi_db``](#markdown-header--ete_ncbi_db)    | None			 |Path to the ete3 NCBI database to be used. Default will use the default location at ~/.etetoolkit/taxa.sqlite.                                               |
| [``-v`` ``--verbose``](#markdown-header--verbose)            | False           | Turn on logging information about OMArk progress.                                                                                                           |


## Input data

The standard input for the OMArk pipeline is a proteome - a FASTA file where each gene is represented by only one protein sequence.  
As described in the [Example](#example) section below, the first step of the pipeline is to run the OMAmer software on this FASTA file in order to obtain an OMAmer search result file.  
This OMAmer search file will be the main input of the OMArk software itself.

#### Splicing isoforms

If your proteome file contains multiple isoforms per gene, you can still use it as an input from the OMArk pipeline but it will require an additonal step.  

As before, you can run an OMAmer search on the whole proteome.  
When running OMArk, however, you must provide it *with a* ``.splice`` file via the ``--isoform_file`` option.  

A splice file is a text file where:
* Items on a line are protein identifiers (the same as in the associated FASTA file) separated by semi-colons (;).
* All proteins on the same line are products of the same gene.
* There are as many lines as there are genes in the proteome, with genes with only one product being represented by a single protein identifier.

Here is an extract of .splice file generated for Danio rerio RefSeq proteome, recapitulating protein isoforms of three genes:
```
NP_001258730.1;XP_005166105.1;XP_017211994.1;XP_009300826.1;XP_017211995.1
NP_001334620.1
NP_001300751.1;NP_571866.2;XP_005166949.1
```
OMArk will choose, for each gene, the isoform with the best OMAmer mapping as a representative for computing its metrics.

## Output

A default OMAmer output consists of 4 files with the same name but different extensions.

OMArk outputs the main results of the analysis in two complementary files: a machine-readable one, identified by its .sum extension, and a human-readable summary ending with ``_detailed_summary.txt``.
These commented files reports:
* The reference lineage that was used for quality assessment
* The number of conserved Hierarchical Orthologous Groups (HOGs) used for completeness assessment
* The completeness assessment results (Single, Duplicated, Missing)
* The whole proteome quality assessment results (Consistent placements, Inconsistent Placements, Contaminants, Missing genes)
* The species and contaminants detected in the proteome

The file with the .pdf extension is a graphical representation of the completeness and whole proteome quality assesment.

The file with the .tax extension indicates: the closest taxonomic lineage in the OMA database and the selected reference lineage.

The file with the .omq extension recapitulates the identifiers of the HOGs used in the completeness analysis, and the category to which they were attributed.

The file with the .ump extensions recapitulates the identifiers for all proteins by the category in which they were placed.

## Example

You can run OMArk on an example files set stored inside the example\_data folder. Remember to download an OMAmer database as indicated in the installation section.

First: you can run OMAmer on the proteome FASTA. (For more documentation about installing OMAmer: see its [Github](https://github.com/DessimozLab/omamer).
This step should take less than 15 minutes.

	omamer search --db  LUCA.h5 --query example_data/UP000005640_9606.fasta --out example_data/UP000005640_9606.omamer

Then, use OMArk (Should take less than 10 minutes) after creating an empty output folder:

	mkdir example_data/omark_output

	omark -f example_data/UP000005640_9606.omamer -d LUCA.h5 -o example_data/omark_output

You can now explore OMArk results in the ``omark_output`` folder.

## Visualising multiple OMArk results

We include a script for visualising OMArk results of multiple datasets. This script is available at the utils folder of the repository under the name: ``plot_all_results.py``.  

You can also use an interactive version of this script as a Jupyter Notebook (``Explore_data.ipynb``). Note the Notebook needs the ``plot_all_results.py`` as dependency and should be downloaded alongside it.

You can use the plotting script with the following command:


	plot_all_results.py -i example_data/omark_output -o fig.png
 
Note that by default it will use the prefix of the .sum file present in the OMArk folders as label. You can override this behaviour and provide more data by providing a mapping file (-m) option formatted as the ``mapping_template.txt`` file in the same folder.  
There you can choose to provide a Species name and a NCBI identifier for each result (Note: The filename column must be equal to the prefix of the .sum file for the mapping to be successful). 
If you do so and provide taxonomic information as part of the additional data, you can also choose to order the data according to the NCBI taxonomy using the ``-t`` option. 

## Webserver

Omark is available as a public webserver at <https://omark.omabrowser.org/home/> where users are free to upload a proteome and run the OMArk pipeline. OMArk results can then be viewed side-to-side with precomputed results on closely related species, as is the recommended use case for OMArk. Precomputed results available on the webserver are based on UniProt Reference Proteomes.  
