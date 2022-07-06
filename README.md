# OMArk

OMArk is a software of proteome (protein-coding gene repertoire) quality assessment. It provides measure of proteome completeness, characterize all protein coding genes in the light of existing homologs, and identify the presence of contamination from other species.
OMArk rely on the OMA orthology database, from which it exploits orthology relationships, and on the OMAmer software for fast placement of all proteins into gene families.

## Installation

You can use OMArk by cloning this repository and installing it manually with your Python installer.

Example command from the git directory:
``python setup.py install``
or
``pip install .``

You can then use it on your Python environment by calling it as a command line tool.
OMArk rely on an OMAmer database to run. For all OMArk features to work correctly, it is better for this database to cover a wide range of species.
We recommend using one constructed from the whole OMA database. You can download one manually on this link : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6462027.svg)](https://doi.org/10.5281/zenodo.6462027)File :  OMAmerDB.tar.gz


##Usage

Required arguments: ``-f (--file)``, ``-d (--database)``

    usage :python omark.py [-h] [-f FILE] [-d DATABASE] [-o OUTPUTFOLDER] [-t TAXID] [-of OG_FASTA]

## Arguments
| Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [``-f`` ``--file``](#markdown-header--file)                  |                 | Path to an OMAmer search output file                                                                                                                        |
| [``-d`` ``--db``](#markdown-header--database)                |                 | Path to an OMAmer database                                                                                                                                  |
| [``-o`` ``--outputFolder``](#markdown-header--outputFolder)  | ./omark_output/ | Path to the folder into which OMArk results will be output. OMArk will create it if it does not exist.                                                      |
| [``-t`` ``--taxid``](#markdown-header--taxid)                | None            | NCBI taxid corresponding to the input proteome (Optional).                                                                                                  |
| [``-of`` ``--og_fasta``](#markdown-header--og_fasta)         | None            | The original proteomes file. Provide if you want optional FASTA file to be outputted by OMArk (Sequences by categories, sequences by detected species, etc) |
| [``-i``, ``--isoform_file``](#markdown-header--isoform_file) | None            | A semi-colon separated file, listing all isoforms of each genes, with one gene per line. Use if your input proteome include more than one protein per gene. |
| [``-v`` ``--verbose``](#markdown-header--verbose)            | False           | Turn on logging information about OMArk progress.                                                                                                           |

## Output

A default OMAmer output consists of 4 files with the same name but different extensions.

OMArk output the main results of the analysis in two complementary files: a machine-readable one, identified by its .sum extension, and a human-readable summary ending with ``_detailed_summary.txt``.
These commented files reports:
* The reference lineage that was used for quality assessment
* The number of conserved Hierarchical Orthologous Groups (HOGs) used for completeness assessment
* The completeness assessment results (Single, Duplicated, Missing)
* The whole proteome quality assessment results (Consistent placements, Inconsistent Placements, Contaminants, Missing genes)
* The species and contaminant detected in the proteome

The file with the .tax extension indicate: the closest taxonomic lineage in the OMA database and the selected reference lineage.

The file with the .omq extension recapitulates the HOGs identifier used in the completeness analysis, and the category to which they were attributed.

The file with the .ump extensions recapitulates the identifier for all proteins that were not mapped in OMAmer.

##Example

You can run OMArk on an example files stored on the example\_data folder. Remember to download an OMAmer databqse as indicated in the installation section.

First: you can run OMAmer on the proteome FASTA. (For more documentation about installing OMAmer: see its [Github](https://raw.githubusercontent.com/DessimozLab/omamer)
This step should take less than 15 minutes.

	omamer search --db  LUCA.h5 --query example_data/UP000005640_9606.fasta  --score sensitive --out example_data/UP000005640_9606.omamer

Then, use OMArk (Should take less than 10 minutes) after creating an empty output folder:

	mkdir example_data/omark_output

	omark -f example_data/UP000005640_9606.omamer -d LUCA.h5 -o example_data/omark_output

You can now explore OMArk results in the ``omark_output`` folder
