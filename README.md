# OMArk

OMArk is a software of proteome (protein-coding gene repertoire) quality assesment. It provides measure of proteome completeness, characterize all protein coding genes in the light of existing homologue,and identify the presence of contaminatino from other species.
OMArk rely on the OMA orthology database, from which it exploits orthology relationships, and on the OMAmer software for fast placement of all proteins into gene families.

## Installation

You can install OMArk by cloning this repertory and installing manually.

##Usage

Required arguments: ``-f (--file)``, ``-d (--database)``, ``-o (--output)`` 

    usage : omark [-h] [-f FILE] [-d DATABASE] [-o OUTPUTFOLDER] [-t TAXID] [-of OG_FASTA]

## Arguments
| Flag                 | Default                | Description |
|:--------------------|:----------------------|:-----------|
| [``-f`` ``--file``](#markdown-header--file)||Path to an OMAmer search output file
| [``-d`` ``--db``](#markdown-header--database)||Path to an OMAmer database
| [``-o`` ``--outputFolder``](#markdown-header--outputFolder)||Path to an (existing) folder into which output OMArk results.
| [``-t`` ``--taxid``](#markdown-header--taxid)|None| NCBI taxid corresponding to the input proteome (Optionnal).
| [``-of`` ``-og_fasta``](#markdown-header--og_fasta)|None| The original proteomes file. Provide if you want optionnal FASTA file to be outputted by Omark (Sequences by categories, sequences by detected species, etc)

## Output

A default OMAmer output consists of 4 files with the same name but different extensions.

OMArk output the main results of the analysis in a file with the .sum extension in the output folder
This commented file reports:
* The reference lineage that was used for quality assesment
* The number of conserved Hierarchical Orthologous Groups (HOGs) used for completeness assesment
* The completeness assesment results (Single, Duplicated, Missing)
* The whole proteome quality assesment results (Accurate placements, Erroneous Placements, Contaminants, Missing genes)
* The species and contaminant detected in the proteome

The file with the .tax extension indicate: the closest taxonomic lineage in the OMA database and the selected reference lineage.

The file with the .omq extension recapitules the HOGs identifier used in the compleness analysis, and the catefory to which they were attributed.

The file with the .ump extensions recapitules the identifier for all proteins that were not mapped in OMAmer.

