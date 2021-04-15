## What is this script good for ?
Consensus files of TE sequences often need a tedious pre-processing : LTR transposable elements are splitted between their intern and LTR parts. If you want to use those TE sequences for an alignment, you have to restore the real sequence of these TEs by flanking the intern part with the LTR part (like this : LTR|Intern|LTR ).

The goal of this script is to simplify this task and to provide a correctly reconstructed new fasta file.

## How to use it ?

This script need two files : a classic **consensus fasta file** (an example being the Dfam families.fa file present in the folder) and a **dictionnary in tsv format** providing correspondance between intern and LTR parts.

1) Generating the dictionnary

You can choose to write the dictionnary yourself (TE_name, TE_intern_part, TE_LTR_part, tab_separated):

dictionnary example :  
`Copia  Copia_I    Copia_LTR`  
`Roo   Roo-I_DM   Roo-LTR_DM`

OR generate it automatically using 'generate_dictionnary.py' script and the default list of suffixes ("standard_suffixes.tsv"):

`./generate_dictionnary.py -f families.fa -o families.dictionnary.tsv`

You can also specify a custom list of suffixes in a tsv file, with the first line being the list of intern suffixes, and the second line being the list of LTR suffixes, separated with tabs.

`./generate_dictionnary.py -f families.fa --suffix custom_suffixes.tsv -o families.dictionnary.tsv`

example of suffix tsv file :  
`_I -I_DM`  
`_LTR   -LTR_DM`

**It is advised to have a look at the generated dictionnary and manually cure it if needed.**

2) Get the new consensus fasta file using this command:

`./generate_new_consensus_fasta.py -f families.fa -d families.dictionnary -o new_families.fa`

Notes :

All TE sequences that are present in the fasta must be represented in the dictionnary.
TODO : Maybe add the option to directly output TE sequences that are not written in the dictionnary.
This can be easily implemented and tested by removing a random TE from the dictionnary...

