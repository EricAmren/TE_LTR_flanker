#!/usr/bin/python3
## Dictionnary generator

## Dependencies
import argparse
import logging
from Bio import SeqIO

## Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="The consensus TE fasta file", type=str)
parser.add_argument("--suffix", default="standard_suffixes.tsv", help="A tsv file enumerating each possible suffix for the intern part (in first line) and LTR part (in second line) of a TE.", type=str)
parser.add_argument("-o", "--output", help="File where the dictionnary will be outputted", type=str)

args = parser.parse_args()

def get_TE_ID(record):
    return record.split()[-1]

def generate_matching_pairs_from_suffixes(consensus_fasta, suffix_file, output_file):
    with open(suffix_file, 'r') as input:
        I_suffix_list = input.readline().split()
        LTR_suffix_list = input.readline().split()
    I_suffix_set = set(I_suffix_list)
    LTR_suffix_set = set(LTR_suffix_list)

    def is_LTR_part(name):
        for suffix in LTR_suffix_set:
            if name.endswith(suffix):
                return True
        return False

    def is_I_part(name):
        for suffix in I_suffix_set:
            if name.endswith(suffix):
                return True
        return False

    def remove_suffix(seq_ID):
        suffix_list = list(LTR_suffix_set) + list(I_suffix_set)
        for suffix in suffix_list:
            if seq_ID.endswith(suffix):
                return seq_ID[:-len(suffix)]
        return seq_ID

    I_part_dict = {}
    LTR_part_dict = {}
    not_splitted_TE = []
    detected_conflicts = []

    records = list(SeqIO.parse(consensus_fasta, "fasta"))
    for record in records:
        TE_ID = get_TE_ID(record.description)
        TE_name = remove_suffix(TE_ID)
        if is_LTR_part(TE_ID):
            if TE_name in LTR_part_dict :
                detected_conflicts.append(TE_name)
                continue
            LTR_part_dict[TE_name] = TE_ID
        elif is_I_part(TE_ID):
            if TE_name in I_part_dict :
                detected_conflicts.append(TE_name)
                continue
            I_part_dict[TE_name] = TE_ID
        else:
            not_splitted_TE.append(TE_ID)

    if detected_conflicts :
        logging.warning("Detected conflict, following TEs can be matched with several suffixes : " + ', '.join(set(detected_conflicts)) + "\nYou should manually set it in the dictionnary.")
     
    matching_TE_parts = set(I_part_dict.keys()).intersection(set(LTR_part_dict.keys()))
    unmatched_I_list = [I_part_dict[x] for x in set(I_part_dict.keys()).difference(set(LTR_part_dict.keys()))]
    unmatched_LTR_list = [LTR_part_dict[x] for x in set(LTR_part_dict.keys()).difference(set(I_part_dict.keys()))]
    if len(unmatched_I_list + unmatched_LTR_list) > 0 :
        logging.warning("Some TE names contain a recognized suffix, yet couldn't be matched. While this can be expected, you might want to manually cure some of them in the dictionnary :\n" + ', '.join(unmatched_I_list) + "\n" + ', '.join(unmatched_LTR_list) + "\n")
    TE_dictionnary = ""
    for TE_name in matching_TE_parts:
        TE_dictionnary += "\t".join([TE_name, I_part_dict[TE_name], LTR_part_dict[TE_name]]) + "\n"
    TE_dictionnary += '\n'.join(unmatched_I_list + unmatched_LTR_list + not_splitted_TE)
    with open(output_file, 'w') as output:
        output.write(TE_dictionnary)
    # return TE_dictionnary
    

if __name__ == "__main__":
    generate_matching_pairs_from_suffixes(args.fasta, args.suffix, args.output)
