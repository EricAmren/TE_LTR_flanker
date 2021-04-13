#!/usr/bin/python3
## Dictionnary generator

## Dependencies
import argparse
import logging

## Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("consensus_fasta_file", help="The consensus TE fasta file", type=str)
parser.add_argument("--suffix", default="standard_suffixes.tsv", help="A tsv file enumerating each possible suffix for the intern part (in first line) and LTR part (in second line) of a TE.", type=str)

args = parser.parse_args()

def generate_matching_pairs_from_suffixes(consensus_fasta, suffix_file):
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

    with open(consensus_fasta, 'r') as input:
        for line in input:
            if line.startswith('>'):
                seq_ID = line.split()[-1].replace('>', '')
                TE_name = remove_suffix(seq_ID)
                if is_LTR_part(seq_ID):
                    LTR_part_dict[TE_name] = seq_ID
                elif is_I_part(seq_ID):
                    I_part_dict[TE_name] = seq_ID
                else:
                    not_splitted_TE.append(seq_ID)
    
    matching_TE_parts = set(I_part_dict.keys()).intersection(set(LTR_part_dict.keys()))
    unmatched_I_list = [I_part_dict[x] for x in set(I_part_dict.keys()).difference(set(LTR_part_dict.keys()))]
    unmatched_LTR_list = [LTR_part_dict[x] for x in set(LTR_part_dict.keys()).difference(set(I_part_dict.keys()))]
    if len(unmatched_I_list + unmatched_LTR_list) > 0 :
        logging.warning("Some TE names contain a recognized suffix, yet couldn't be matched. While this can be expected, you might want to manually cure some of them in the dictionnary :\n" + ', '.join(unmatched_I_list) + "\n" + ', '.join(unmatched_LTR_list) + "\n")
    TE_dictionnary = ""
    for TE_name in matching_TE_parts:
        TE_dictionnary += "\t".join([TE_name, I_part_dict[TE_name], LTR_part_dict[TE_name]]) + "\n"
    TE_dictionnary += '\n'.join(unmatched_I_list + unmatched_LTR_list + not_splitted_TE)
    return TE_dictionnary

if __name__ == "__main__":
    print(generate_matching_pairs_from_suffixes(args.consensus_fasta_file, args.suffix))
