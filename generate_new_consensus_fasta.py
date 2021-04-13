#!/usr/bin/python3

## Generating new fasta with flanked intern parts using dictionnary

## Dependencies
import argparse
import logging
from Bio import SeqIO

## Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="The consensus TE fasta file", type=str)
parser.add_argument("-d", "--dict", help="TSV file generated using `generate_dictionnary.py`. See README for more information.", type=str)
parser.add_argument("-o", "--output", help="File where the new fasta will be outputted", type=str)

args = parser.parse_args()


def import_dictionnary(TE_dictionnary_file):
    TE_dict = dict()
    with open(TE_dictionnary_file, 'r') as dictionnary:
        for line in dictionnary :
            if len(line.split()) == 3:
                TE_name, I_part, LTR_part = line.split()
                TE_dict[TE_name] = [I_part, LTR_part]
            elif len(line.split()) == 1:
                TE_name = line.strip()
                TE_dict[TE_name] = []
            else :
                raise ValueError("Invalid dictionnary line : \n" + line + "\nEach line should contains either 3 columns corresponding to the TE_name, TE_intern_part, TE_LTR_part OR only contains 1 column corresponding to the TE_name.")
    return TE_dict

def get_TE_ID(record):
    return record.split()[-1]

def generating_LTR_flanked_fasta_file(consensus_fasta, TE_dictionnary_file, output_file):
    TE_dict = import_dictionnary(TE_dictionnary_file)
    I_dict = {v[0] : k for k, v in TE_dict.items() if v}
    LTR_dict = {v[1] : k for k, v in TE_dict.items() if v}

    record_iterator = SeqIO.parse(consensus_fasta, "fasta")
    reconstructed_records = dict()
    new_records = []
    for record in record_iterator:
        TE_ID = get_TE_ID(record.description)
        if TE_ID in TE_dict: # If TE_ID directly in keys, then it means its a non_splitted TE
            record.id = TE_ID
            record.description = ""
            # new_records.append(record)
        elif TE_ID in I_dict:
            if I_dict[TE_ID] not in reconstructed_records :
                reconstructed_records[I_dict[TE_ID]] = [record, None]
            else :
                reconstructed_records[I_dict[TE_ID]][0] = record
        elif TE_ID in LTR_dict:
            if LTR_dict[TE_ID] not in reconstructed_records :
                reconstructed_records[LTR_dict[TE_ID]] = [None, record]
            else :
                reconstructed_records[LTR_dict[TE_ID]][1] = record
        else:
            raise ValueError
    for merged_TE_ID, TE_parts_list in reconstructed_records.items() :
        I_part_record, LTR_part_record = TE_parts_list
        new_record = I_part_record
        new_record.seq = LTR_part_record.seq + I_part_record.seq + LTR_part_record.seq
        new_record.id = merged_TE_ID
        new_record.description = ""
        new_records.append(new_record)
    SeqIO.write(new_records, output_file, "fasta")

if __name__ == "__main__":
    generating_LTR_flanked_fasta_file(args.fasta, args.dict, args.output)
