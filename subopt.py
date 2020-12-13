import RNA
import sys
import os
import pandas as pd
from Bio import SeqIO
from io import StringIO


def build_subopt_result(structure, energy, data):
    """
        Build the subopt result dictionary
    """
    dot_bracket_str = ""
    
    if not structure == None:
        # use the same number of structures as num_of_structures 
        if data['counter'] < data['num_of_structures']:
            # print (">subopt %d" % data['counter'])
            # print ("%s" % data['sequence'])
            # print ("%s [%6.2f]" % (structure, energy))
            id_subopt = ">%s#%s#subopt%d" % (data['family_name'], data['id'], data['counter'])
            data['rna_dict']['%s' % id_subopt] = structure

            dot_bracket_str = dot_bracket_str + (id_subopt + '\n' +
                            data['rna_dict']['sequence'] + '\n' +
                            structure + '\n'
                            )

            data['rna_sequence'].append(dot_bracket_str)

            # increase structure counter
            data['counter'] = data['counter'] + 1
        else:
            pass


def compute_structure(family_name, num_of_structures, num_of_sequences):
    rna_sequence = []
    dict_RNA, rna_sequence = read_fasta_file(family_name, num_of_sequences, num_of_structures, rna_sequence)

    # seq = "ATAATGATACTTCCGTCCAGTCACAGTCCGAGCGTGAAGCGGCAAGCCTGCCGCAATCCGCCGCAGGCAATGAGGGCGGCCCGGTGAGATCCACGGGGTGGTGAAATTCCAATGA"

    # # Set global switch for unique ML decomposition
    # RNA.cvar.uniq_ML = 1

    # for id in dict_RNA:
    #     seq = dict_RNA[id]['sequence']
        
    #     subopt_data = { 'counter' : 0, 'rna_dict': dict_RNA[id], 'num_of_structures': num_of_structures, 'id': id, 'rna_sequence': rna_sequence, 'family_name': family_name}

    #     # Create a 'fold_compound' for our sequence
    #     a = RNA.fold_compound(seq)

    #     # Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
    #     # the MFE and print each structure using the function above
    #     a.subopt_cb(500, build_subopt_result, subopt_data)
        
    return dict_RNA, "\n".join(rna_sequence)


def read_fasta_file(family_name, num_of_sequences, num_of_structures, rna_sequence):
    """
        This function allows to read RNA sequences from a fasta file

        Parameters
        ----------

        family_name, num_of_sequences, num_of_structures

        Returns
        -------
        dict_RNA, rna_sequence


    """
    input_fasta_file = "Input/" + family_name + ".fa/" + family_name + ".fa"
    dict_RNA = {}
    counter = 0
    
    # Set global switch for unique ML decomposition
    RNA.cvar.uniq_ML = 1
    
    with open(input_fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # read sequences
            if counter >= num_of_sequences:
                break
            id = record.id
            dict_RNA[id] = {}
            dict_RNA[id]['sequence'] = str(record.seq)

            seq = str(record.seq)
            subopt_data = { 'counter' : 0, 'rna_dict': dict_RNA[id], 'num_of_structures': num_of_structures, 'id': id, 'rna_sequence': rna_sequence, 'family_name': family_name}
            # Create a 'fold_compound' for our sequence
            a = RNA.fold_compound(seq)
            # Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
            # the MFE and print each structure using the function above
            a.subopt_cb(150, build_subopt_result, subopt_data)
        
            counter = counter + 1
    return dict_RNA, rna_sequence