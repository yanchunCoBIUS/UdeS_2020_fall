import RNA
import sys
import os
import pandas as pd
from Bio import SeqIO
from io import StringIO


def read_fasta_file(input_fasta_file, num_of_sequences):
    """
        This function allows to read RNA sequences from a fasta file

        Parameters
        ----------

        input_fasta_file
            RNA sequence file path

        Returns
        -------
        id_seq, sequence, dict_RNA


    """
    dict_RNA = {}
    counter = 0
    with open(input_fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # read sequences
            if counter >= num_of_sequences:
                break
            dict_RNA[record.id] = {}
            # dict_RNA[record.id]['id'] = record.id
            dict_RNA[record.id]['sequence'] = str(record.seq)
            counter = counter + 1
    return dict_RNA


def build_subopt_result(structure, energy, data):
    """
        Build the subopt result dictionary
    """
    if not structure == None:
        # use the same number of structures as num_of_clusters 
        if data['counter'] < data['num_of_clusters']:
            # print (">subopt %d" % data['counter'])
            # print ("%s" % data['sequence'])
            # print ("%s [%6.2f]" % (structure, energy))

            data['structure'][">subopt %d" % data['counter']] = structure

            # increase structure counter
            data['counter'] = data['counter'] + 1
        else:
            pass


def compute_structure(family_name, num_of_clusters, num_of_sequences):
    dict_RNA = read_fasta_file("Input/" + family_name + ".fa/" + family_name + ".fa", num_of_sequences)
    
    # seq = "ATAATGATACTTCCGTCCAGTCACAGTCCGAGCGTGAAGCGGCAAGCCTGCCGCAATCCGCCGCAGGCAATGAGGGCGGCCCGGTGAGATCCACGGGGTGGTGAAATTCCAATGA"

    # Set global switch for unique ML decomposition
    RNA.cvar.uniq_ML = 1

    for id in dict_RNA:
        seq = dict_RNA[id]['sequence']
        
        subopt_data = { 'counter' : 0, 'structure': dict_RNA[id], 'num_of_clusters': num_of_clusters}

        # Create a 'fold_compound' for our sequence
        a = RNA.fold_compound(seq)

        # Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
        # the MFE and print each structure using the function above
        a.subopt_cb(500, build_subopt_result, subopt_data)
        
    return dict_RNA

