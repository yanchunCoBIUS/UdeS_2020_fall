from super_n_motif.project import *
from subopt import *
import json


def generate_cluster(family_list, num_of_structures, num_of_sequences):
    input_for_SNM = ""

    for family_name in family_list:
        dict_RNA, rna_sequence = compute_structure(family_name, num_of_structures, num_of_sequences)
        input_for_SNM = input_for_SNM + '\n' + rna_sequence
            
    runSNM(input_for_SNM)
