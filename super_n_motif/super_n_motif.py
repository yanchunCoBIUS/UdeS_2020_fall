import os
import subprocess
import numpy as np
import csv
import shutil
import platform
import uuid


def runSNM(dot_bracket_str, uuid):
    """ Run the super-n-motif exec.

    Reconstruct the RNA sequence file that super-n-motif program take as input.
    The super-n-motif program output a dissimilarity matrix and a matrix with the occurrence of motifs for
    each RNA sequence.

    """
    if platform.system() == 'Windows':
        path_to_prog = os.path.join(os.getcwd(), "super_n_motif", "lib",  "superNMotifWin")
    else:
        path_to_prog = os.path.join(os.getcwd(), "super_n_motif", "lib", "superNMotif")

    #Save to file
    filename = os.path.join(path_to_prog, uuid+'_RNA.db')
    # if not os.path.isdir(path_to_prog):
    #     os.mkdir(path_to_prog)
    
    with open(filename, 'w') as file:
        file.write(dot_bracket_str)

    arg_list = [os.path.join(path_to_prog, "supernmotifs"), "-p", "1", "-i", filename, "-o", os.path.join(path_to_prog, uuid + "_resultMatrix")]

    subprocess.check_call(arg_list, cwd=path_to_prog)
    
    # remove file of rna sequence once done with it
    # os.remove(filename)

