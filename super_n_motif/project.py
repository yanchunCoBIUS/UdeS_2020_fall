import os
import subprocess
import numpy as np
import csv
import shutil
from super_n_motif.rnaTree import RnaTree
import platform

class Project:
    """

    Parameters
    ----------
    args[0]: json object
        Containing the RNA sequences
    args[1]:  str
        user id

    Attributes
    ----------
    rna_sequence: json object
        Containing the RNA sequences
    tree: rnaTree object
        class containing the RNA sequence as a hierarchical clustering (tree based clustering )
    user_id: str
        id to avoid two or more user accessing the same data
    """
    def __init__(self, rna_sequence, user_id, num_iteration, family_name):
        self.rna_sequence = rna_sequence
        self.userId = user_id
        self.num_iteration = num_iteration
        self.family_name = family_name
        self.tree = RnaTree(self.rna_sequence, self.num_iteration)

    def runSNM(self):
        """ Run the super-n-motif exec.

        Reconstruct the RNA sequence file that super-n-motif program take as input.
        The super-n-motif program output a dissimilarity matrix and a matrix with the occurrence of motifs for
        each RNA sequence.

        """
        # create a dot bracket format string
        dot_bracket_str = ""
        for i in range(0, 20):
            rna_sequence_object = self.rna_sequence.iloc[:, i]
            dot_bracket_str = dot_bracket_str + ('>' + rna_sequence_object.name + '\n' +
                            rna_sequence_object['sequence'] + '\n' +
                            rna_sequence_object['>subopt %d' % self.num_iteration] + '\n'
                            )
        # for x in self.rna_sequence:
        #     dot_bracket_str = dot_bracket_str + ('>' + x['header'] + '\n' +
        #                     x['sequence'] + '\n' +
        #                     x['structure'] + '\n'
        #                     )
        if platform.system() == 'Windows':
            path_to_prog = os.path.join(os.getcwd(), "super_n_motif", "lib",  "superNMotifWin")
        else:
            path_to_prog = os.path.join(os.getcwd(), "super_n_motif", "lib", "superNMotif")
        #Save to file
        dir_path = os.path.join(path_to_prog, self.family_name)
        filename = os.path.join(dir_path, 'rna'+ str(self.userId) +'.db')
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
        with open(filename, 'w') as file:
            file.write(dot_bracket_str)
        FNULL = open(os.devnull, 'w')

        arg_list = [os.path.join(path_to_prog, "supernmotifs"), "-p", "5", "-i", filename, "-o", os.path.join(path_to_prog, self.family_name, "resultMatrix"+ str(self.userId))]

        subprocess.check_call(arg_list, cwd=path_to_prog)
        
        # remove file of rna sequence once done with it
        # os.remove(filename)


    def buildTree(self, **kwargs):
        """

        Parameters
        ----------
        option: array of parameters
            may contain various options to modify the tree structure. See the note in RnaTree class

        Returns
        -------
        tree: rnaTree class object

        """

        # Get the dissimilarity matrix between every RNA sequence
        dissim_matrix = self.readDissimMat()

        # Get the motif occurrence matrix
        if platform.system() == 'Windows':
            path_SNM = os.path.join(os.getcwd(), "super_n_motif", 'lib', 'superNMotifWin')
        else:
            path_SNM = os.path.join(os.getcwd(), "super_n_motif", 'lib', 'superNMotif')
        path_matrix = os.path.join(path_SNM, self.family_name, 'resultMatrix' + str(self.userId), 'matNmRepRaw_SSbyNm.csv')
        
        # change 'rb' to 'r'
        SSbyNm_matrix = np.array(list(csv.reader(open(path_matrix, "r"), delimiter=",")))

        self.tree.build(dissim_matrix, SSbyNm_matrix, **kwargs)

        #Delete file to avoid bloating server
        # dir_path = os.path.join(path_SNM, "resultMatrix" + str(self.userId))
        # shutil.rmtree(dir_path)

    def readDissimMat(self):

        """ Read a csv file of a dissimilarity matrix

        Note
        ----
        The output by super-n-motif program is a non standard csv. This mean it can't be read correctly by library like
        pandas or numpy. We read the csv line by line and convert it to a matrix.

        Returns
        -------
        matrix: float matrix
            The dissimilarity matrix for a given RNA sequence stored as a csv

        """
        if platform.system() == 'Windows':
            path_SNM = os.path.join(os.getcwd(), "super_n_motif", 'lib', 'superNMotifWin')
        else:
            path_SNM = os.path.join(os.getcwd(), "super_n_motif", 'lib', 'superNMotif')

        path_matrix = path_matrix = os.path.join(path_SNM, self.family_name, 'resultMatrix' + str(self.userId), 'matDissim_SSbySS.csv')
        
        # change 'rb' to 'r'
        reader = csv.reader(open(path_matrix, "r"), delimiter=",", quoting=csv.QUOTE_NONNUMERIC)

        dissim_matrix = np.array([], dtype=float)
        #skip first line of file which is empty
        next(reader, None)
        for x in reader:
            #skip last element in array which is empty. Convert array to numpy array
            line = x[0:len(x) - 1]
            dissim_matrix= np.append(dissim_matrix, line)

        return dissim_matrix
