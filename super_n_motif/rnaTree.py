import scipy.cluster.hierarchy as hier
import numpy as np
import copy
import json

class RnaTree:
    """ Tree based clustering of RNA based on secondary RNA structure

    Parameters
    ----------
    rna_sequence: json object
        Contain the RNA sequences

    Attributes
    ----------
    tree: array
        Array containing the RNA tree
    """
    def __init__(self, rna_sequence, num_iteration):
        self.rna_sequence = rna_sequence
        self.num_iteration = num_iteration
        self.tree = []

    def build(self, dissim_matrix, SSbyNm_matrix, **kwargs):

        """
        Note
        ----
        More options could be added to the tree.build() function. Options could be:
         - Which library to use and the method to calculate the distance
           as 'euclidian' or 'ward'. We use scipy as default, with the 'single' method.
         - Which score to use to measure the motif importance for node. We use variance and percentage as default.

        Parameters
        ----------
        dissim_matrix: 2D array float
            Array of the dissimilary between each RNA sequence
        SSbyNm_matrix: 2D array int
            Array
        **kwargs : dictionary
            Various options to change the tree construction. See note

        """
        # Linkage option build from the dissimilarity matrix. Use scipy and single algorithm
        link = hier.linkage(dissim_matrix, 'single')

        # Build tree
        # Initiate the leaves

        # sequence_count = len(self.rna_sequence)
        sequence_count = 20

        motif_list = SSbyNm_matrix[0][1:]

        occurance_matrix = self.extract_motif_occurrence(SSbyNm_matrix)

        for k in range(sequence_count):
            rna_sequence_object = self.rna_sequence.iloc[:, k]
            # We ignore every motif that has no occurrence for a given sequence
            ind = (occurance_matrix[k]).nonzero()

            # List of base motif that the super-n-program can output
            base_motif_list = {'B':[0, {}], 'H':[0, {}], 'I':[0, {}], 'M':[0, {}], 'P':[0, {}], 'S':[0, {}], 'E3':[0, {}], 'E5':[0, {}], 'G4':[0, {}]}

            # Iterate over each leaf and organize their motif list in a more lisible way for the json
            for i in range(len(ind[0])):
                current_motif = motif_list[ind[0][i]]
                current_motif_occurance = occurance_matrix[k][ind[0][i]]

                if current_motif in base_motif_list:
                    base_motif_list[current_motif] = [int(current_motif_occurance), {}]
                else:
                    name = {current_motif : int(current_motif_occurance)}
                    base_motif_curr = current_motif[0]
                    # Check if the base motif is E or G and complete them
                    if base_motif_curr == 'E' or base_motif_curr == 'G':
                        base_motif_curr = current_motif[:2]
                    base_motif_list[base_motif_curr][1].update(name)

            leaf = { "name": rna_sequence_object.name, "id": k + 1, "isLeaf": True, "Model": base_motif_list,
                   "sequence": rna_sequence_object['sequence'], "structure":rna_sequence_object['>subopt %d' % self.num_iteration], "children":[] }
            self.tree.append(leaf)

        # Build the node of the tree
        for pairing in link:
            # Calculate the motif occurance in node
            dico = {}
            ind1, ind2 = int(pairing[0]), int(pairing[1])
            children = self.find_child([ind1, ind2])

            # Find all motif associated with children to calculate their presence in percentage
            for child in children:
                for key in self.tree[child]['Model']:
                    if self.tree[child]['Model'][key][1]:
                        for subKey in self.tree[child]['Model'][key][1]:
                            if subKey in dico:
                                dico[subKey] += 1
                            else:
                                dico[subKey] = 1

            # Which score to use to calculate the motif importance
            if "score" in kwargs:
                if kwargs["score"] == "variance":
                    motif_var = self.score_variance(children)

            # Calculate the presence of a motif by its percentage
            print (children)
            nb_child = len(children)
            for key in dico:
                dico[key] = dico[key] * 100 / nb_child
            
            # Remove motif which are under a chosen percentage
            dico = {k: v for k, v in dico.items() if v  >= int(kwargs['percentage'])}

            node = {"name": motif_var[0][0]+"_"+str(motif_var[0][1]), "id": len(self.tree) + 1, "isLeaf": False, "score": {'Percentage': dico, 'Variance': motif_var},
                    "sequence": "", "structure": "", "children": [ind1, ind2]}
            self.tree.append(node)

    def extract_motif_occurrence(self, matrix):
        """
        Parameters
        ----------
        matrix: 2D array of string
            Contain the super-n-motif matrix extracted from matNmRepRaw_SSbyNm.csv
            The header is composed of every motif found in the given RNA sequence list
            The first column is every RNA sequence
            The rest of the matrix is the number of occurrence for a specific motif for each RNA sequence

        Returns
        -------
        SSbyNm: int matrix
            RMatrix with only the number of occurrence for each motif
        """

        # Remove header
        matrix = np.delete(matrix, 0, 0)
        SSbyNm = []
        for x in matrix:
            line = x[1:len(x) - 1]
            SSbyNm.append(line)

        # Convert to numpy int array
        SSbyNm = np.array(SSbyNm).astype(int)
        return SSbyNm

    def score_variance(self, children):
        """

        Parameters
        ----------
        children: array of indices
            The incices correspond to specific children of a node in the tree

        Returns
        -------
        variance: float array
            Array of ordered variance for each base motif
        """
        motif_var = {'B':0, 'H': 0, 'I': 0, 'M': 0, 'P': 0, 'S': 0,'E3':0, 'E5': 0, 'G4': 0}
        motif_to_remove = []
        for key in motif_var:
            motif_occurance = []
            for child in children:
                motif_occurance.append(self.tree[child]['Model'][key][0])
            if np.count_nonzero(motif_occurance) == 0:
                motif_to_remove.append(key)
            else:
                motif_var[key] = np.around(np.var(np.array(motif_occurance)), 2)
        for motif in motif_to_remove:
            motif_var.pop(motif)
        motif_var = sorted(motif_var.items(), key=lambda tup: tup[1])
        return motif_var

    def find_child(self, index_array):
        """ Return the children for an index in a tree

        Parameters
        ----------
        index_array: int
            Array of two indices that are either the children of a node

        Returns
        -------
        array: int array
            All the children for a given array of indices
        """

        children = []
        # search children for left index
        if (not self.tree[index_array[0]]['isLeaf']):
            children = self.find_child(self.tree[index_array[0]]['children'])
        else:
            children = [index_array[0]]

        # search children for right index
        if (not self.tree[index_array[1]]['isLeaf']):
            children = children + self.find_child(self.tree[index_array[1]]['children'])
        else:
            children = children + [index_array[1]]
        return children

    def tree_format(self, format, *filename):
        """
        Note
        ----
        Only json format is supported, but it could be extented to support newick or xml.

        Parameters
        ----------
        format: str
            Indicate the format the tree should take.
        filename: str, optional
            Indicate the name of the save file. If empty, no file is saved.

        Returns
        -------

        """

        format_str = ""
        if format == 'json':
            format_str = self.format_json_recurs(len(self.tree)-1)
            if filename:
                with open(filename[0], 'w') as outfile:
                    json.dump(format_str, outfile)
        return format_str

    def format_json_recurs(self, ind):
        """ Recursively iterate the tree and convert it to a string of the chosen format

        Parameters
        ----------
        ind: int
            Tree node or leaf indice

        Returns
        -------
        data: string
            String with the part of the tree it finished copying
        """
        if not self.tree[ind]['isLeaf']:
            ind1, ind2 = self.tree[ind]['children']
            str1 = self.format_json_recurs(ind1)
            str2 = self.format_json_recurs(ind2)
            data = copy.deepcopy(self.tree[ind])
            data['children'] = [str1, str2]
        else:
            data = copy.deepcopy(self.tree[ind])
        return data