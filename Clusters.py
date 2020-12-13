import multiprocessing
from super_n_motif.super_n_motif import *
from subopt import *
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import numpy as np
import os
import uuid
from sklearn.cluster import SpectralClustering, MeanShift, KMeans, Birch, AffinityPropagation, AgglomerativeClustering, DBSCAN, OPTICS, MiniBatchKMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import *


class Clusters:
    """
    Parameters
    ----------
    
    Attributes
    ----------
    dataFrameSNM:
    family_list:  
    num_of_structures: 
    num_of_sequences:
    """
    def __init__(self, family_list, num_of_structures, num_of_sequences):
        self.family_list = family_list
        self.num_of_structures = num_of_structures
        self.num_of_sequences = num_of_sequences
        self.dict_allFamilies = dict()
        self.dataFrameSNM = self.transform_SNM_to_DataFrame()
        self.score_list = []                

    def get_score_list(self):
        return self.score_list


    def splitTrainAndTest(self, trainSize, testSize):
        """
        Create two dataFrames:
        :param trainSize: training size
        :param testSize:  test size
        :return: split data
        """
        x, y = self.getFeaturesAndTarget()

        # Create two new dataframes, one with the training rows, one with the test rows
        xTrain, xTest, yTrain, yTest = train_test_split(x, y, train_size=trainSize, test_size=testSize,
                                                        shuffle=True, stratify=y)
        return xTrain, xTest, yTrain, yTest

    def getFeaturesAndTarget(self):
        # features
        x = self.dataFrameSNM.iloc[:, 1:-1]
        # target
        y = self.dataFrameSNM.iloc[:, [-1]]
        return x, y
    
    def use_allFamilies_runSNM(self, uid):
        """
        Run SUM
        """
        input_for_SNM = ""

        for family_name in self.family_list:
            dict_RNA, rna_sequence = compute_structure(family_name, self.num_of_structures, self.num_of_sequences)

            # dictionary contains all information
            self.dict_allFamilies[family_name] = dict_RNA

            input_for_SNM = input_for_SNM + '\n' + rna_sequence
                
        runSNM(input_for_SNM, uid)
        
        # pool = multiprocessing.Pool()
        # results = pool.map(self.run_subopt, self.family_list)
        # pool.close()
        # pool.join()
        
        # input_for_SNM = ""
        # for i in results:
        #     input_for_SNM = input_for_SNM + i
            
        # runSNM(input_for_SNM, uid)
    
    
    # def run_subopt(self, family_name):
    #     dict_RNA, rna_sequence = compute_structure(family_name, self.num_of_structures, self.num_of_sequences)

    #     # dictionary contains all information
    #     self.dict_allFamilies[family_name] = dict_RNA

    #     return rna_sequence


    def get_dataFrameSNM(self):
        """
        return SNM's DataFrame
        """
        return self.dataFrameSNM

    def get_dict_allFamilies(self):
        """
        return dictionary contains all information
        """
        return self.dict_allFamilies

    def transform_SNM_to_DataFrame(self):
        """
        Read Super_n_motif(SNM) result: matSnmRep_SSbySnm.csv as a DataFrame
        """
        uid = "family{}_structure{}_sequence{}".format(len(self.family_list), self.num_of_structures, self.num_of_sequences)
        self.use_allFamilies_runSNM(uid)

        df = pd.read_csv(os.path.join(os.getcwd(), "super_n_motif", "lib", "superNMotif", uid+"_resultMatrix", "matSnmRep_SSbySnm.csv"))
        df['Family'] = [self.family_list.index(i.split('#')[0]) for i in df.iloc[:, 0]]
        # family_column = []
        # for i in range(len(self.family_list)):
        #     single_family = np.full(self.num_of_structures * self.num_of_sequences, i)
        #     family_column.extend(single_family)
        # try:
        #     df['Family'] = family_column
        # except AttributeError as e:
        #     pass
        return df

    def k_means(self, x, y):
        model = KMeans(n_clusters=len(self.family_list), random_state=0)
        
        label = model.fit_predict(x)
        # print("Kmeans")
       
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label


    def agglomerativeClustering(self, x, y):
        """
        linkage: {“ward”, “complete”, “average”, “single”}, default=”ward”
            ward: minimizes the variance of the clusters being merged.
            average: uses the average of the distances of each observation of the two sets.
            complete:  uses the maximum distances between all observations of the two sets.
            single uses the minimum of the distances between all observations of the two sets.
        affinity: str or callable, default=’euclidean’
        """
        model = AgglomerativeClustering(n_clusters=len(self.family_list))
        label = model.fit_predict(x)
        # print("Agglomerative")
        
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)
        return label


    
    def birch(self, x, y):
        """
        """
        model = Birch(n_clusters=len(self.family_list))
        label = model.fit_predict(x)
        # print("Birch")
        
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label


    def dbscan(self, x, y):
        """
        metric: string, or callable, default=’euclidean’
            The metric to use when calculating distance between instances in a feature array.
        epsfloat: default=0.5
        min_samplesint: default=5
        """
        model = DBSCAN()
        label = model.fit_predict(x)
        # print("DBSCAN")
        
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label

    
    def miniBatchKMeans(self, x, y):
        """
        """
        model = MiniBatchKMeans(n_clusters=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        # print("MiniBatchKMeans")
        
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label


    def meanShift(self, x, y):
        """
        """
        model = MeanShift()
        label = model.fit_predict(x)
        # print("MeanShift")
        
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label


    def spectralClustering(self, x, y):
        """
        """
        model = SpectralClustering(n_clusters=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        # print("SpectralClustering")
       
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)
        return label

    
    def gaussianMixture(self, x, y):
        """
        """
        model = GaussianMixture(n_components=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        # print("GaussianMixture")
       
        score = adjusted_rand_score(np.array(y['Family']), label)
        self.score_list.append(score)
        # print("\t ARI MEASURE: %0.3f" % score)

        return label
