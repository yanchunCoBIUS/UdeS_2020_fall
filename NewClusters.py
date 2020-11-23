from super_n_motif.super_n_motif import *
from subopt import *
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import numpy as np
import os
from sklearn.cluster import SpectralClustering, MeanShift, KMeans, Birch, AffinityPropagation, AgglomerativeClustering, DBSCAN, OPTICS, MiniBatchKMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import accuracy_score


class NewClusters:
    """
    Parameters
    ----------
    
    Attributes
    ----------
    family_list:  
    num_of_structures: 
    num_of_sequences:
    """
    def __init__(self, family_list, num_of_structures, num_of_sequences):
        self.family_list = family_list
        self.num_of_structures = num_of_structures
        self.num_of_sequences = num_of_sequences
        # self.df = new_df
        self.dataFrameSNM = self.transform_SNM_to_DataFrame()
                            

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

    def get_dataFrameSNM(self):
        """
        return SNM's DataFrame
        """
        return self.dataFrameSNM
    

    def transform_SNM_to_DataFrame(self):
        """
        11/17 Need to improve: change name and class, because is reuse the Clusters class.
        """
        # df = self.df
        df = pd.read_csv(os.path.join(os.getcwd(), "super_n_motif", "lib", "superNMotif", "resultMatrix", "matSnmRep_SSbySnm.csv"))
        
        family_column = []
        for i in range(len(self.family_list)):
            single_family = np.full(self.num_of_sequences, i)
            family_column.extend(single_family)
        df['Family'] = family_column
        return df


    def k_means(self, x, y):
        model = KMeans(n_clusters=len(self.family_list), random_state=0)
        # yTrain_label = kmeans.labels_
        # print(yTrain_label)
        # print(np.array(yTrain['Family']))
        label = model.fit_predict(x)
        print("NEW: Kmeans")
        # print("Predict: \n", yTest_predict)
        # print("Actual: \n", np.array(yTest['Family']))

        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
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
        print("NEW: Agglomerative")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label


    
    def birch(self, x, y):
        """
        """
        model = Birch(n_clusters=len(self.family_list))
        label = model.fit_predict(x)
        print("NEW: Birch")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
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
        print("NEW: DBSCAN")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label

    
    def miniBatchKMeans(self, x, y):
        """
        """
        model = MiniBatchKMeans(n_clusters=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        print("NEW: MiniBatchKMeans")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label


    def meanShift(self, x, y):
        """
        """
        model = MeanShift()
        label = model.fit_predict(x)
        print("NEW: MeanShift")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label


    def spectralClustering(self, x, y):
        """
        """
        model = SpectralClustering(n_clusters=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        print("NEW: SpectralClustering")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label

    
    def gaussianMixture(self, x, y):
        """
        """
        model = GaussianMixture(n_components=len(self.family_list), random_state=0)
        label = model.fit_predict(x)
        print("NEW: GaussianMixture")
        score = accuracy_score(np.array(y['Family']), label)
        print("\t accuracy score:", score)
        return label
