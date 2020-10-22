from super_n_motif.super_n_motif import *
from subopt import *
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import numpy as np
import os


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
                
    # def knn(self, xNewTrain, yNewTrain):
    #     """

    #     :param xNewTrain: features
    #     :param yNewTrain: target
    #     :return: trained classifier
    #     """
    #     # Use kNN classifier to predict
    #     # divide training into two dataFrames, one for training, one for find best key
    #     xKNNTrain, xKNNTest, yKNNTrain, yKNNTest = self.splitTrainAndTest(xNewTrain, yNewTrain, 0.2, 0.8)
    #     Classifier = findBestKforKNN(xKNNTrain, yKNNTrain)
    #     # Classifier = KNeighborsClassifier(n_neighbors=3)
    #     Classifier.fit(xKNNTest, yKNNTest)
    #     return Classifier

    # def rf(self, xNewTrain, yNewTrain):
    #     """

    #     :param xNewTrain: features
    #     :param yNewTrain: target
    #     :return: trained classifier
    #     """
    #     # Use Random Forest classifier to predict
    #     # k-fold cross validation score
    #     # Classifier = kFoldCrossValidation(xNewTrain, yNewTrain, randomForest(xNewTrain, yNewTrain))[0]
    #     Classifier = RandomForestClassifier(random_state=0)
    #     start_time = time.time()
    #     Classifier.fit(xNewTrain, yNewTrain)
    #     rfTrainingTime = (time.time() - start_time)
    #     # print("\n---Training RF %s seconds ---\n" % rfTrainingTime)
    #     return Classifier, rfTrainingTime

    # def nn(self, xNewTrain, yNewTrain):
    #     """

    #     :param xNewTrain: features
    #     :param yNewTrain: target
    #     :return: trained classifier
    #     """
    #     # Create a Neural Network Classifier
    #     # Classifier = neuralNetwork(xNewTrain, yNewTrain)
    #     Classifier = MLPClassifier(random_state=0)
    #     start_time = time.time()
    #     Classifier.fit(xNewTrain, yNewTrain)
    #     nnTrainingTime = (time.time() - start_time)
    #     # print("\n---Training NN %s seconds ---\n" % nnTrainingTime)
    #     return Classifier, nnTrainingTime

    # def adaboost(self, xNewTrain, yNewTrain):
    #     """

    #     :param xNewTrain: features
    #     :param yNewTrain: target
    #     :return: trained classifier
    #     """
    #     # Create a AdaBoost Classifier
    #     Classifier = AdaBoostClassifier()
    #     Classifier.fit(xNewTrain, yNewTrain)
    #     return Classifier

    def splitTrainAndTest(self, trainSize, testSize):
        """
        Create two dataFrames:
        :param trainSize: training size
        :param testSize:  test size
        :return: split data
        """
        # features
        x = self.dataFrameSNM.iloc[:, 1:-1]
        # target
        y = self.dataFrameSNM.iloc[:, [-1]]
        # Create two new dataframes, one with the training rows, one with the test rows
        xTrain, xTest, yTrain, yTest = train_test_split(x, y, train_size=trainSize, test_size=testSize,
                                                        shuffle=True, stratify=y)
        return xTrain, xTest, yTrain, yTest
    
    def use_allFamilies_runSNM(self):
        """
        Run SUM
        """
        input_for_SNM = ""

        for family_name in self.family_list:
            dict_RNA, rna_sequence = compute_structure(family_name, self.num_of_structures, self.num_of_sequences)

            # dictionary contains all information
            self.dict_allFamilies[family_name] = dict_RNA

            input_for_SNM = input_for_SNM + '\n' + rna_sequence
                
        runSNM(input_for_SNM)
        

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
        self.use_allFamilies_runSNM()

        df = pd.read_csv(os.path.join(os.getcwd(), "super_n_motif", "lib", "superNMotif", "resultMatrix", "matSnmRep_SSbySnm.csv"))
        
        family_column = []
        for i in range(len(self.family_list)):
            single_family = np.full(self.num_of_structures * self.num_of_sequences, i)
            family_column.extend(single_family)
        df['Family'] = family_column
        return df