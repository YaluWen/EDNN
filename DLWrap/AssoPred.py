#!/usr/bin/env python
import multiprocessing
from MLPScreenFinal import MLPScreen
from MLPScreenFinal import MLPScreenParallel
from MLPScreenFinal import MLPScreenWrap
from MLPPrediction import MLPPrediction
from MLPPrediction import MLPPredictionWrap
import os
import numpy
import pandas
import pathlib
from functools import partial

class AssoPred:
    def __init__(self, verbose=False):
        self.verbose = verbose;
        pass
    
    def AssoPredFunc(self, GeneticDataTrain, GeneticDataTest, AssociationDir, train_y_input, geneindexFile, outputPredFile, numThreads, binary_outcome, test_y_input=None, alpha=0.1, vfold=5, nperm=20, reg_tune=1, reg_set=0.2, nunit1=50, nunit2=10, seed_value=0, output_level=1, adjusted=True, assc_only=False):
        pathlib.Path(AssociationDir).mkdir(parents=True, exist_ok=True) 
        genes=numpy.genfromtxt(geneindexFile,dtype=str)
        traindata=numpy.genfromtxt(GeneticDataTrain,dtype=str)
        if assc_only==False:
            testdata=numpy.genfromtxt(GeneticDataTest,dtype=str)
        if genes.shape[0]!=traindata.shape[0]:
            print("Genetic Training Data and the genes are not consistent!")
            raise SystemExit(0)
        if assc_only==False:
            if testdata.shape[0]!=traindata.shape[0]:
                print("Genetic Training Data and Testing Data are not consistent!")
                raise SystemExit(0)
        
        # load ytrain #
        tmpdata = pandas.read_csv(train_y_input,sep='\t|,',engine='python')
        tmpdata = tmpdata.drop(tmpdata.columns[[0, 1]], axis=1);
        yall=tmpdata.to_numpy()
        inputs=numpy.column_stack((traindata,genes))
        
        # create a "pool" of 3 processes to do the calculations
        pool = multiprocessing.Pool(processes=numThreads//2)
        numThreadspara=numThreads//(numThreads//2)
        param_fixed={"outputfolder":AssociationDir, "train_y_input":train_y_input, "binary_outcome":binary_outcome,"vfold":vfold, "nperm":nperm, "reg_tune":reg_tune,"reg_set":reg_set,"nunit1":nunit1,"nunit2":nunit2, "seed_value":seed_value, "numThreads":numThreadspara,"output_level":output_level,"adjusted":adjusted}
        pool_fixed=partial(MLPScreenParallel, fixed_param=param_fixed)
        pool.map(pool_fixed, inputs)
        pool.close()
        pool.join()

        if assc_only==False:
            pred=MLPPredictionWrap(seed_value=seed_value, numThreads=numThreads, binary_outcome=binary_outcome,GeneticDataTrain=GeneticDataTrain, GeneticDataTest=GeneticDataTest, AssociationDir=AssociationDir, geneindexFile=geneindexFile, train_y_input=train_y_input, outputPredFile=outputPredFile, alpha=alpha)

        if test_y_input is not None:
            tmp=MLPPrediction()
            result=tmp.evaluate(pred, ytrue=None, ytruefile=test_y_input, ifbinary=binary_outcome)
            numpy.savetxt(outputPredFile+"summary",result)



