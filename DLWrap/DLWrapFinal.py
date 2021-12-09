from MLPScreenFinal import MLPScreen
from MLPScreenFinal import MLPScreenWrap
from MLPPrediction import MLPPrediction
from MLPPrediction import MLPPredictionWrap
from AssoPred import AssoPred
from myGBLUP import myGBLUP

from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import roc_auc_score
#import cPickle
import time
import sys
import numpy
import argparse
from numpy import loadtxt
from pathlib import Path
from itertools import zip_longest
import pickle
import os
import pandas
import multiprocessing

numpy.set_printoptions(precision=4, linewidth=200)

if __name__ == '__main__':
    numpy.random.RandomState(1234567890)
    numpy.random.seed(6780)
    parser = argparse.ArgumentParser()
    numThreads = len(os.sched_getaffinity(0))
    #numThreads = 20
    print("The number of existing threads is "+str(numThreads))
    
   
    
    ###################################################################
    ##################### Procedure parameters ########################
    ###################################################################
    parser.add_argument("--TestSelection", metavar="TestSelection", type=int, help="Analysis to be performed.=1: Perform screening analysis for one gene; =2: Perform screening analysis for multiple genes; =3: Perform prediction analysis given screening results; =4: Perform both screening and prediction analyses; and =5: Perform prediction with selected genes", required=True, choices=[1,2,3,4,5])
                        
    ######################################################################
    ##################### Data-related parameters ########################
    ######################################################################
    
    ### Phenotypes ##
    parser.add_argument('--train_y_input', metavar='train_y_input', default=None, help='Phenotype of training set: FID, IID, Y and it should have header', required=True)
    parser.add_argument('--test_y_input', metavar='test_y_input', default=None, help='Phenotype of testing set. It is in the same format as the training set')
    parser.add_argument('--binary_outcome', metavar='binary_outcome', default=0, type=int, help='Whether the outcome is binary: =1 yes binary; =0 no continuous. Default is continuous (i.e., =0)')
    
    ###### Genotypes ######
    ## only for one gene ##
    parser.add_argument('--train_x_input', metavar='train_x_input', default=None, help='Genotype of training set,the comma delimited version with the frist 6 columns being FID, IID, PAT, MAT, SEX, PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header. Valid when only one gene is to be screened.')
    parser.add_argument('--gene_name', metavar="gene_name",default=None, help="The name of the gene. If not provided, it will be 0. Only valid when screening_index=1. ")
    parser.add_argument('--test_groups', metavar='test_groups', default=None, help='Groups of tests that need to be performed. If not provided, then all genes will be grouped together. Valid when only one gene is to be screened.')
    
    ## multiple genes ##
    parser.add_argument('--GeneticDataTrain', metavar='GeneticDataTrain', default=None, help='The list of training genotype files. Each file contains genotypes of training set, that is the comma delimited version with the frist 6 columns being FID,IID,PAT,MAT,SEX,PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header.')
    parser.add_argument('--GeneticDataTest', metavar='GeneticDataTest', default=None, help='The list of testing genotype files. Each file contains genotypes of testing set, that is the comma delimited version with the frist 6 columns being FID,IID,PAT,MAT,SEX,PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header. The files should be ordered according to the same order in the GeneticDataTrain (i.e., in the same gene order)')
    parser.add_argument('--geneindexFile', metavar='geneindexFile', default=None, help='The list of gene_names, and it should be in the same order as the GeneticDataTrain folder.')
    parser.add_argument('--test_groupsFile', metavar='test_groupsFile', default=None, help='The list of groups of test files. It should be in the same order as the GeneticDataTrain folder. Each file contains the grouping for each gene.-----NOT SUPPORTED FOR THE MOMENT!')  
    
    ######################################################################
    ##################### Output-related parameters ########################
    ######################################################################
    parser.add_argument('--AssociationDir', metavar='AssociationDir', default=None, help='The folder where screening results are saved.')
    parser.add_argument('--outputPredFile', metavar='outputPredFile', default=None, help='The predicted value output')
    
    ##################################################################
    #################### modules for MLP parameters ##################
    ##################################################################
    # model parameters #
    parser.add_argument('--seed_value', metavar='seed_value', type=int, default=0, help='The random_seed to be set. Default 0. ')  
    parser.add_argument('--nunit1', metavar='nunit1', type=int, default=50, help='The number of hidden units in the first hidden layer. The default is 50. Valid only when screening is performed.')
    parser.add_argument('--nunit2', metavar='nunit2', type=int, default=10, help='The number of hidden units in the second hidden layer. The default is 10. Valid only when screening is performed.')
    parser.add_argument('--adjusted', metavar='adjusted', type=bool, default=True, help='Whether the number of hidden units should be adjusted, default is True. Valid only when screening is performed.')
    parser.add_argument('--vfold', metavar='vfold', type=int, default=10, help='The number of cross-validation. Default is 10. Valid only when screening is performed.')
    parser.add_argument('--nperm', metavar='nperm', type=int, default=100, help='The number of permutation. The default is 100. Valid only when screening is performed.')    
    parser.add_argument('--reg_tune', metavar='reg_tune', type=int, default=1, help='Whether the dropout percentage should be tuned. 1=Yes. 0=No. The default is 1. Valid only when screening is performed.')
    parser.add_argument('--reg_set', metavar='reg_set', type=float, default=0.2, help='The dropout percentage. The default is 0.2. If the dropout percentage is to be tuned, then this parameter is ignored. Valid only when screening is performed.')      
    parser.add_argument('--output_level',metavar='output_level',type=int,default=0,help='The level of details of the output. Default is 0. =0:  only useful output is saved. For screening, the p-value for each gene and the model are saved. =1: intermediate level output is saved. For screening, the p-value for each gene, the model, the score are saved. Valid only when screening is performed.')
    
    #####################################################################
    ############# modules for prediction-related parameters #############
    #####################################################################   
    parser.add_argument('--alpha', metavar='alpha', type=float, default=0.05, help='the significance level for selecting predictive genes. Valid only prediction is performed.')   
    parser.add_argument('--pre_selected_gene_List', metavar='pre_selected_gene_List', default=None, help='The list of pre-selected genes')
    parser.add_argument('--pre_selected_model_list', metavar='pre_selected_model_list', default=None, help='The list of models locations for the pre-selected genes.')   
    args = parser.parse_args()
    
    if args.TestSelection not in [1,2,3,4,5]:
        raise("Set the augument of TestSelection to select the tests that you wish to perform")

    #####################################################################
    ################ Screening for one gene procedure ###################
    #####################################################################
    if (args.TestSelection==1):
        if args.AssociationDir is None:
            raise("The output folder (AssociationDir) for screening results must be specified!")
        if args.train_x_input is None:
            raise("The genotype file (train_x_input) must be specifed.")
        MLPScreenWrap(outputfolder=args.AssociationDir, train_x_input=args.train_x_input, train_y_input=args.train_y_input, gene_name=args.gene_name, binary_outcome=args.binary_outcome, vfold=args.vfold, nperm=args.nperm, reg_tune=args.reg_tune, reg_set=args.reg_set, nunit1=args.nunit1, nunit2=args.nunit2, seed_value=args.seed_value, group_index_file=args.test_groups, numThreads=numThreads, output_level=args.output_level, adjusted=args.adjusted)
    

    #####################################################################
    ####################  Screening multiple genes  #####################
    #####################################################################
    if args.TestSelection==2:
        if args.GeneticDataTrain is None:
            raise("A file (GeneticDataTrain) containing all genotypes files must be specified!")
        if args.AssociationDir is None:
            raise("The output folder (AssociationDir) for screening results must be specified!")
        if args.geneindexFile is None:
            raise("A file (geneindexFile) containing all gene names must be specified!")          
        predasso=AssoPred()
        print(args.GeneticDataTrain)
        predasso.AssoPredFunc(GeneticDataTrain=args.GeneticDataTrain, GeneticDataTest=args.GeneticDataTest, AssociationDir=args.AssociationDir, train_y_input=args.train_y_input, geneindexFile=args.geneindexFile, outputPredFile=args.outputPredFile, numThreads=numThreads, binary_outcome=args.binary_outcome, test_y_input=args.test_y_input, alpha=args.alpha, vfold=args.vfold, nperm=args.nperm, reg_tune=args.reg_tune, reg_set=args.reg_set, nunit1=args.nunit1, nunit2=args.nunit2, seed_value=args.seed_value, output_level=args.output_level, adjusted=args.adjusted, assc_only=True)
        
    #####################################################################
    ########################  Prediction only  ##########################
    #####################################################################
    if (args.TestSelection==3):
        if args.GeneticDataTrain is None:
            raise("A file (GeneticDataTrain) containing all genotypes files must be specified!")
        if args.GeneticDataTest is None:
            raise("A file (GeneticDataTest) containing all testing genotypes files must be specified!")
        if args.AssociationDir is None:
            raise("The folder (AssociationDir) containing all screening results must be specified!")
        if args.geneindexFile is None:
            raise("A file (geneindexFile) containing all gene names must be specified!")  
        if args.outputPredFile is None:
            raise("A file (outputPredFile) to saved predicted value must be specified!")  
        pred=MLPPredictionWrap(seed_value=args.seed_value, numThreads=numThreads, binary_outcome=args.binary_outcome, GeneticDataTrain=args.GeneticDataTrain, GeneticDataTest=args.GeneticDataTest, AssociationDir=args.AssociationDir, geneindexFile=args.geneindexFile, train_y_input=args.train_y_input, outputPredFile=args.outputPredFile, alpha=args.alpha, test_y_input=args.test_y_input)
        
    
    #####################################################################
    ####################  Associaiton + Prediction  #####################
    #####################################################################
    if args.TestSelection==4:
        if args.GeneticDataTrain is None:
            raise("A file (GeneticDataTrain) containing all genotypes files must be specified!")
        if args.GeneticDataTest is None:
            raise("A file (GeneticDataTest) containing all testing genotypes files must be specified!")
        if args.AssociationDir is None:
            raise("The output folder (AssociationDir) for screening results must be specified!")
        if args.geneindexFile is None:
            raise("A file (geneindexFile) containing all gene names must be specified!")  
        if args.outputPredFile is None:
            raise("A file (outputPredFile) to saved predicted value must be specified!")  
        predasso=AssoPred()
        predasso.AssoPredFunc(GeneticDataTrain=args.GeneticDataTrain, GeneticDataTest=args.GeneticDataTest, AssociationDir=args.AssociationDir, train_y_input=args.train_y_input, geneindexFile=args.geneindexFile, outputPredFile=args.outputPredFile, numThreads=numThreads, binary_outcome=args.binary_outcome, test_y_input=args.test_y_input, alpha=args.alpha, vfold=args.vfold, nperm=args.nperm, reg_tune=args.reg_tune, reg_set=args.reg_set, nunit1=args.nunit1, nunit2=args.nunit2, seed_value=args.seed_value, output_level=args.output_level, adjusted=args.adjusted, assc_only=False)
        
    if args.TestSelection==5:
        if args.pre_selected_gene_List is None:
            raise("A file containing all selected genes must provided!")
        if args.pre_selected_model_list is None:
            raise("A file containing all models of the selected genes must provided!")


