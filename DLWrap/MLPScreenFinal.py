#!/usr/bin/env python
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import roc_auc_score
from statistics import stdev
from statistics import mean
from numpy import loadtxt
from keras import layers
from keras import regularizers
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
#import matplotlib.pyplot as plt
import numpy
import scipy
from scipy.stats import wilcoxon
from numpy.lib.scimath import log
import tensorflow as tf
from pathlib import Path
import pandas
import pickle




import os
import random
#import tensorflow
def MLPScreenParallel(inputs, fixed_param):
    outputfolder=fixed_param["outputfolder"]
    train_y_input=fixed_param["train_y_input"]
    binary_outcome=fixed_param["binary_outcome"] 
    vfold=fixed_param["vfold"]
    nperm=fixed_param["nperm"] 
    reg_tune=fixed_param["reg_tune"]
    reg_set=fixed_param["reg_set"]
    nunit1=fixed_param["nunit1"]
    nunit2=fixed_param["nunit2"]
    seed_value=fixed_param["seed_value"]
    numThreads=fixed_param["numThreads"]
    output_level=fixed_param["output_level"]
    adjusted=fixed_param["adjusted"]
    train_x_input=inputs[0]
    gene_name=inputs[1]
    MLPScreenWrap(outputfolder=outputfolder, train_x_input=train_x_input, train_y_input=train_y_input, binary_outcome=binary_outcome, vfold=vfold, nperm=nperm, reg_tune=reg_tune, reg_set=reg_set, nunit1=nunit1, nunit2=nunit2, seed_value=seed_value, group_index_file=None, numThreads=numThreads,output_level=output_level, adjusted=adjusted, gene_name=gene_name)

    

def MLPScreenWrap(outputfolder, train_x_input, train_y_input, binary_outcome, vfold, nperm, reg_tune, reg_set, nunit1, nunit2, seed_value, group_index_file, numThreads,output_level, adjusted=False, gene_name=None):
    # raise exception #
    if (outputfolder is None):
        raise Exception('outputfolder must be supplied for screening')
    if (train_x_input is None):
        raise Exception('train_x_input must be supplied for screening')
    if (train_y_input is None):
        raise Exception('train_y_input must be supplied for screening')
    
    # load data for gene x to be screened#
    x = pandas.read_csv(train_x_input,sep='\t|,',engine='python')
    tmpdata = x.drop(x.columns[[0, 1, 2, 3, 4, 5]], axis=1);
    x=tmpdata.to_numpy()
    y = pandas.read_csv(train_y_input,sep='\t|,',engine='python')
    tmpdata =y.drop(y.columns[[0, 1]], axis=1);
    y=tmpdata.to_numpy()
    
    # adjust nunit1 and nunit2 if necessary #
    if adjusted:
        print("Adjusting nunits for each hidden layers")
        nunit1=min(x.shape[1]//2,nunit1);
        nunit2=min(x.shape[1]//4,nunit2); 
        nunit1=max(200,nunit1); #cat(nunit2,nunit1)
        nunit2=max(20,nunit2);
        print("nunit1="+str(nunit1)+";nunit2="+str(nunit2))

    if (gene_name is None):
        geneindex=0
    else:
        geneindex=gene_name
        
    print("screening for gene "+str(geneindex))
    
    outputfile=outputfolder+'/result'+geneindex+'.npy'
    outputfilemodel=outputfolder+'/result_model'+geneindex
        
    if binary_outcome==1:
        loss='binary_crossentropy'
        cross_cri='auc'
    else:
        loss='mse'
        cross_cri='cor'

    ifbinary=binary_outcome;

    ### Variable Screening part ###   
    if group_index_file is not None:
        if os.path.exists(group_index_file):
            txt = Path(group_index_file).read_text()
            txt = txt.replace(',', '\t')
            group_index = []
            for line in txt.splitlines():
                group_index.append([int(i) for i in line.split()])
        else:
            print('Groups of tests that need to be performed are not provided. All genetic variants will be grouped together! ')
            group_index = [[a+1 for a in range(x.shape[1])]]
            group_index.append([1])
            if (x.shape[1]>1):
                group_index.append([2])
    else:
        print('Groups of tests that need to be performed are not provided. All genetic variants will be grouped together! ')
        group_index= [[a+1 for a in range(x.shape[1])]]
        group_index.append([1])
        if (x.shape[1]>1):
            group_index.append([2])
    #print(group_index)
        
    ## Running screening models ##
    mlpscreen = MLPScreen(seed_value=seed_value,numThreads=numThreads)
    print("Screening: start model fitting for gene"+geneindex)
    model,pvalue,scorereturn,scorelistreturn,reg = mlpscreen.deeptestfunc(xdata=x, ydata=y, group_index=group_index, nperm=nperm, vfold=vfold, cross_cri=cross_cri, loss=loss, nunit1=nunit1, nunit2=nunit2, reg=reg_set, reg_tune=reg_tune,ifbinary=ifbinary,geneindex=geneindex)
    
    if(output_level==1):
        my_pvaluescore={'pvalue':pvalue, 'score':scorereturn,'scorelist':scorelistreturn}
    else:
        my_pvaluescore={'pvalue':pvalue}
    
    # Save all items in a tuple
    items_to_save = (my_pvaluescore)
    
    #print(outputfile)
    with open(outputfile, 'wb') as filehandler:
        pickle.dump(items_to_save, filehandler)

    # save model
    #print(outputfilemodel)
    model.save(outputfilemodel)


        
        

class MLPScreen:
    def __init__(self, verbose=False, seed_value=0,numThreads=-9):
        self.verbose = verbose
        if (seed_value!=0):
            # Seed value seed_value= 0
            # 1. Set the `PYTHONHASHSEED` environment variable at a fixed value
            os.environ['PYTHONHASHSEED']=str(seed_value)
        
            # 2. Set the `python` built-in pseudo-random generator at a fixed value
            random.seed(seed_value)

            # 3. Set the `numpy` pseudo-random generator at a fixed value
            numpy.random.seed(seed_value)
        
            # 4. Set the `tensorflow` pseudo-random generator at a fixed value
            tf.compat.v1.set_random_seed(seed_value) # later version
        
            # 5. Configure a new global `tensorflow` session
            inter_op_parallelism_threads=1
            if numThreads>0:
                numIntraOpThreads = numThreads // inter_op_parallelism_threads
            else:
                numIntraOpThreads = 4
            os.environ['OMP_NUM_THREADS'] = str(numIntraOpThreads)
            import tensorflow
            #from keras import backend as K
            
            session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=numIntraOpThreads, inter_op_parallelism_threads=inter_op_parallelism_threads)
            sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
            tf.compat.v1.keras.backend.set_session(sess)
            print(tf.config.list_physical_devices())
        pass
    
    def regtunemodel(self,x,y,x_test,y_test,nunit1,nunit2,xshape,activation,reg,ifbinary,loss,optimizer,metrics,cross_cri):
        model = Sequential();
        model.add(Dense(nunit1, input_dim=xshape, activation=activation))
        if reg!=0:
            model.add(Dropout(reg))
        model.add(layers.BatchNormalization())
        model.add(Dense(nunit2, activation=activation))
        model.add(layers.BatchNormalization())
        if ifbinary==0:
            model.add(Dense(1))
        if ifbinary==1:
            model.add(Dense(1, activation="sigmoid"))
        model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
        model.fit(x, y, epochs=100, validation_split = .2, verbose=0)
        predictions =model.predict(x_test)   
        if cross_cri=="mse":
            Resulttmp=numpy.array([reg,mean_squared_error(y_test,predictions)]);
        elif cross_cri=="mae":
            Resulttmp=numpy.array([reg,mean_absolute_error(y_test,predictions)]);
        elif cross_cri=="cor": 
            if numpy.std(predictions)==0:
                Resulttmp=numpy.array([reg,-3]);
            else:
                r = numpy.corrcoef(y_test[:,0],predictions[:,0]);
                Resulttmp=numpy.array([reg,-r[0,1]]);
        elif cross_cri=="auc":
            Resulttmp=numpy.array([reg,1-roc_auc_score(y_test,predictions)]);
        return Resulttmp
        
        
    def deepmodelfunc(self,xall,yall,reg=0.2,split=0.8,cross_cri="cor",loss='mse',optimizer = "adam",nunit1=50,nunit2=10,metrics=["mae","mse","accuracy"], reg_tune=1, activation="sigmoid", ifbinary=0, geneindex="model"):
        # find regularization parameters
        if reg_tune:
            x, x_test, y, y_test = train_test_split(xall, yall, test_size=0.20, shuffle=True, random_state=1234);
            xshape=x.shape[1]
            # No dropout model;
            reg=0;
            Resulttmp=self.regtunemodel(x=x, y=y, x_test=x_test, y_test=y_test, nunit1=nunit1, nunit2=nunit2, xshape=xshape, activation=activation, reg=reg, ifbinary=ifbinary, loss=loss, optimizer=optimizer, metrics=metrics,cross_cri=cross_cri)
            Result=numpy.copy(Resulttmp);
       
            # Low dropout model;
            reg=0.2
            Resulttmp=self.regtunemodel(x=x, y=y, x_test=x_test, y_test=y_test, nunit1=nunit1, nunit2=nunit2, xshape=xshape, activation=activation, reg=reg, ifbinary=ifbinary, loss=loss, optimizer=optimizer, metrics=metrics,cross_cri=cross_cri)
            Result=numpy.vstack((Result,Resulttmp));      
        
            # Mid dropout model
            reg=0.5
            Resulttmp=self.regtunemodel(x=x, y=y, x_test=x_test, y_test=y_test, nunit1=nunit1, nunit2=nunit2, xshape=xshape, activation=activation, reg=reg, ifbinary=ifbinary, loss=loss, optimizer=optimizer, metrics=metrics,cross_cri=cross_cri)
            Result=numpy.vstack((Result,Resulttmp));      
        
            # Large dropout model
            reg=0.7
            Resulttmp=self.regtunemodel(x=x, y=y, x_test=x_test, y_test=y_test, nunit1=nunit1, nunit2=nunit2, xshape=xshape, activation=activation, reg=reg, ifbinary=ifbinary, loss=loss, optimizer=optimizer, metrics=metrics,cross_cri=cross_cri)
            Result=numpy.vstack((Result,Resulttmp));      
            
            # Find optimal regs
            index=numpy.argmin(Result[:,1]) #    print(Result1[index,:])
            reg=Result[index,0];
            #print("selected dropout rate=", reg)
            
        # build model with suggested regs
        print("building models with dropout rate=",reg);
        model = Sequential(name=geneindex);
        model.add(Dense(nunit1, input_dim=xall.shape[1], activation=activation,name=geneindex+"_dense0"))
        if reg!=0:
            model.add(Dropout(reg,name=geneindex+"_drop0"))
        model.add(layers.BatchNormalization(name=geneindex+"_norm0"))
        model.add(Dense(nunit2, activation=activation,name=geneindex+"_dense1"))
        model.add(layers.BatchNormalization(name=geneindex+"_norm1"))
        if ifbinary==0:
            model.add(Dense(1,name=geneindex+"_dense2"))
        if ifbinary==1:
            model.add(Dense(1, activation="sigmoid",name=geneindex+"_dense2"))    
        model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
        model.fit(xall, yall, epochs=100, validation_split = .2, verbose=0)        
        return model, reg



    def scorefuncmean(self,model,loss,y_test,x_test,group_index,nperm):
    #model=mydnn_model;
        ypred=model.predict(x_test);
        if loss=="mse"or loss=="mean_squared_error":
            g1=(y_test[:,0]-ypred[:,0])*(y_test[:,0]-ypred[:,0]);
        if loss=="mae" or loss=="mean_absolute_error":
            g1=abs(y_test[:,0]-ypred[:,0])
        if loss=="binary_crossentropy":
            g1=-y_test[:,0]*log(ypred[:,0])-(1-y_test[:,0])*log(1-ypred[:,0])

        score=numpy.zeros((x_test.shape[0], len(group_index),nperm));
        for perm in range(0, nperm):
            num=0;
            for index in group_index:
                index=numpy.array(index);
                #print(index)
                xtmp=numpy.copy(x_test);
                indexrow=numpy.random.permutation(x_test.shape[0])
                xtmp[:,index[:]-1]=numpy.take(xtmp[:,index[:]-1], indexrow, axis=0)
                ytmp=model.predict(xtmp)
                if loss=="mse"or loss=="mean_squared_error": 
                    g2=(y_test[:,0]-ytmp[:,0])*(y_test[:,0]-ytmp[:,0]);
                if loss=="mae" or loss=="mean_absolute_error":
                    g2=abs(y_test[:,0]-ytmp[:,0]);
                if loss=="binary_crossentropy":
                    g2=-y_test[:,0]*log(ytmp[:,0])-(1-y_test[:,0])*log(1-ytmp[:,0]);  
                score[:,num,perm]=g1-g2;
                num=num+1;
            if perm==0:
                scorefinal=score[:,:,perm]
            else:
                scorefinal=scorefinal+score[:,:,perm]
        scorefinal=scorefinal/nperm;
        return scorefinal, score

    def deeptestfunc(self, xdata, ydata, group_index, nperm=20, reg=0.2, split=0.8, vfold=10, cross_cri="cor", loss='mse',optimizer = "adam", nunit1=50, nunit2=10, metrics=["mae","mse","accuracy"], reg_tune=1, activation="sigmoid", verbose=0, ifbinary=0, geneindex="model"):
    # For significance test 
        print("starting significance test")
        kf = KFold(n_splits=vfold,random_state=1234,shuffle=True); myresult=list()
        for train, test in kf.split(xdata):
            xall, x_test, yall, y_test = xdata[train], xdata[test], ydata[train], ydata[test];
            mydnn_model,mydnn_reg=self.deepmodelfunc(xall, yall, reg=reg, split=split, cross_cri=cross_cri, loss=loss, optimizer = optimizer, nunit1=nunit1, nunit2=nunit2, metrics=metrics, reg_tune=reg_tune, activation=activation, ifbinary=ifbinary);
            score,scorelist=self.scorefuncmean(model=mydnn_model,loss=loss, y_test=y_test, x_test=x_test, group_index=group_index, nperm=nperm)
            list1=[mydnn_model,mydnn_reg,score,scorelist]
            #list1=[xall, x_test, yall, y_test ]
            myresult.append(list1)
            #print(bcolors.WARNING +"fold=",len(myresult)," finished!\n"+bcolors.ENDC)
            print("running fold=",len(myresult)," finished!\n")
        
        ## rbind score and scorelist for each sample #
        scorereturn=numpy.copy(myresult[0][2])
        scorelistreturn=numpy.copy(myresult[0][3]);
        if reg_tune:
            regarray=numpy.zeros((vfold,1))
            regarray[0,0]=myresult[0][1]
        if vfold>1:
            for fold in range(1,vfold):
                scorereturn= numpy.concatenate((scorereturn, myresult[fold][2]), axis=0)
                scorelistreturn=numpy.concatenate((scorelistreturn, myresult[fold][3]), axis=0)
                if reg_tune:
                    regarray[fold,0]= myresult[fold][1];
        print("calculating p-values")
    
        ## pvalue new ###
        pvalue=numpy.ones((3,len(group_index))); ## 3*hypotheses
        # without permutations #
        indexg=20;
        for index in range(0,len(group_index)):
            if stdev(scorelistreturn[:,index,indexg-1])!=0:
                tmpmean=mean(scorelistreturn[:,index,indexg-1]);    #a1=apply(results[[2]]$scorelist[,,kk],2,mean);
                tmpstd=stdev(scorelistreturn[:,index,indexg-1]);    #b1=apply(results[[2]]$scorelist[,,kk],2,sd);
                sqrtntraintmp=numpy.lib.scimath.sqrt(xdata.shape[0]);              #pp1=pnorm((a1/b1)*sqrt(ntrain))
                statstmp=tmpmean/tmpstd*sqrtntraintmp;
                pp1=scipy.stats.norm.cdf(statstmp);
                pvalue[0,index]=pp1   
    
        for index in range(0,len(group_index)):
            if stdev(scorereturn[:,index])!=0:
                tmpmean=mean(scorereturn[:,index]);
                tmpstd=stdev(scorereturn[:,index]);
                sqrtntraintmp=numpy.lib.scimath.sqrt(xdata.shape[0]);              #pp1=pnorm((a1/b1)*sqrt(ntrain))
                statstmp=tmpmean/tmpstd*sqrtntraintmp;
                pp3=scipy.stats.norm.cdf(statstmp);
                pvalue[2,index]=pp3
    
        for index in range(0,len(group_index)):
            tmpmean=mean(scorereturn[:,index]);
            scoretmp=scorelistreturn[:,index,:];
            stscoretmp = numpy.lib.scimath.sqrt(sum(numpy.var(scoretmp, axis=1)))/xdata.shape[0];
            statstmp=tmpmean/stscoretmp;
            pp2=scipy.stats.norm.cdf(statstmp);
            pvalue[1,index]=pp2
              
        # Final models
        if reg_tune:
            reg=numpy.median(regarray);
        print("build the final model with reg="+str(reg)) 
        model = Sequential(name=geneindex);
        model.add(Dense(nunit1, input_dim=xdata.shape[1], activation=activation,name=geneindex+"l0dense"))
        if reg!=0:
            model.add(Dropout(reg,name=geneindex+"l0drop"))
        model.add(layers.BatchNormalization(name=geneindex+"l0norm"))
        model.add(Dense(nunit2, activation=activation,name=geneindex+"l1dense"))
        model.add(layers.BatchNormalization(name=geneindex+"l1norm"))
        if ifbinary==0:
            model.add(Dense(1,name=geneindex+"l2"))
        if ifbinary==1:
            model.add(Dense(1, activation="sigmoid",name=geneindex+"l2"))
        model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
        model.fit(xdata, ydata, epochs=100, validation_split = .2, verbose=verbose)
        print(tf.config.list_physical_devices('GPU'))

        return model,pvalue,scorereturn,scorelistreturn,reg
