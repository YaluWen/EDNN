#!/usr/bin/env python
import pickle;
import numpy;
import os
import tensorflow;
import keras;
from keras import Sequential, Model
from keras import layers
from keras.layers import Dense, concatenate, Input, Dropout
import pandas
from sklearn.metrics import roc_auc_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from statistics import stdev
import random
import tensorflow as tf
from myGBLUP import myGBLUP


def MLPPredictionWrap(binary_outcome,GeneticDataTrain,GeneticDataTest, AssociationDir, geneindexFile, train_y_input, outputPredFile, alpha=0.05, seed_value=0,numThreads=-9, test_y_input=None):
    mlppred = MLPPrediction(seed_value=seed_value,numThreads=numThreads)
    if binary_outcome==1:
        loss="binary_crossentropy";
    else:
        loss="mse";
    [pred1,selectgenes, blup]=mlppred.prediction(GeneticDataTrain=GeneticDataTrain, GeneticDataTest=GeneticDataTest, AssociationDir=AssociationDir, geneindexFile=geneindexFile, ytrainFile=train_y_input, outputPredFile=outputPredFile, alpha=alpha, loss=loss,numThreads=numThreads,binary_outcome=binary_outcome)
    print(selectgenes)
    numpy.savetxt(outputPredFile+"summary_genes", selectgenes, fmt="%s")
    if test_y_input is not None:
        tmpdata = pandas.read_csv(test_y_input,sep='\t|,',engine='python')
        tmpdata = tmpdata.drop(tmpdata.columns[[0, 1]], axis=1);
        ytrue=tmpdata.to_numpy()
        results=mlppred.evaluate(pred=pred1, ytrue=ytrue, ytruefile=None, ifbinary=binary_outcome)
        results1=mlppred.evaluate(pred=blup, ytrue=ytrue, ytruefile=None, ifbinary=binary_outcome)
        print("accuracy="+str(results)+":"+str(results1))
        print(type(results))
        numpy.savetxt(outputPredFile+"summary_accuracy", [results, results1])
    return pred1

class MLPPrediction:
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
            from keras import backend as K        
            session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=numIntraOpThreads, inter_op_parallelism_threads=inter_op_parallelism_threads)
            sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
            tf.compat.v1.keras.backend.set_session(sess)
        pass

    def prediction(self, GeneticDataTrain, GeneticDataTest, AssociationDir, geneindexFile, ytrainFile, outputPredFile=None, alpha=0.1, optimizer = "adam", activation="sigmoid", metrics=["mae","mse","accuracy"], loss="mse", numThreads=1,binary_outcome=0):
        
        ## Load association results and find genes that are predictive given alpha ##
        genes=numpy.genfromtxt(geneindexFile,dtype=str)
        traindata=numpy.genfromtxt(GeneticDataTrain,dtype=str)
        testdata=numpy.genfromtxt(GeneticDataTest,dtype=str)
        if genes.shape[0]!=traindata.shape[0]:
            print("Genetic Training Data and the genes are not consistent!")
            raise SystemExit(0)
        if testdata.shape[0]!=traindata.shape[0]:
            print("Genetic Training Data and Testing Data are not consistent!")
            raise SystemExit(0)
        
        n_genes=genes.shape[0]
        pvalues=numpy.ones(n_genes,dtype=float);
        includegene=numpy.zeros(n_genes,dtype=int);
        
        for index in range(0, n_genes):
            files=AssociationDir+"result"+genes[index]+".npy";
            if os.path.exists(files):
                with open(files, 'rb') as handle:
                    rs=pickle.load(handle)
                try:
                    pvalues[index]=rs['pvalue'][1,0]
                except KeyError:
                    pvalues[index]=rs[1]['pvalue'][1,0]
                if pvalues[index] < alpha:
                    includegene[index]=1
        #print("All pvalues")
        print(pvalues)
        
        # load ytrain #
        tmpdata = pandas.read_csv(ytrainFile,sep='\t|,',engine='python')
        tmpdata = tmpdata.drop(tmpdata.columns[[0, 1]], axis=1);
        yall=tmpdata.to_numpy()
    
        ## Build a GBLUP model ##
        blue=myGBLUP()
        grm=blue.GRMparallel(inputs=traindata, genotestfile=testdata, outputfile=outputPredFile+"GRM", numThreads=numThreads)
        blup_pred=blue.myBLUP(y=yall,K=grm)
        blup_pred_tr=blup_pred[0:len(yall)]
        blup_pred_test=blup_pred[len(yall):len(blup_pred)]
        print("gBLUP is done!")
        
        ## Start to load models ##
        if sum(includegene)==0:
            print("No genes are selected!\nHere are the p-values for each gene.\nYou may choose to loose the criteria!")
            print("BLUP results will be used!")
            print(numpy.vstack((genes,pvalues)).T);
            predictions=blup_pred_test;
            selectgenes=0;
            #raise SystemExit(0)
        
        else:
            print("The following genes are included for the prediction model!")
            print(genes[includegene==1])           
            selectgenes=genes[includegene==1]
            traindata=traindata[includegene==1]
            testdata=testdata[includegene==1]
            models=numpy.ndarray(sum(includegene==1),dtype=object);
            modelsconcat=numpy.ndarray(sum(includegene==1),dtype=object);
            modelsinput=numpy.ndarray(sum(includegene==1),dtype=object);
            xtrain_array=numpy.ndarray(sum(includegene==1),dtype=object);
            xtest_array=numpy.ndarray(sum(includegene==1),dtype=object);
            
            for index in range(0,sum(includegene==1)):
                #print(AssociationDir+"result_model"+selectgenes[index])
                model1=keras.models.load_model(AssociationDir+"result_model"+selectgenes[index])
                model1.trainable=False
                models[index]=model1
                modelsconcat[index]=model1.get_layer(model1.layers[-2].name).output
                modelsinput[index]=model1.input;
        
                tmpdata = pandas.read_csv(traindata[index],sep='\t|,',engine='python')
                tmpdata = tmpdata.drop(tmpdata.columns[[0, 1, 2, 3, 4, 5]], axis=1);
                xalltmp=tmpdata.to_numpy()
                xtrain_array[index]=xalltmp;
        
                tmpdata= pandas.read_csv(testdata[index],sep='\t|,',engine='python')
                TestID=tmpdata.iloc[:,[0, 1, 2, 3, 4, 5]]
                tmpdata =tmpdata.drop(tmpdata.columns[[0, 1, 2, 3, 4, 5]], axis=1);
                xtesttmp=tmpdata.to_numpy()
                xtest_array[index]=xtesttmp

            blupinput = Input((1), name='blupinput')
            model_blup=Dense(1,input_dim=1,activation=activation)(blupinput)
            #model_blup=layers.BatchNormalization()(model_blup)
            tmp=modelsconcat.tolist()
            #tmp.append(model_blup)
            
            #model_concat = concatenate(modelsconcat.tolist(), axis=-1)
            model_concat = concatenate(tmp, axis=-1)
            model_concat = Dense(50, activation=activation,name="combine1")(model_concat)
            model_concat = Dense(10, activation=activation,name="final_layer1")(model_concat)
            #model_concat = Dense(5, activation=activation,name="final_layer3")(model_concat)
            model_concat = concatenate([model_concat, model_blup], axis=-1)
            model_concat = Dense(100, activation=activation,name="final_layer2")(model_concat)
            model_concat = Dropout(0.2,name="final_drop")(model_concat)
            model_concat = Dense(30, activation=activation,name="final_layer")(model_concat)
            if binary_outcome==1:
                model_concat = Dense(1, activation=activation,name="final_layera")(model_concat)
            else:
                model_concat = Dense(1, activation="linear",name="final_layera")(model_concat)

            tmp=modelsinput.tolist()
            tmp.append(blupinput)
            model = Model(inputs=tmp, outputs=model_concat,name="final")
            #model = Model(inputs=modelsinput.tolist(), outputs=model_concat,name="final")
            model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
            
            # Fit model #
            xtrain_all=xtrain_array.tolist()
            xtrain_all.append(blup_pred_tr)
            model.fit(xtrain_all, yall, epochs=100, validation_split = .2, verbose=0) 
            xtest_all=xtest_array.tolist()
            xtest_all.append(blup_pred_test)
            predictions =  model.predict(xtest_all)


        #if sum(includegene)==1:
        #    selectgenes=genes[includegene==1]
        #    traindata=traindata[includegene==1]
        #    testdata=testdata[includegene==1]
            
        #    modelsconcat=numpy.ndarray(sum(includegene==1),dtype=object);
        #    modelsinput=numpy.ndarray(sum(includegene==1),dtype=object);
        #    xtrain_array=numpy.ndarray(sum(includegene==1),dtype=object);
        #    xtest_array=numpy.ndarray(sum(includegene==1),dtype=object);
            
        #    tmpdata= pandas.read_csv(testdata[0],sep='\t|,',engine='python')
        #    TestID=tmpdata.iloc[:,[0, 1, 2, 3, 4, 5]]
        #    tmpdata =tmpdata.drop(tmpdata.columns[[0, 1, 2, 3, 4, 5]], axis=1);
        #    x_test=tmpdata.to_numpy()
            
        #    model=keras.models.load_model(AssociationDir+"result_model"+selectgenes[0])
        #    predictions =  model.predict(x_test.tolist()) 
            tmpdata= pandas.read_csv(testdata[index],sep='\t|,',engine='python')###myself
            TestID=tmpdata.iloc[:,[0, 1, 2, 3, 4, 5]]###myself
            print(TestID.shape)
            print(predictions.shape)    
            FinalPred=numpy.concatenate((TestID,predictions),axis=1)
            if (outputPredFile is not None):
                numpy.savetxt(outputPredFile,FinalPred)

        return [predictions, selectgenes,blup_pred_test]
            #tensorflow.keras.utils.plot_model(model,to_file='demo.png',show_shapes=True)
            #print(model.summary())
            #print(model.trainable_weights)
            #print("weights:", len(model.weights))
            #print("trainable_weights:", len(model.trainable_weights))
            #print("non_trainable_weights:", len(model.non_trainable_weights))

    def evaluate(self, pred, ytrue=None, ytruefile=None, ifbinary=0):
        if ytrue is None and ytruefile is None:
            print("The true value is not provided! Please provide the true value!")
        if ytrue is not None and ytruefile is not None:
            print("The true value provided will be used!")
        if ytrue is None and ytruefile is not None:
            tmpdata = pandas.read_csv(ytruefile,sep='\t|,',engine='python')
            tmpdata = tmpdata.drop(tmpdata.columns[[0, 1]], axis=1);
            ytrue=tmpdata.to_numpy()
        if ifbinary==1:
            result=[roc_auc_score(ytrue,pred),0]
        
        if ifbinary==0:
            if numpy.std(pred)==0:
                r=-3
            else:
                r = numpy.corrcoef(ytrue[:,0],pred[:,0]);
                r = r[0,1]
            mse = mean_squared_error(ytrue,pred)
            result=[r, mse]
        return result
