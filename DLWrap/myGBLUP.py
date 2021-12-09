import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize
import multiprocessing
from functools import partial
import pandas

class myGBLUP:
    def __init__(self, verbose=False):
        self.verbose = verbose
        pass
    
    def f_REML(self,lambda_a , argstuple):
        n_p = argstuple[0]
        theta=argstuple[1]
        theta=theta.reshape(theta.shape[0],1)
        omega_sq=argstuple[2]
        return n_p * np.log (sum (omega_sq/( theta + lambda_a ) ) ) + sum( np.log ( theta + lambda_a ) ) 

    def myBLUP(self,y,K,X=None):
        pi = 3.14159
        n = len(y)
        nall = K.shape[0]
        if X is None:
            X=np.ones((nall,1))
        p=X.shape[1]
        Z=np.zeros((nall,nall),float); np.fill_diagonal(Z,1)
    
        # Separate training from testing #
        Z=Z[range(0,n),:]
        Xall=X
        Xtest=Xall[range(n,nall),:]
        X=Xall[range(0,n),:]
    
        XtX=np.matmul(X.T,X) #t(X) %*% X, XtX <- crossprod(X, X)
        if LA.matrix_rank(XtX) < p:
            raise("X not full rank")    
        XtXinv=LA.inv(XtX) #XtXinv <- solve(XtX)
        tmp=np.zeros((n,n),float); np.fill_diagonal(tmp,1)
        S=tmp-np.matmul(np.matmul(X,XtXinv),X.T)
    
        offset = np.sqrt(n)
        Hb = np.matmul(np.matmul(Z,K),Z.T) + offset*tmp

        Hb_system_value, Hb_system_vector = LA.eigh(Hb)
        phi = Hb_system_value - offset
        phi = phi.reshape(phi.shape[0],1)
    
        if (min(phi) < -1e-06): 
            raise("K not positive semi-definite.")
    
        U = Hb_system_vector
        SHbS = np.matmul(np.matmul(S , Hb), S)
        SHbS_system_value, SHbS_system_vector = LA.eigh(SHbS)
        theta=SHbS_system_value[p:n]- offset
        theta=theta.reshape(theta.shape[0],1)
        Q=SHbS_system_vector[:,p:n]
    
        omega = np.matmul(Q.T, y) #t(x) %*% y 
        omega=omega.reshape(omega.shape[0],1)
        omega_sq = omega * omega;
        ## REML ##
        argstuple=[n-p,theta,omega_sq]
        soln = minimize(self.f_REML,[np.std(y)], args=argstuple, bounds=[(6.82001e-05,1e+09)])
        lambda_opt=soln.x
        df=n-p
    
    
        # Prediction #
        Vu_opt = sum(omega_sq/(theta + lambda_opt))/df
        Ve_opt = lambda_opt * Vu_opt
        Hinv = np.matmul(U ,(U.T/(phi+lambda_opt)))
    
        W = np.matmul(X.T,np.matmul(Hinv,X))
        beta = LA.solve(W,np.matmul(X.T,np.matmul(Hinv,y)))


        KZt = np.matmul(K,Z.T)
        KZt_Hinv = np.matmul(KZt, Hinv)
        u = np.matmul(KZt_Hinv, (y - np.matmul(X,beta)))#+np.matmul(Xall,beta)
        return u

    def GRM(self,GenoFile,param_fixed):
        print("GRM for "+GenoFile)
        ifreturn=param_fixed["ifreturn"]
        outputfile=param_fixed["outputfile"]
        allinput=param_fixed["allinput"]
        index=(allinput.tolist()).index(GenoFile)
        allgenotest=param_fixed["allgenotest"]
        GenoFile_test=allgenotest[index]
        sumpq=0
        Geno=pandas.read_csv(GenoFile)
        Geno=Geno.drop(Geno.columns[[0, 1, 2, 3, 4, 5]], axis=1);
        Geno=Geno.to_numpy();
        Geno_test=pandas.read_csv(GenoFile_test)
        Geno_test=Geno_test.drop(Geno_test.columns[[0, 1, 2, 3, 4, 5]], axis=1);
        Geno_test=Geno_test.to_numpy();
        Geno=np.concatenate((Geno, Geno_test), axis=0);
       
        freq=np.mean(Geno, axis=0)/2
        P=2*(freq-0.5);
        sumpq=sumpq+sum(freq*(1-freq))
        Geno=Geno - 1 - P 
        G=np.matmul(Geno,Geno.T)
        if ifreturn:
            return [G, sumpq]
        else:
            np.savetxt(outputfile+"_GRM_"+str(index),G,delimiter=",")
            np.savetxt(outputfile+"_GRM_"+str(index)+"sumpq",[sumpq])

    def GRMparallel(self,inputs,genotestfile,outputfile,numThreads): 
        if numThreads//2 > 10: 
            pool = multiprocessing.Pool(processes=max(1,numThreads//2))
            numThreadspara=numThreads//(numThreads//2)
            param_fixed={"outputfile":outputfile, "allinput":inputs, "allgenotest":genotestfile, "ifreturn": False}
            pool_fixed=partial(self.GRM, param_fixed=param_fixed)
            pool.map(pool_fixed, inputs)
            pool.close()
    
            for i in range(0,len(inputs)):
                grm=outputfile+"_GRM_"+str(i)
                if i==0:
                    G=pandas.read_csv(grm,sep='\t|,',engine='python',header=None)
                    G=G.to_numpy()
                    sumpq=pandas.read_csv(grm+"sumpq", sep='\t|,', engine='python', header=None)
                else:
                    G=G+(pandas.read_csv(grm,sep='\t|,',engine='python',header=None)).to_numpy()
                    sumpq=sumpq+pandas.read_csv(grm+"sumpq",sep='\t|,',engine='python',header=None)
            G=G/sumpq.to_numpy()/2
        else:
            param_fixed={"outputfile":outputfile, "allinput":inputs, "allgenotest":genotestfile, "ifreturn": True}
            for i in range(0,len(inputs)):
                grm=myGBLUP()
                if i==0:
                    [G, sumpq]=grm.GRM(inputs[i],param_fixed)
                else:
                    [G1, sumpq1]=grm.GRM(inputs[i],param_fixed)
                    G=G+G1
                    sumpq=sumpq+sumpq1
            G=G/sumpq/2
        return G
    
    