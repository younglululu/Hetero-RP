import time
import math
import sys
import os

import numpy as np
import pandas as pd
import argparse

from scipy import stats
from scipy.sparse import coo_matrix,find
from scipy.optimize import nnls
from sklearn.neighbors import kneighbors_graph
from itertools import repeat

import diptest #from diptest import diptest
from cvxopt import solvers, matrix, spdiag


def arguments():
    parser = argparse.ArgumentParser(description='Heterogeneity Rescaling Pursuit (Hetero-RP): Towards Enhanced and Interpretable Clustering/Classification in Integrative Genomics')
    group1 = parser.add_argument_group('Input_Options', 'Options related to parsing the input data.')
    group1.add_argument('-skip_row_header', action='store_true', dest='skip_row_header',default=False,  help=("Whether the first column indicates the header. Default behavior is not set if there is no header explicitly given in the first column. Explicitly set -skip_row_header to be able to skip the header for downstream analysis"))
    group1.add_argument('-skip_column_header', action='store_true', dest='skip_column_header', default=False, help=("Whether the first row indicates the header. Default behavior is not set if there is no header explicitly given in the first row. Explicitly set -skip_column_header to be able to skip the header for downstream analysis"))
    group1.add_argument('-object_in_row', action='store_true', dest='object_in_row', default=False, help=("Whether each row represents an object and each column represents a feature. By default each row represents a feature and each column represents an object. Explicitly set -object_in_row to let each row represents an object and each column represents a feature."))
    group1.add_argument('-sep', dest='sep', default=',', help=("Delimiter to use. Default is ','."))
    group1.add_argument('-data_files', nargs='+', dest='data_files', help=("The data files, each containing a table where each row correspond to a feature, and each column correspond to an object. The row number can vary but the column number must be consistent. All values are separated by sep."))
    group1.add_argument('-signed_graph_file', dest='signed_graph_file', help=("The signed graph encoding either the positive-links and negative-links information, one row for one edge in the format: object_A<sep>object_B<sep>weight. The edge is undirected. The index of objects starts from 1. The weight is 1 for positive-links and -1 for negative links."))
    
    group2 = parser.add_argument_group('Preprocessing_Options', 'Options related to preprocessing the features.')
    group2.add_argument('-check_feature_validity', action='store_true', dest='check_feature_validity',default=False,  help=("Check whether the feature is informative or not. Default behavior is not set if all features are believed informative. Explicitly set -check_feature_validity to be able to check features one by one."))
    group2.add_argument('-label_file', dest='label_file', help=("The label file used to check the feature validity (optional). Each row is an integer number indicating the group of corresponding object, the total number of objects must be consistent with the data files."))
    group2.add_argument('-p_value', type=float, default=0.05, dest='p_value', help=("The p-value used in multi-modality test or t-test depending on the availability of the label file."))
    #group2.add_argument('-label_as_signed_graph', action='store_true', dest='label_as_signed_graph',default=False,  help=("Explicitly set -label_as_signed_graph to build auxiliary signed graph based on the label information, which can be override by -signed_graph_file. "))
    
    group3 = parser.add_argument_group('Program_Options', 'Options related to Hetero-RP.')
    group3.add_argument('-use_knn', action='store_true', dest='use_knn',default=False,  help=("When combine the auxiliary knowledge signed graph with the k-nearest-neighbor graph constructed from the input data. Default behavior is not set if there is sufficient auxiliary knowledge. Explicitly set -use_knn to be able to combine the auxiliary knowledge signed graph with the k-nearest-neighbor graph."))
    group3.add_argument('-k', type=int, default=0, dest='k', help=("k in k-nearest-neighbor graph. Default is sqrt(n) where n is the number of objects."))
    group3.add_argument('-gammaVal', type=float, default=0, dest='gammaVal', help=("The parameter controls the trade-off between the auxiliary knowledge signed graph and the k-nearest-neighbor graph constructed from the input data. Default value is provided to achieve the balanced contribution between two."))
    group3.add_argument('-lambdaVal', type=float, default=0, dest='lambdaVal', help=("The parameter shrinks weight towards unit and towards each other. If not set, the parameter will be automatically chosen."))
    
    group4 = parser.add_argument_group('Output_Options', 'Options related to the output.')
    group4.add_argument('-output_dir', dest='output_dir', help=("The output folder, storing the output result. If not specified, the output folder is set as the current directory."))
    group4.add_argument('-viz', action='store_true', dest='viz',default=False,  help=("Whether the visualized plot is generated or not. Default behavior is not set. Explicitly set -viz to be able to generate the visualized plot."))
    
    args  = parser.parse_args()
    print("skip_row_header:\t"+str(args.skip_row_header));
    print("skip_column_header:\t"+str(args.skip_column_header));
    print("object_in_row:\t"+str(args.object_in_row));
    print("sep:\t"+str(args.sep));
    print("data_files:\t"+str(args.data_files));
    print("signed_graph_file:\t"+str(args.signed_graph_file));
    print("check_feature_validity:\t"+str(args.check_feature_validity));
    print("label_file:\t"+str(args.label_file));
    print("p_value:\t"+str(args.p_value));
    print("use_knn:\t"+str(args.use_knn));
    print("k:\t"+str(args.k));
    print("gammaVal:\t"+str(args.gammaVal));
    print("lambdaVal:\t"+str(args.lambdaVal));
    print("output_dir:\t"+str(args.output_dir));
    print("viz:\t"+str(args.viz));
    
    if len(args.data_files) < 1:
        parser.error("You must give at least one data file!")
        sys.exit(0)
    
    return args

def maximum (A, B):
    BisBigger = A-B
    BisBigger.data = np.where(BisBigger.data < 0, 1, 0)
    return A - A.multiply(BisBigger) + B.multiply(BisBigger)
    
def numpy_to_cvxopt_matrix(A):
    if isinstance(A, np.ndarray):
        if A.ndim == 1:
            return matrix(A, (A.shape[0], 1), 'd')
        else:
            return matrix(A, A.shape, 'd')
    else:
        return A

def buildKnnGraph(X, k):
    A = kneighbors_graph(X, n_neighbors=k, mode='distance', n_jobs=-1);
    nzIdx = np.nonzero(A); print(nzIdx)
    varX = np.sum(np.square(np.std(X,axis=0))); print(varX)
    A[nzIdx[0],nzIdx[1]] = np.exp(-np.power(A[nzIdx[0],nzIdx[1]],2)/(2.2472*varX* np.power(n,-0.4))); #del nzIdx;
    print(A.max())
    return maximum( A, A.T ); 

def getY(X, graphMat):
    [n,p]=X.shape; minVal = graphMat.min(); 
    decay = 0.9; maxIter = 10; graphMat1 = graphMat.copy(); Y = np.zeros(p,);
    
    for iter in range(maxIter):
        L = coo_matrix((np.array(np.sum(graphMat1, axis=0)).flatten(), (np.arange(n), np.arange(n))), shape=(n,n)).tocsr() - graphMat1; 
        Y = np.zeros(p,);
        for featIdx in range(p): Y[featIdx]=np.dot(X[:,featIdx].T*L,X[:,featIdx]);
        selectedIdx = np.where(Y>0)[0]; selectedLen = len(selectedIdx); del L,selectedIdx;
        
        if selectedLen < 0.8*p:
            graphMat1[graphMat1<0] *= decay; 
        else:
            break;
    Y[Y<0]=0; del graphMat1;
    return Y;

def heteroRP(X, graphMat, lambdaVal):
    [n,p]=X.shape;
    Y = getY(X, graphMat); selectedIdx = np.where(Y>0)[0]; selectedLen = len(selectedIdx);
    
    P = spdiag(numpy_to_cvxopt_matrix(Y[selectedIdx]+lambdaVal))
    q = numpy_to_cvxopt_matrix(Y[selectedIdx])
    G = spdiag(numpy_to_cvxopt_matrix(-1*np.ones((selectedLen,1))))
    h = numpy_to_cvxopt_matrix(np.ones((selectedLen,1)))
    A = numpy_to_cvxopt_matrix(np.ones((1,selectedLen)))
    b = matrix(np.array([0]), (1,1), 'd')
    
    sol = solvers.qp(P,q,G,h,A,b); newCoeff=np.array(sol['x']).flatten();
    weightVec = np.zeros(p,); 
    weightVec[selectedIdx]=newCoeff; 
    weightVec = weightVec+1;
    print(weightVec);
    return weightVec;

def heteroRPIter(X, graphMat):
    [n,p]=X.shape; 
    Y = getY(X, graphMat); selectedIdx = np.where(Y>0)[0]; selectedLen = len(selectedIdx);
    r = min(n,selectedLen); currSigma = 10; B = stats.t.ppf(1-np.sqrt(selectedLen)/(2*r*np.log(r)),selectedLen-1); lambda_0 = B / np.sqrt(selectedLen-1+B*B);
    tol = 1e-5; iter1 = 0; maxIter = 20; diffVal = float('Inf'); currCoeff = np.zeros((selectedLen,1));
    
    P = spdiag(numpy_to_cvxopt_matrix(Y[selectedIdx]+2*lambda_0*currSigma))
    q = numpy_to_cvxopt_matrix(Y[selectedIdx])
    G = spdiag(numpy_to_cvxopt_matrix(-1*np.ones((selectedLen,1))))
    h = numpy_to_cvxopt_matrix(np.ones((selectedLen,1)))
    A = numpy_to_cvxopt_matrix(np.ones((1,selectedLen)))
    b = matrix(np.array([0]), (1,1), 'd')
    
    while diffVal > tol and iter1 < maxIter:
        iter1 = iter1 + 1;
        sol = solvers.qp(P,q,G,h,A,b); newCoeff=np.array(sol['x']).flatten(); #newCoeff=newCoeff.reshape(selectedLen,1)
        diffVal = np.sum(np.abs(newCoeff-currCoeff)); print("iteration="+str(iter1)+"\t diff="+str(diffVal));
        currSigma = Y.T.dot(((newCoeff+1)**2)); currSigma = np.sqrt(currSigma*selectedLen); 
        currCoeff = newCoeff; P = spdiag(numpy_to_cvxopt_matrix(Y+2*lambda_0*currSigma)); 
    weightVec = np.zeros(p,); 
    weightVec[selectedIdx]=newCoeff;
    weightVec = weightVec+1;
    print(weightVec);
    return weightVec;

if __name__ == '__main__':  
    args = arguments()
    
    dataMat = []; dataNameList = []; dataDimList = []; featNameList = [];
    for dataIdx in range(len(args.data_files)):
        print("loading data file... "+args.data_files[dataIdx]);
        dataHeader = pd.read_csv(args.data_files[dataIdx],sep=args.sep,nrows=1); 
        
        colStartIdx = 0; headerVal = None;
        if args.skip_row_header: colStartIdx=1; 
        if args.skip_column_header: headerVal='infer'; 
        currDataMat = pd.read_csv(args.data_files[dataIdx],sep=args.sep,header=headerVal, usecols=range(colStartIdx, dataHeader.shape[1])).as_matrix(); del dataHeader;
        if args.object_in_row: currDataMat = currDataMat.T;
        
        if dataIdx == 0: 
            dataMat = currDataMat;
        else: 
            dataMat = np.vstack((dataMat,currDataMat)); 
            
        dataNameList.append(os.path.splitext(os.path.basename(args.data_files[dataIdx]))[0]);
        dataDimList.append(currDataMat.shape[0]);
    
    [p,n] = dataMat.shape; print("Success loading data. p="+str(p)+" n="+str(n)); 
    featureValidFlags = np.ones(p,); 
    weightVec = np.zeros(p,);

    if args.check_feature_validity:
        print("Checking feature validity..."); 
         
        if not args.label_file:
            print("Turn to unsupervised mode because the label file doesn't exist...")
             
            for featIdx in range(p):
                [dipVal, pVal] = diptest.diptest(dataMat[featIdx,:]); 
                if pVal > (1.0-args.p_value): featureValidFlags[featIdx] = 0;
                else: featureValidFlags[featIdx] = 1;
        else:
            print("Turn to supervised mode because the label file exists...")
            labelVec = pd.read_csv(args.label_file, header=None, usecols=[0]).as_matrix().flatten();
            uniqLabelVec = np.unique(labelVec); 
             
            for featIdx in range(p):
                parameterList = []; dataVec=dataMat[featIdx,:]
                for labelIdx in range(len(uniqLabelVec)):
                    parameterList.append(dataVec[np.where(labelVec==uniqLabelVec[labelIdx])]);
                output=getattr(stats, "f_oneway")(*parameterList); print(str(output.pvalue))
                if output.pvalue > (1.0-args.p_value): featureValidFlags[featIdx] = 0;
                else: featureValidFlags[featIdx] = 1;
                del dataVec,parameterList,output;
            del labelVec,uniqLabelVec;
         
        print(featureValidFlags)
    else:
        print("Skip checking feature validity.")
     
    validFeatIdx = np.array(np.where(featureValidFlags>0)).flatten();
     
    if args.signed_graph_file:
        print("loading auxiliary signed graph file... "+args.signed_graph_file);
         
        edgelist = pd.read_csv(args.signed_graph_file,sep=args.sep,header=None).as_matrix();
        
        A_sup = coo_matrix((edgelist[:,2],(edgelist[:,0]-1, edgelist[:,1]-1)), shape=(n,n)).tocsr(); 
        A_sup = A_sup+A_sup.T; del edgelist;
        A = A_sup;
    else:
        print("Skip auxiliary signed graph file.");
     
    if args.use_knn:
        print("constructing k-nearest-neighbor graph from the input data...");
        k = int(round(np.sqrt(n))); 
        if args.k > 0: k = args.k; 
        print("k is set as: "+str(k));
        A_unsup = buildKnnGraph(dataMat[validFeatIdx,:].transpose(), k);
         
        if args.signed_graph_file:
            print("val1: "+str(np.sum(A_sup)));
            print("val2: "+str(np.sum(A_unsup)));

            gammaVal = np.sum(np.sum(A_sup*A_unsup,axis=0)) / np.sum(np.sum(A_unsup*A_unsup,axis=0)); 
            #if np.isnan(gammaVal): gammaVal = np.sum(np.sum(A_sup,axis=0)) / np.sum(np.sum(A_unsup,axis=0)); 
            print("gamma is set as: "+str(gammaVal));
            
            A = A_sup + gammaVal*A_unsup; 
        else:
            A = A_unsup;
    else:
        print("Skip constructing k-nearest-neighbor graph from the input data.");
     
    if 'A_sup' in locals(): del A_sup;
    if 'A_unsup' in locals(): del A_unsup;
    if not ('A' in locals()): 
        print("Error! A should exists!!");
        sys.exit(0)
     
    if args.lambdaVal == 0:
        weightVec1 = heteroRPIter(dataMat[validFeatIdx,:].transpose(), A);
    else:
        weightVec1 = heteroRP(dataMat[validFeatIdx,:].transpose(), A, args.lambdaVal);
        
    weightVec[validFeatIdx] = weightVec1;
    weightVec = weightVec.reshape(p,1);

    
    startIdx=0; 
    for dataIdx in range(len(args.data_files)):
        print("Writing weights to files...");
        outputDir = os.path.dirname(os.path.abspath(args.data_files[dataIdx]));
        if args.output_dir: outputDir = args.output_dir;
            
        featURL = os.path.join(outputDir, "feat_"+dataNameList[dataIdx]+".csv");
        np.savetxt(featURL, weightVec[startIdx:startIdx+dataDimList[dataIdx]], delimiter=",",fmt='%.6f',header="weight");
        startIdx = startIdx + dataDimList[dataIdx];
    
    if args.viz:
        print("Generating visualized plot of feature weights...");
        outputTmpFile = "visual.html_part2"; outputVizFile = "visual.html";
        if args.output_dir: 
            outputTmpFile = os.path.join(outputDir, outputTmpFile);
            outputVizFile = os.path.join(outputDir, outputVizFile);
        
        with open(outputTmpFile, "w") as text_file:
            text_file.write("var nameArr = [");
            for dataIdx in range(len(args.data_files)):
                if dataIdx > 0:  text_file.write(',');
                text_file.write(','.join(list(repeat('"'+dataNameList[dataIdx]+'"',dataDimList[dataIdx]))));
            text_file.write("];\n");
            
            text_file.write("var startArr = [");
            for dataIdx in range(len(args.data_files)):
                if dataIdx > 0:  text_file.write(',');
                text_file.write(','.join(str(x) for x in range(dataDimList[dataIdx])));
            text_file.write("];\n");
            
            text_file.write("var valArr = [");
            text_file.write(','.join("{:.6f}".format(float(x)) for x in weightVec));
            text_file.write("];\n");
            
            text_file.write("var tagArr = [");
            for dataIdx in range(len(args.data_files)):
                if dataIdx > 0:  text_file.write(',');
                text_file.write(','.join('"'+dataNameList[dataIdx]+'_'+str(x+1)+'"' for x in range(dataDimList[dataIdx])));
            text_file.write("];\n");
            
            text_file.write("var lenArr = [");
            text_file.write(','.join(str(x) for x in dataDimList))
            text_file.write("];\n");
            
            text_file.write("var labelArr = [");
            text_file.write(','.join('"'+str(x)+'"' for x in dataNameList))
            text_file.write("];\n");
            
            text_file.close();
        os.system("cat visual.html_part1 "+outputTmpFile+" visual.html_part3 > "+outputVizFile);
        os.system("rm "+outputTmpFile)
    else:
        print("Skip generating visualized plot of feature weights.");
        
    print("Done.");
    
