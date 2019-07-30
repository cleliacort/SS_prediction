import argparse
from gor_train import *
import os
from sklearn import svm
import pickle, gzip

parser = argparse.ArgumentParser(description='Traning a model for Protein Secondary structure prediction using the SVM.')
parser.add_argument('-ids','--idsList',type=str, metavar='', required=True, help='The file containing the ids list')
parser.add_argument('-dssp','--folderOfDssp',type=str, metavar='', required=True, help='The folder name containing the dssp files')
parser.add_argument('-prf','--folderOfProfiles',type=str, metavar='', required=True, help='The folder containing the profile files')
parser.add_argument('-w','--window',type=int, metavar='', required=True, help='The length of the window to extract')
parser.add_argument('-md','--modelName',type=str, metavar='', required=True, help='The name of the svm model')
parser.add_argument('-o','--outputName',type=str, metavar='', default='svm_pred', required=False, help='Specify the name for the output file')
args=parser.parse_args()

w=args.window
codf={'H':1,'E':2,'C':3,'-':3}

try:
    predFile=open(args.outputName+"_all.txt","w+")
except:
    print ("Error: Pass a valid name for the output file (no the path).")


mySVC = pickle.load(gzip.open(args.modelName, 'r'))
ids=open(args.idsList)
for i in ids:
    X_test=[]
    namep=""
    chain=""
    pathnmDSSP=""
    pathnmPRF=""

    #make the right name on the bases of the ids input (PDB ID format or jpred based on the domain)
    if ":" in i:
        #manipulate the name to trasfor in the rigth form
        namep=i.split(':')[0].rstrip().lower()
        chain=i.split(':')[1].rstrip()
        #create the path where take the input file
        pathnmDSSP=os.path.join(args.folderOfDssp,namep+"."+chain+".dssp")
        pathnmPRF=os.path.join(args.folderOfProfiles,namep+"."+chain+".prf")
    else:
         namep=i.rstrip()
         #create the path where take the input file
         pathnmDSSP=os.path.join(args.folderOfDssp,namep+".dssp")
         pathnmPRF=os.path.join(args.folderOfProfiles,namep+".prf")
   
    #trasform the profile in a list and padding it
    P=mklistpr(pathnmPRF)
    Pnew=padding(w,P)

    #make a unique big sublist with all the possible window for each position
    ss=open(pathnmDSSP).readlines()[-1].strip()
    for s in range (0, len(ss)):
        xt=[]
        Windo=Pnew[s:w+s]
        for pos in range (len(Windo)):
            for aa in range (len(Windo[pos])):
                xt.append(Windo[pos][aa])
    #upload this sublist to the final one
        X_test.append(xt)

##################################################################################################

    #MAKE THE PREDICTION 
    y_pred=mySVC.predict(X_test)
    #print the result in the correct form 
    listss=list(ss) #to avoid problem in the indeces between list and string 
    for p in range (len(y_pred)):
        if listss[p]=='-':listss[p]='C'
        predFile.write(i.rstrip()+"\t"+listss[p]+"\t"+list(codf.keys())[list(codf.values()).index(y_pred[p])]+"\n")
