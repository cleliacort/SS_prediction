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
parser.add_argument('-o','--outputName',type=str, metavar='', default='svm_model', required=False, help='Specify the output name')
args=parser.parse_args()

w=args.window
X_train=[]
Y_train=[]
codf={'H':1,'E':2,'C':3,'-':3}

ids=open(args.idsList)
for i in ids:
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

    P=mklistpr(pathnmPRF)
    Pnew=padding(w,P)

    #make a unique big sublist with all the possible window for each position
    ss=open(pathnmDSSP).readlines()[-1].strip()

    for i in range (0, len(ss)):
        xt=[]
        Windo=Pnew[i:w+i]
        for pos in range (len(Windo)):
            for aa in range (len(Windo[pos])):
                xt.append(Windo[pos][aa])
    #upload this sublist to the final one
        X_train.append(xt)
    #update the unique big list for the codification of the dssp
    for s in ss:
        Y_train.append(codf[s])
#print (len(Y_train))
#print (len(X_train))
#######################################################################################

#SVM part
#mySVC = svm.SVC(C=8.0, kernel='rbf', gamma=0.3)
#mySVC = svm.SVC(C=2.0, kernel='rbf', gamma=0.3)
#mySVC = svm.SVC(C=2.0, kernel='rbf', gamma=2.0)
mySVC = svm.SVC(C=8.0, kernel='rbf', gamma=2.0)
mySVC.fit(X_train,Y_train)

try:
   out= gzip.open(args.outputName+".gz",'w')
except:
   print ("Error: Pass a valid name for the output file (no the path).")

pickle.dump(mySVC, out)

