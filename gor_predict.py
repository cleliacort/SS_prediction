#WORK:compute the prediction of the SS using the GOR model previous trained
#USAGE: Follow the help message
#EXAMPLE:python gor_predict.py -ids ids.txt -m gor_model.txt -nf profile -w 17
import pickle
import argparse
from gor_train import *
import os

parser = argparse.ArgumentParser(description='Compute the Prediction of Protein Secondary Structure using a GOR model.')
parser.add_argument('-ids','--idsList',type=str, metavar='', required=True, help='The file containing the ids list')
parser.add_argument('-dssp','--folderOfDssp',type=str, metavar='', required=True, help='The folder name containing the dssp files')
parser.add_argument('-nf','--folderOfProfiles',type=str, metavar='', required=True, help='Specify the name of the folder containing the profiles')
parser.add_argument('-m','--GORmodel',type=str, metavar='', required=True, help='The file containing the GOR model')
parser.add_argument('-w','--window',type=int, metavar='', required=True, help='The length of the window using to train the model')
parser.add_argument('-o','--outputName',type=str, metavar='', default='gor_prediction', required=False, help='Specify the output name')
args=parser.parse_args()

InfoMa=pickle.load(open(args.GORmodel,"rb"))
try:
	all_pred=open(args.outputName+"_all.txt","w+")
except:
	print ("Error: Pass a valid name for the output file (no the path).")

#all_pred=open("all_predicted_ss.txt","w+")

ids=open(args.idsList)
for i in ids:
    ###ADJUST THE NAME AND THE PATHS
    namep=""
    chain=""
    pathnmPRF=""
    pathnmDSSP=""
    c=0

    if ":" in i:
        #manipulate the name to trasfor in the rigth form
        namep=i.split(':')[0].rstrip().lower()
        chain=i.split(':')[1].rstrip()
        #create the path where take the input file
        pathnmPRF=os.path.join(args.folderOfProfiles,namep+"."+chain+".prf")
        pathnmDSSP=os.path.join(args.folderOfDssp,namep+"."+chain+".dssp")
        c=1
    else:
        #manipulate the name to trasfor in the rigth form
        namep=i.rstrip()
        #create the path where take the input file
        pathnmPRF=os.path.join(args.folderOfProfiles,namep+".prf")
        pathnmDSSP=os.path.join(args.folderOfDssp,namep+".dssp")

    ###TRASFORM IN THE RIGHT FORM THE NEW PROFILE ON WHICH DO THE PREDICTION
    P=mklistpr(pathnmPRF)
    len_seq=len(P)#i save the len of P at this point because with the padding i will change P
    padding(len(list(InfoMa.values())[0]),P)

    ###MAKE THE PREDICTION AND WRITE THE OUT FILE IN A SPECIFI FORMAT [IDS] [ORIGNAL_SS] [PREDICTED_SS]
    dssp=open(pathnmDSSP).readlines()[-1].strip()
    for pos in range (len_seq):
        dic={}
        Windo=P[pos:args.window+pos]
        for singleSS in InfoMa:#cicla sui dizionari
            s=0
            for row in range (len(InfoMa[singleSS])):
                for col in range (len(InfoMa[singleSS][row])):
                    s=s+InfoMa[singleSS][row][col]*Windo[row][col]
            dic[singleSS]=s

        pred_ss=max(dic, key=dic.get)

        orig_dssp=dssp[pos]
        if dssp[pos]=="-":
            orig_dssp="C"

        all_pred.write(i.rstrip()+"\t"+orig_dssp+"\t"+pred_ss+"\n")
