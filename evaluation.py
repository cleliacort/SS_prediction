import sys
import os
from indices import *
import argparse
import numpy as np
import matplotlib.pyplot as plt

#directory=sys.argv[1]
#list_of_ids=sys.argv[2]

parser = argparse.ArgumentParser(description='Compute the Evaluation of a model.')
parser.add_argument('-fd','--predictionFolder',type=str, metavar='', required=False, default='.', help='The folder name in which there are the predictions output')
parser.add_argument('-rn','--numberOfRun',type=int, metavar='', required=True, help='The number of set for cross validation')
args=parser.parse_args()

if args.predictionFolder:
    os.chdir(args.predictionFolder)

path="."
files=os.listdir(path)
lista=['H','C','E']

crossN=args.numberOfRun
run=0
resQ=[]
resIndxSingl={x: np.zeros([crossN,3]) for x in lista}
resIndxAll={x: np.zeros([3])  for x in lista}

resSOVAll={x:0 for x in lista}
resSOVSingl={x: np.zeros([crossN])  for x in lista}
#print resSOVSingl
numSeq={x:0 for x in lista}
#print resIndx
#print resIndx

for name in files:
    if "all.txt" in name:
        run+=1
        dic=mk_dic_from_pred(name)
        resQ.append(Q3_overall(dic))

        #compute the indices (MCC,PPV,SEN) and update the resul indices dictionary
        Indx=indices(dic)
        #print Indx
        for key in Indx:
            for ind in range(len(Indx['H'])):
                resIndxSingl[key][run-1][ind]=resIndxSingl[key][run-1][ind]+Indx[key][ind]

        #compute the SOV and update the resul indices dictionary of SOV and of num of seq to be considered
        sov=SOV(name)
        for ss in sov:
            resSOVSingl[ss][run-1]=sov[ss]

##RESULT
#compute the final accuracy Q
avgQ=np.mean(resQ)
SDq=np.std(resQ)
if SDq!=0:
	print ("The Three-class accuracy is\t"+str(avgQ)+"+/-"+str(round(SDq,2))+"%\n")
else:
	print ("The Three-class accuracy is\t"+str(avgQ)+"%\n")

#normalize the indices (MCC,PPV,SEN)by dividing for the different run of cross validation
SDIndx={x: np.zeros([3])  for x in lista}
resIndxTot={x: np.zeros([3])  for x in lista}
#print resIndxSingl
for ss in resIndxSingl:
     for ind in range(len(resIndxSingl['H'][0])):
         SDIndx[ss][ind]=np.std(resIndxSingl[ss][:,ind])
         resIndxTot[ss][ind]=np.mean(resIndxSingl[ss][:,ind])
#print resIndxTot

I=['SEN','PPV','MCC']
for res in resIndxTot:
	print ("For\t"+res+":")
	for ind in range(len(resIndxTot[res])):
		if SDIndx[res][ind]!=0:
			print ("the "+str(I[ind])+" is\t"+str(round(resIndxTot[res][ind],2))+"+/-"+str(round(SDIndx[res][ind],3)))
		else:
			print ("the "+str(I[ind])+" is\t"+str(round(resIndxTot[res][ind],2)))

print ("\n")
#normalize the SOV by dividing for the number of sequences in which the secondary structure appears
SDSov={x:0 for x in lista}
SovTot={x:0 for x in lista}
for run in resSOVSingl:
    SovTot[run]=round(np.mean(resSOVSingl[run]),2)
    SDSov[run]=round(np.std(resSOVSingl[run]),2)

for ss in SovTot:
	if SDSov[ss]!=0:
		print ("The SOV for "+ss+" is\t"+str(SovTot[ss])+"+/-"+str(SDSov[ss]))
	else:
		print ("The SOV for "+ss+" is\t"+str(SovTot[ss]))
