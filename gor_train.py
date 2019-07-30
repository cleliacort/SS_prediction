#WORK: Do prediction of secondary structures using the GOR method
#USAGE: Follow the help message
#EXAMPLE:python gor.py -ids ids.txt -dssp dssp -prf profile -w 17
import os
import argparse
import pickle
import math
import sys
import numpy as np


###trasform a profile file into a matrix
def mklistpr(profile):
    try:
        pr=open(profile)
    except:
        print ("File not found. Please check the correctness of the folder contaning profiles.")
        sys.exit()
    big = []
    for line in pr:
        N=line.rstrip()
        big.append([float(x) for x in N.split('\t')])
    return big

###to add list of zero at the begin and at the end[MAKING A NEW MATRIX]
def padding(window,matrixpr):
    #window=int(window)
    for i in range (0,int(window/2)):
        matrixpr.insert(0,[0.0]*20)
        matrixpr.append([0.0]*20 )
    return matrixpr

###function that update the dic of the matrix of SS
def training(paddingpr,window,dssp,dic,R,allSS):
    ss=open(dssp).readlines()[-1].strip()
    l=[]
    for i in ss:
        if i not in l :
            if i=="-":l.append("C")
            else:l.append(i)

    window=int(window)

    for i in range (0, len(ss)):
        Sec=""
        if ss[i]=="-":
            Sec="C"
        else:
            Sec=ss[i]
        allSS[Sec]+=1
        Windo=paddingpr[i:window+i]
        for j in range(window):
            for k in range (20):
                dic[Sec][j][k]=dic[Sec][j][k]+Windo[j][k]
                R[j][k]=R[j][k]+Windo[j][k]
    return dic,R,allSS

###function to rescale the matrix of ss and the one for the residues
def rescaleMa(matrix,lista):
    for i in range (len(matrix)):
        for j in range(len(matrix[0])):
            matrix[i][j]=matrix[i][j]/lista[i]
    return matrix

###make a dic with the Information values for using directly it in the prediction
def infomation(model,residues,ss,dic):
    #residue=matrix of the frequence of residue among the window in all the Sequences
    #ss=dictionary contaning the total amount of the secondary structures found
    for i in model:
        for j in range (len(model[i])):
            for k in range(len(model[i][0])):
                dic[i][j][k]=math.log(model[i][j][k]/(residues[j][k]*ss[i]))
    return dic

###############################################################################
###############################################################################

if __name__=='__main__':
    ###Manage the arguments
    parser = argparse.ArgumentParser(description='Traning a model for Protein Secondary structure prediction using the GOR method.')
    parser.add_argument('-ids','--idsList',type=str, metavar='', required=True, help='The file containing the ids list')
    parser.add_argument('-dssp','--folderOfDssp',type=str, metavar='', required=True, help='The folder name containing the dssp files')
    parser.add_argument('-prf','--folderOfProfiles',type=str, metavar='', required=True, help='The folder containing the profile files')
    parser.add_argument('-w','--window',type=int, metavar='', required=True, help='The length of the window to extract')
    parser.add_argument('-o','--outputName',type=str, metavar='', default='gor_model.txt', required=False, help='Specify the output name')
    args=parser.parse_args()

    R=np.zeros([args.window,20])
    SSlist=['H','C','E']
    dicOverall={x: np.zeros([args.window,20]) for x in SSlist}
    allSS={x: 0 for x in SSlist}
    dicInfo={x: np.zeros([args.window,20]) for x in SSlist}

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

        #trasfomr the profile input in a matrix and then do the padding on it
        P=mklistpr(pathnmPRF)
        Pnew=padding(int(args.window),P)

        #start the traning using the padding profile and the dssp and join the result in a unique big dic
        training(Pnew,int(args.window),pathnmDSSP,dicOverall,R,allSS)

    #rescale the residue matrix
    SumRes=R.sum(axis=1)
    rescaleMa(R,SumRes)
    #rescale the overall dictionary
    for z in dicOverall:
        dicOverall[z]=rescaleMa(dicOverall[z],SumRes)
    #rescale the dictionary counting the ss in general
    tot=sum(allSS.values())
    for s in allSS:
        allSS[s]=allSS[s]/float(tot)

    #make the dic with the information values
    I=infomation(dicOverall,R,allSS,dicInfo)

    #save the I into an external file
    try:
        outfile=open(args.outputName,"wb")
    except:
        print ("Error: Pass a valid name for the output file (no the path).")
    else:
        pickle.dump(I,outfile)
    # #pickle.load(open("gor_training.txt","rb")) #to open the file in another script
