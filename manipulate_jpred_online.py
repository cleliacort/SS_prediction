#USAGE:
#EXAMPLE:python ./fasta/jpred_site.py -ids Job_ids.txt -dssp dssp -pred fasta
import tarfile
import os
import argparse


parser = argparse.ArgumentParser(description='Manipulate file coming from the result of Jpred online')
parser.add_argument('-ids','--idsList',type=str, metavar='', required=True, help='The file containing the ids list')
parser.add_argument('-dssp','--folderOfDssp',type=str, metavar='', required=True, help='The folder name containing the dssp files')
parser.add_argument('-pred','--folderOfPrediction',type=str, metavar='', required=True, help='Specify the name of the folder containing the profiles')
args=parser.parse_args()
all_pred=open("jpred_online_all.txt","w+")

#open the list of Job ids
ids=open(args.idsList)
for i in ids:
    string=""
    pathnmPred=os.path.join(args.folderOfPrediction,i.rstrip())

    #open the taz.zip folder and search  for the file contaning the prediction and the one with the dssp information coming from dssp
    folder=tarfile.open(pathnmPred+"/"+i.rstrip()+".tar.gz")
    jnet=folder.extractfile(i.rstrip()+".jnet")
    fasta=folder.extractfile(i.rstrip()+".input")

    #extract info from the dssp file
    code=fasta.readlines()[0].rstrip()[1:]
    id_dssp=code.split(":")[0].lower()
    chain_dssp=code.split(":")[1]
    dssp_file_name=id_dssp+"."+chain_dssp
    pathnmDSSP=os.path.join(args.folderOfDssp,dssp_file_name+".dssp")

    dssp=open(pathnmDSSP).readlines()[-1].strip()

    #extract info from the prediction file and make the final string to update on the output file
    for line in jnet.readlines():
        if line.startswith("jnetpred"):

            s=line.index(":")+1
            pred=(line[s:].rstrip()).split(",")
            for ss in range (len(pred)-1):
                if pred[ss]=="-":
                    pred[ss]="C"

                all_pred.write(dssp_file_name+"\t"+dssp[ss]+"\t"+pred[ss]+"\n")
