#WORK:It take as input the DSSP file and return us one file with the fasta sequence and the second file with secondary structure
#USAGE: PYTHON [NAME_SCRIPT_PY] [FOLDER OF DSSP FILES] [IDS LIST]
#EXAMPLE: python manipulate_dssp.py DSSP_blind final_blind_dset.txt
import sys
import os

DSSP=sys.argv[1] #name directory of DSSP files
IDS=sys.argv[2] #IDS of the final blind set

lg={'E':'E','B':'E','H':'H','G':'H','I':'H','T':'C','S':'C',' ':'C'}

os.mkdir("dssp")
os.mkdir("fasta")

#CREATE A LIST WITH THE IDS INSIDE THE FILE PASSED AS INPUT
#ids=[]
ids_list=open(IDS)
#for line in ids_list:
#    ids.append(line.rstrip().lower())


#WORK ON A SINGLE FILE DSSP TO MAKE DUE DIFFERENT FILES WITH INFORMATION COMING
#FROM THE DSSP FILE
for i in ids_list.xreadlines():
    i = i.rstrip()

    namep=i.split(':')[0].rstrip().lower()
    chain=i.split(':')[1].rstrip()
    #APRE FILE CORRISPONDENTE AD OGNI IDS
    pathname=os.path.join(DSSP,namep+".dssp")
    fl=open(pathname)

    #CREA NEW FILE PER OGNI IDS
    fsname=namep+"."+chain+".fasta"
    ssname=namep+"."+chain+".dssp"

    F=open(os.path.join("fasta",fsname),"w+")
    S=open(os.path.join("dssp",ssname),"w+")

    #INSERT THE HEAD FOR THE TWO FILE FOR EACH IDS
    F.write(">"+namep.upper()+":"+chain+"\n")
    S.write(">"+namep.upper()+":"+chain+"\n")

    #RIEMPI I FILES
    stringFS=""
    stringSS=""
    flag=-1
    for line in fl:

        if flag==+1 and line[11]==chain:
            FS=line[13]
            stringFS+= str(FS)

            SS=lg[line[16]]
            stringSS+=str(SS)

        elif "#  RESIDUE" in line:
            flag=+1

    F.write(stringFS)
    S.write(stringSS)
