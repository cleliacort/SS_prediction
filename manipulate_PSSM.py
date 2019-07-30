#WORK:extract from the pssm file the information on the profile
#USAGE:PYTHON [NAME_SCRIPT] [FOLDER OF PSSM] [IDS LIST]
#EXAMPLE:python manipulate_PSSM.py PSSM ids.txt
import sys
import os

PSSM=sys.argv[1]#take the name of the DSSP folder
IDS=sys.argv[2]#take a list of IDS

os.mkdir("profile")

ids_list=open(IDS)

#open the ids list
for i in ids_list.xreadlines():
    namep=""
    chain=""
    profi=""
    pathname=""
    n=0
    #if you use the PDB ID
    if ":" in i:
        #save the ID and the chain
        namep=i.split(':')[0].rstrip().lower()
        chain=i.split(':')[1].rstrip()

        #open the file pssm using the ID and the chain coming from the list
        pathname=os.path.join(PSSM,namep+"."+chain+".pssm")
        profi=namep+"."+chain+".prf"
        n=1
    #if you use the domain IDs
    else:
        #clean the ID from the possible \n at the end
        namep = i.rstrip()

        #open the file pssm
        pathname=os.path.join(PSSM,namep+".pssm")
        profi=namep+".prf"



    try:
        fl=open(pathname)
    except:
        if n==1:
            print "There is no pssm file for the id"+"\t"+namep+chain
        else:
            print "There is no pssm file for the id"+"\t"+namep
    else:



        #open this new file on which we want to write
        P=open(os.path.join("profile",profi),"w+")
        #print P

        #scroll the file pssm and take the information we want
        flag=-1
        start=0
        for line in fl.readlines()[:-6]:
            #to skip the first lines
            if " pseudocounts" in line:
                flag=1

            #to skip the header
            elif flag==1:
                #take the index of A to start from this point since it change for different pssm file
                flag=2
                start= line.rfind("A")-1 #take the index of A to start from this point since it change
                #for different pssm file=> find the first occurence from the rigth

            #to add the lines
            elif flag==2:
                l=line[start:start+81]#select the part of interest
                l = [float(j)/100.0 for j in l.split()]#compute the percentage

            #print in a new file the new values with the percentage
                P.write('\t'.join(str(i) for i in l))
                P.write('\n')

        fl.close()
        P.close()
