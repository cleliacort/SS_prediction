import sys
import re



#file_input=sys.argv[1]
#list_of_ids=sys.argv[2]
#print "hello"
###MAKE THE DICTIONARY OF ALL THE POSSIBLE TUPLE COMBINATION
def mk_dic_from_pred(all_pred):
    lista=['H','C','E']
    dic={}
    for i in lista:
        for j in range(len(lista)):
            dic[(i,lista[j])]=0
###UPDATE THE DICTIONARY WITH THE INFORMATION COMING FROM THE PREDICTION
    unique_fl=open(all_pred)
    for line in unique_fl:
        dic[line.rstrip().split()[1],line.rstrip().split()[2]]+=1
    return dic

##ACCURACY OVERALL
def Q3_overall(dic):
    lista=['H','C','E']
    N=sum(dic.values())
    ptot=0
    for ss in lista:
        ptot+=dic[ss,ss]
    Q=(ptot/float(N))*100
    return round(Q,2)

    ###SENSITIVITY,PPV,MCC
def indices(dic):
    lista=['H','C','E']
    N=sum(dic.values())
    ssIndx={}
    for elm in lista:
        #print elm
        c=0
        u=0
        o=0
        for j in range(len(lista)):
            if elm==lista[j]:
                c=dic[elm,elm]
            else:
                u+=dic[elm,lista[j]]
                o+=dic[lista[j],elm]
        sen=c/float(c+u)
        ppv=c/float(c+o)

        n=N-c-u-o
        mcc=((c*n)-(o*u))/((c+o)*(c+u)*(n+o)*(n+u))**0.5
        ssIndx[elm]=[sen,ppv,mcc]
    return ssIndx

    ###SOV
def SOV(all_pred):
    lista=['H','C','E']
    all_seq_SS={x: 0 for x in lista}
    result={x: 0 for x in lista}

    ##make track of the ids
    unique_fl=open(all_pred)
    ids=[]
    for l in unique_fl:
        code=l.split()[0].rstrip()
        if code not in ids:
            ids.append(code)
    #print ids

    for i in ids:
        unique_fl=open(all_pred)
        dssp=""
        pred=""
        #print i
        for j in unique_fl:
            #print j
            if i.rstrip()==j.split()[0].rstrip():
                dssp=dssp+j.split()[1].rstrip()
                pred=pred+j.split()[2].rstrip()
        # print dssp
        # print pred
    #     # for ss in lista:
        for ss in lista:
            if ss in dssp and ss in pred:
                all_seq_SS[ss]+=1

                pattern=ss+"+"
                sdssp=[]
                edssp=[]
                for match in re.finditer(pattern,dssp):
                    sdssp.append(match.start())
                    edssp.append(match.end())
                # print i.rstrip()ss+str(sdssp)
                #print i.rstrip()+ss+str(edssp)
                #print "\n"

                spred=[]
                epred=[]
                for match in re.finditer(pattern,pred):
                    spred.append(match.start())
                    epred.append(match.end())
                # print i.strip()+ss+str(spred)
                # print i.rstrip()+ss+str(epred)

                #if len(sdssp)!=0 and len(spred)!=0:
                    #all_seq_SS[ss]+=1
                minoverl=[]
                maxoverl=[]
                lenO=[]
                deltaForAll=[]
                N=0
                #numOver=0
                #flag=0
                for z in range(len(sdssp)):
                    flag=0
                    x=range(sdssp[z],edssp[z])
                    #print i.rstrip()+ss+str(x)
                    for k in range(len(spred)):
                        y=range(spred[k],epred[k])
                        minimum=len(list(set(x)&set(y)))
                        #print ss+str(minov)
                        #print ss+str(y)
                        if minimum!=0:
                            flag+=1
                            #N+=len(x)
                            #print N
                            #numOver+=1
                            #make a list with the min information
                            minoverl.append(minimum)
                            #make a list with the max information
                            maxlist=[]
                            maxlist.extend((sdssp[z],edssp[z],spred[k],epred[k]))
                            maxrange=range(min(maxlist),max(maxlist))
                            maxoverl.append(len(maxrange))
                #print ss+str(maxoverl)
                            #make a list with the length of the dssp fragment(observed)
                            lenO.append(len(x))
                #print ss+str(lenO)
                            #make a list with the delta values
                            delta=[]
                            delta.extend((len(maxrange)-minimum,minimum,len(x)/float(2),len(y)/float(2)))
                            #print ss+str(delta)
                            deltaForAll.append(min(delta))
                    if flag==0:
                        N+=len(x)
                    else:
                        N+=len(x)*flag
                #print deltaForAll
                #print ss+str(N)

                SOV=0
                for p in range(len(minoverl)):
                    SOV+=((minoverl[p]+deltaForAll[p])/float(maxoverl[p]))*lenO[p]
                #print ss+str(SOV)
                SOV=100*(1/float(N))*SOV
                #print ss+str(SOV)
                SOV=round(SOV,2)
                #print ss+str(SOV)

                if SOV!=0:
                    result[ss]=result[ss]+SOV
                #     print i.rstrip()+"\t"+ss+"\t"+str(SOV)
                # #print i.rstrip()
    for ssFin in result:
        result[ssFin]=result[ssFin]/float(all_seq_SS[ssFin])
	
    return result


if __name__=='__main__':
    file_input=sys.argv[1]
    dic=mk_dic_from_pred(file_input)
    Q=Q3_overall(dic)
    Indx=indices(dic)
    res,seq=SOV(file_input)
    #print res
