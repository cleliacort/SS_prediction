# Secondary Structure prediction

The knowledge gap between the sequence and structure is an open problem in the study of proteins. In the post-genomics era, this gap is exponentially enlarging due to the decrease of the cost to sequencing.

Although more sequence data are available, the same are completely useless without their functional annotation. To overcame this problem, several bioinformatics methods were developed for trying to extract reliable functional information from the huge "black hole" of sequence data.

In this project, I analyzed two methods to predict the secondary structures, GOR and SVM. The first one is a probabilistic procedure; while the second one is a machine learning technique. Depending also on their background mathematical theory, they present advantages and disadvantages. 

The script in this repository allow to execute the training and predicting procedures for both GOR and SVM. Moreover, there is also a script to evaluate the performance of the prediction which in practice compute the main scoring indexes (accuracy,MCC,PPV,SEN). 
In the same directory, I upload also two script to respectively manipulate the DSSP and PSSM file.
