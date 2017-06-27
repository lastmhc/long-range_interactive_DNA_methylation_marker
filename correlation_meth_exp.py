#Python code to caculate Pearson correaltion of selected CpG-long range gene pairs

#usage: python correlation_meth_exp.py methFile.txt expFile.txt meth_exp_listFile.txt output.txt

#methFile.txt
#Probe_ID	TCGA-A7-A0D9-01A	TCGA-A7-A0D9-11A	TCGA-A7-A13E-01A	TCGA-A7-A13E-11A	TCGA-A7-A13F-01A    ...
#cg25831075	0.6015	0.2208	0.9311	0.3093	0.9121  ...
#cg12875720	0.3286	0.8968	0.9196	0.9145	0.5265  ...
#cg18557131	0.7803	0.0164	0.8267	0.0381	0.5666  ...
#                           .
#                           .

#expFile.txt
#GeneSymbol	TCGA-A7-A0D9-01A	TCGA-A7-A0D9-11A	TCGA-A7-A13E-01A	TCGA-A7-A13E-11A	TCGA-A7-A13F-01A    ...
#OPA3	4.6091	4.1437	4.49	4.7574	4.1193	4.5463	...
#CHST3	3.8905	6.9281	5.602	6.4923	4.254	7.0067	...
#HOXA10	2.1971	5.4742	0.5471	4.788	0.8207	4.3945	...
#                           .
#                           .


#meth_exp_listFile.txt
#cgProbe    nearGene   longRangeGene 
#cg11376305	SPHK2	FAM83E
#cg16754364	SPN	QPRT
#cg04245057	TBX2	PFKFB3
#               .
#               .

####################################################################
## Step 0 | Install the functions
## Need to install scipy and numpy
####################################################################
import sys
from scipy.stats import *
from numpy import *

####################################################################
## Step 1 | Get methylation values of CpGs
####################################################################
methFile = sys.argv[1]              #methFile.txt
expFile = sys.argv[2]               #expFile.txt
meth_exp_listFile = sys.argv[3]     #meth_exp_listFile.txt
outFile = sys.argv[4]

f_meth = open(methFile)

header_meth = f_meth.readline()
headerSplit_meth = header_meth.strip().split('\t')
headerSplit_meth.pop(0)

dicMeth = {}    #initiate dictionary (keys: CpG probe ID, values: meth value list)

for line_meth in f_meth:
    #parse methFile.txt
    lineSplit_meth = line_meth.strip().split('\t')
    cgID = lineSplit_meth[0]
    lineSplit_meth.pop(0)
    listMeth = map(float,lineSplit_meth)
    
    #get meth value list of CpG
    dicMeth[cgID] = listMeth

f_meth.close()

####################################################################
## Step 2 | Get expression values of genes
####################################################################
f_exp = open(expFile)

header_exp = f_exp.readline()
headerSplit_exp = header_exp.strip().split('\t')
headerSplit_exp.pop(0)

dicExp = {}     #initiate dictionary (keys: gene symbol, values: exp value list)

for line_exp in f_exp:
    #parse expFile.txt
    lineSplit_exp = line_exp.strip().split('\t')
    cgID = lineSplit_exp[0]
    lineSplit_exp.pop(0)
    listExp = map(float,lineSplit_exp)
    
    #get exp value list of gene symbol
    dicExp[cgID] = listExp

f_exp.close()

####################################################################
## Step 3 | Caclualte the correaltion of CpG-long range gene pairs
####################################################################
f_list = open(meth_exp_listFile)
fw = open(outFile, 'w')
header_list = f_list.readline()

for line_list in f_list:
    #parse meth_exp_listFile.txt
    lineSplit_list = line_list.strip().split('\t')
    cgProbe = lineSplit_list[0]
    nearGene = lineSplit_list[1]
    longRangeGene = lineSplit_list[2]

    #calculate the correaltion and p-value
    if(dicMeth.has_key(cgProbe) and dicExp.has_key(longRangeGene)):
        pearson = pearsonr(dicMeth[cgProbe],dicExp[longRangeGene])
        cor = str(pearson[0])
        p_val = str(pearson[1])
        wStr = line_list.strip() + "\t" + cor + "\t" + p_val + "\n"
        fw.write(wStr)      #print

fw.close()
f_list.close()
