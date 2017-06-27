#Python code for ChIA-PET reads anchored on promoter regions

#usage: python ChiaAnchorPromoter.py input.txt ouput.txt

#refGene.sql and refGene.txt.gz files are downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
#refGene table in MySQL database

#input.txt (downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/)
#Example: wgEncodeGisChiaPetMcf7CtcfInteractionsRep1.txt
#591     chr1    805223  887355  chr1:805223..805726-chr1:886758..887355,2       200     .       805223  887355  16711680        2       503,597 0,81535
#591     chr1    854905  873967  chr1:854905..855680-chr1:873268..873967,2       200     .       854905  873967  16711680        2       775,699 0,18363
#                                                                               .
#                                                                               .

####################################################################
## Step 0 | Install and define the functions
####################################################################

import MySQLdb,sys

#Function for parsing genomic postion from strings in input.txt 
def ParsePosition(posString):
    posStringSplit = posString.strip().split(':')
    posChromosome = posStringSplit[0]
    posStart = long(posStringSplit[1].split('..')[0])
    posEnd = long(posStringSplit[1].split('..')[1])
    return [posStart,posEnd,posChromosome]

#Function for overlap beteewn two genomic positions
def ChechOverlapPlusInfo(aList,bList):
    if aList[0] <= bList[0] and aList[1] >= bList[0]:
        overlap = 1
    elif aList[0] <= bList[1] and aList[1] >= bList[1]:
        overlap = 1
    elif aList[0] >= bList[0] and aList[1] <= bList[1]:
        overlap = 1
    elif aList[0] <= bList[0] and aList[1] >= bList[1]:
        overlap = 1
    else:
        overlap = 0
    return [overlap,aList,bList]

####################################################################
## Step 1 | Get every promoter positions
####################################################################

#database connection
conn=MySQLdb.connect(host="localhost", user="USER", passwd="PASSWD", db="hg19")
cur=conn.cursor()

#query and make chromosome list
table = "refGene"
sql1 = "SELECT DISTINCT chrom FROM " + table + ""
cur.execute(sql1)
res1 = cur.fetchall()
chr_list = [i[0] for i in res1]

dicPos_promoter = {}    #initiate dictionary (keys: chromosomes, values: promoter position list)

for chrom in chr_list:

    #query promoter positions from refGene table
    sql = "select strand,txStart,txEnd,name,name2 from " + table + " where chrom=" +"\'" + chrom + "\';"
    cur.execute(sql)
    recs = cur.fetchall()

    posPromoterList = []
    
    #get every promoter positions and append to dicPos_promoter
    for rec in recs:
        strand = rec[0]
        txStart = rec[1]
        txEnd = rec[2]
        nmNum = rec[3]
        geneSymbol = rec[4]

        if(strand == "+"):
            posPromoterList.append([txStart-2000,txStart+500,nmNum,geneSymbol])
        
        if(strand == "-"):
            posPromoterList.append([txEnd-500,txEnd+2000,nmNum,geneSymbol])

    dicPos_promoter[chrom] = posPromoterList

cur.close()
conn.close()

print "done...loading Promoter postion"

####################################################################
## Step 2 | Get ChIA-PET reads overlaped with promoters
####################################################################

inFile = sys.argv[1]    #inputFile
outFile = sys.argv[2]   #ouputFile

#Checking overlap
f_chia = open(inFile)
fw = open(outFile,'w')

for line in f_chia:

    #parse input.txt
    lineSplit = line.strip().split('\t')
    posPair = lineSplit[3]
    posPairSplit = posPair.strip().split('-')
    pos1 = ParsePosition(posPairSplit[0])
    pos2 = ParsePosition(posPairSplit[1].split(',')[0])

    #Chech overlap with promoters
    if((dicPos_promoter.has_key(pos1[2]))):
        for posPromoter in dicPos_promoter[pos1[2]]:
            overlapPlusInfo = ChechOverlapPlusInfo(pos1,posPromoter)
            if(overlapPlusInfo[0] == 1):
                nmNumber = overlapPlusInfo[2][2]
                symbol = overlapPlusInfo[2][3]
                wStr = posPairSplit[0] + "\t" + (posPairSplit[1].split(',')[0]) + "\t" + nmNumber + "\t" + symbol + "\n"
                fw.write(wStr)  #print
            
    if((dicPos_promoter.has_key(pos2[2]))):
        for posPromoter in dicPos_promoter[pos2[2]]:
            overlapPlusInfo = ChechOverlapPlusInfo(pos2,posPromoter)
            if(overlapPlusInfo[0] == 1):
                nmNumber = overlapPlusInfo[2][2]
                symbol = overlapPlusInfo[2][3]
                wStr = (posPairSplit[1].split(',')[0]) + "\t" + posPairSplit[0] + "\t" + nmNumber + "\t" + symbol + "\n"
                fw.write(wStr)  #print
    
fw.close()
f_chia.close()
