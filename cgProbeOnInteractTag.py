#Python code to CpGs (450k probe) on interaction tags of ChIA-PET reads 

#Usage: python cgProbeOnInteractTag.py ChiaAnchorPromoter_output.txt CpG_position_file.txt ouput.txt

#ChiaAnchorPromoter_output.txt 
#chr10:8095767..8096378	chr1:8586064..8586898	NR_104329	GATA3-AS1
#chr10:8095767..8096378	chr1:8586064..8586898	NM_001002295	GATA3
#                                        .
#                                        .

#CpG_position_file.txt
#probe name	gene symbol	chromosome	position
#cg03123289	A1BG	19	58858635
#cg22286978	A1BG	19	58858806
#cg10734734	A1BG	19	58859021
#                 .
#                 .


####################################################################
## Step 0 | Install and define the functions
####################################################################

import sys

#Function for parsing genomic postion from strings in ChiaAnchorPromoter_output.txt
def ParsePosition(posString):
    posStringSplit = posString.strip().split(':')
    posChromosome = posStringSplit[0]
    posStart = long(posStringSplit[1].split('..')[0])
    posEnd = long(posStringSplit[1].split('..')[1])
    return [posStart,posEnd,posChromosome]

####################################################################
## Step 1 | Get every postions of interaction tags
####################################################################
interactome_file = sys.argv[1]  #ChiaAnchorPromoter_output.txt 
cgProbe_file = sys.argv[2]      #CpG_position_file.txt
outFile = sys.argv[3]           #output file

f_inter = open(interactome_file)

dicPosInter = {}    #initiate dictionary (keys: chromosomes, values: interaction tag position list)

for line in f_inter:
    #parse ChiaAnchorPromoter_output.txt
    lineSplit = line.strip().split('\t')
    posAnchor = ParsePosition(lineSplit[0])
    posInter = ParsePosition(lineSplit[1])
    nmNum = lineSplit[2]
    symbol = lineSplit[3]

    #get every promoter positions and append to dicPos_promoter
    listPosInter = [posInter[0],posInter[1],posInter[2],posAnchor,nmNum,symbol]

    if(dicPosInter.has_key(posInter[2])):
        dicPosInter[posInter[2]].append(listPosInter)
    else:
        dicPosInter[posInter[2]] = []
        dicPosInter[posInter[2]].append(listPosInter)

f_inter.close()

####################################################################
## Step 2 | Get CpGs overlaped with interaction tags
####################################################################
f_cg = open(cgProbe_file)
fw = open(outFile,'w')

header = f_cg.readline()

for line in f_cg:

    #parse CpG_position_file.txt
    lineSplit = line.strip().split('\t')
    probeID = lineSplit[0]
    nearGene = lineSplit[1]
    chrom = "chr" + lineSplit[2]
    posCG = long(lineSplit[3])

    #Chech overlap with interaction tags
    if(dicPosInter.has_key(chrom)):
        for posInteract in dicPosInter[chrom]:
            if((posCG >= posInteract[0]) and (posCG <= posInteract[1])):
                wStr = probeID + "\t" + nearGene + "\t" + chrom + "\t" + str(posCG)\
                       + "\t" + posInteract[2] + ":" + str(posInteract[0]) + "-" + str(posInteract[1])\
                       + "\t" + posInteract[4] + "\t" + posInteract[5]\
                       + "\t" + posInteract[3][2] + ":" + str(posInteract[3][0]) + "-" + str(posInteract[3][1])\
                       + "\n"
                fw.write(wStr)  #print

fw.close()
f_cg.close()
