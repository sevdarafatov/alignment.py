## For command line arguments, I imported sys module.

import sys

##For pairwise local alignment, I have installed the package Biopython, and used local alignment function,
## and formal_alignment function to be able to show alignment.
##https://github.com/biopython/biopython

from Bio import pairwise2

from Bio.pairwise2 import format_alignment



## The sequences are taken from txt files, which are given as inputs in terminal etc. after setting the directory.

G1f=open("%s" %sys.argv[1],"r")
G2f=open("%s" %sys.argv[2],"r")
RSf=open("%s" %sys.argv[3],"r")

G1=G1f.readline()
G2=G2f.readline()
RS=RSf.readline()


G1f.close()
G2f.close()
RSf.close()


## Generating VRS (containing all 64 triplets) and counting triplets in RS:


VRS = {
    'ATA':0, 'ATC':0, 'ATT':0, 'ATG':0,
    'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0,
    'AAC':0, 'AAT':0, 'AAA':0, 'AAG':0,
    'AGC':0, 'AGT':0, 'AGA':0, 'AGG':0,
    'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0,
    'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0,
    'CAC':0, 'CAT':0, 'CAA':0, 'CAG':0,
    'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0,
    'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0,
    'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0,
    'GAC':0, 'GAT':0, 'GAA':0, 'GAG':0,
    'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0,
    'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0,
    'TTC':0, 'TTT':0, 'TTA':0, 'TTG':0,
    'TAC':0, 'TAT':0, 'TAA':0, 'TAG':0,
    'TGC':0, 'TGT':0, 'TGA':0, 'TGG':0,
    }

for i in range (len(RS)-2):

    VRS[RS[i:i+3]]+=1



## Calculating Manhattan distance between VRS and VG1:


VG1 = {
    'ATA':0, 'ATC':0, 'ATT':0, 'ATG':0,
    'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0,
    'AAC':0, 'AAT':0, 'AAA':0, 'AAG':0,
    'AGC':0, 'AGT':0, 'AGA':0, 'AGG':0,
    'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0,
    'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0,
    'CAC':0, 'CAT':0, 'CAA':0, 'CAG':0,
    'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0,
    'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0,
    'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0,
    'GAC':0, 'GAT':0, 'GAA':0, 'GAG':0,
    'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0,
    'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0,
    'TTC':0, 'TTT':0, 'TTA':0, 'TTG':0,
    'TAC':0, 'TAT':0, 'TAA':0, 'TAG':0,
    'TGC':0, 'TGT':0, 'TGA':0, 'TGG':0,
    }


dG1=[]

M1={}

## First window of RS-size Manhattan distance:

for i in range (len(RS)-2):
    
    VG1[G1[i:i+3]]+=1


for key in VRS.keys():

    M1[key]= abs(VRS[key] - VG1.get(key))


M1t= sum(M1.values())

dG1.append([0, M1t])

## Updating VG1 and adding Manhattan distances to dG1 list:

for k in range (1,len(G1)-len(RS)+1):
    
    VG1[G1[k:k+3]]+=-1
    VG1[G1[k+len(RS)-3:k+len(RS)]]+=1

    for key in VRS.keys():

        M1[key]= abs(VRS[key] - VG1.get(key))


    M1t= sum(M1.values())

    dG1.append([k, M1t])
    

##print (dG1)


a=min(x[1] for x in dG1)

G1_min=[]

for i in range(len(dG1)):

    if dG1[i][1]== a:
        G1_min.append([i,a])

print (G1_min)

## Same procedure applied for G2:

VG2 = {
    'ATA':0, 'ATC':0, 'ATT':0, 'ATG':0,
    'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0,
    'AAC':0, 'AAT':0, 'AAA':0, 'AAG':0,
    'AGC':0, 'AGT':0, 'AGA':0, 'AGG':0,
    'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0,
    'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0,
    'CAC':0, 'CAT':0, 'CAA':0, 'CAG':0,
    'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0,
    'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0,
    'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0,
    'GAC':0, 'GAT':0, 'GAA':0, 'GAG':0,
    'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0,
    'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0,
    'TTC':0, 'TTT':0, 'TTA':0, 'TTG':0,
    'TAC':0, 'TAT':0, 'TAA':0, 'TAG':0,
    'TGC':0, 'TGT':0, 'TGA':0, 'TGG':0,
    }


dG2=[]

M2={}

## First window of RS-size Manhattan distance:

for i in range (len(RS)-2):
    
    VG2[G2[i:i+3]]+=1


for key in VRS.keys():

    M2[key]= abs(VRS[key] - VG2.get(key))


M2t= sum(M2.values())

dG2.append([0, M2t])


## Updating VG2 and adding Manhattan distances to dG2 list:

for k in range (1,len(G2)-len(RS)+1):
    
    VG2[G2[k:k+3]]+=-1
    VG2[G2[k+len(RS)-3:k+len(RS)]]+=1

    for key in VRS.keys():

        M2[key]= abs(VRS[key] - VG2.get(key))


    M2t= sum(M2.values())

    dG2.append([k, M2t])
    


    
##Finding the reference sequence-sized windows with minimum distance and storing them in their
##respective lists:


G2_min=[]

b=min(y[1] for y in dG2)

for n in range(len(dG2)):

    if dG2[n][1]== b:
        G2_min.append([n,b])

print (G2_min)

## Recalling the best windows from G1 and G2:
##(I took the first best windows even if they had multiple windows with min distance)
##(I could also have recalled random one)

G1_loc=G1[G1_min[0][0]:G1_min[0][0]+len(RS)]

G2_loc=G2[G2_min[0][0]:G2_min[0][0]+len(RS)]

##The pairwise local alignment between G1 and G2 windows:


alignments= pairwise2.align.localms(G1_loc, G2_loc, 4, -3, -2, -2)

##There may be multiple best alignments with the same score, so I used for loop
##to show all of them:

for a in pairwise2.align.localms(G1_loc, G2_loc, 4, -3, -2, -2):

    print (format_alignment(*a))
