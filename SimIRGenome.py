#!/usr/bin/python
from functools import reduce
import os

import numpy as np
from numpy.random import choice

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def biased_DNA_generator(GCRatio,size):
    DNAstr="".join(np.random.choice(["A","T","C","G"], size, p=[0.5-GCRatio/2,0.5-GCRatio/2,GCRatio/2,GCRatio/2]))
    return Seq(DNAstr, generic_dna)


def SimulateGenome(args):
    GCRatio = args.GCRatio
    IDRNum = args.IDRNum
    GSize = args.GSize 
    pre = args.pre
    np.random.seed(args.seed)
    name = os.path.basename(pre)

    # Randomly generate IR length from 13 to 45 base pair
    IRSize=np.random.randint(low=13,high=45, size=IDRNum)

    # Randomly generate the inverted DNA Region length from 31 to 
    IDRSize=np.random.randint(low=31,high=750, size=IDRNum)

    # Randomly generate the base pairs between inverted DNA regions
    NonIRSize=GSize-sum(IRSize)*2 - sum(IDRSize)
    a = np.sort(np.random.choice(range(NonIRSize), IDRNum,replace=False))
    LastSize=NonIRSize-a[-1]
    LastSeq= biased_DNA_generator(GCRatio,LastSize)
    IDRInsertSize=np.insert(np.diff(a), 0, a[0])

    # Randomly generate DNA sequences for all three types of regions
    IRList = [biased_DNA_generator(GCRatio,x) for x in IRSize]
    IRRevList = [x.reverse_complement() for x in IRList]
    IDRList = [biased_DNA_generator(GCRatio,x) for x in IDRSize]
    IDRRevList = [x.reverse_complement() for x in IDRList]
    IDRInsertList = [biased_DNA_generator(GCRatio,x) for x in IDRInsertSize]
    # Concatenate the regions and generate the forward one
    SeqList = [ IDRInsertList[i]+IRList[i]+IDRList[i]+IRRevList[i] for i in range(IDRNum)]
    product = reduce((lambda x, y: x + y), SeqList) +LastSeq
    productRecord = SeqRecord(product, id=name, description="")
    SeqIO.write(productRecord, pre+".fwd.fasta", "fasta")

    # Concatenate the regions and generate the reverse one
    SeqRevList = [ IDRInsertList[i]+IRList[i]+IDRRevList[i]+IRRevList[i] for i in range(IDRNum)]
    productRev = reduce((lambda x, y: x + y), SeqRevList) + LastSeq
    productRevRecord = SeqRecord(productRev, id=name+"-rev", description="")
    SeqIO.write(productRevRecord, pre+".rev.fasta", "fasta")

    # Report the position info
    f = open(pre+".pos.tab", "w")
    cumpos=0
    for i in range(IDRNum):
        cumpos+=IDRInsertSize[i]
        print >>f, name+"\t" +'\t'.join(str(x) for x in [cumpos,cumpos+IRSize[i], cumpos+IRSize[i]+IDRSize[i], cumpos+IRSize[i]+IDRSize[i]+IRSize[i]])
        cumpos+=IRSize[i]+IDRSize[i]+IRSize[i]

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Simulate DNA sequence with inverted repeats')
    parser.add_argument('-l',"--length", type=int,  help='Length of the genome',required=True, dest="GSize",metavar='')
    parser.add_argument('-n',"--number", type=int,  help='Number of inverted DNA regions',required=True, dest="IDRNum",metavar='')
    parser.add_argument('-r',"--gcratio", type=float,  help='GC Ratio',required=True, dest="GCRatio",metavar='')
    parser.add_argument('-p', '--prefix',type=str, help='Prefix of output',required=True, dest='pre', metavar='')
    parser.add_argument('-s',"--seed", type=int,  help='Seed for random number generator',required=False, dest="seed",metavar='')
    args = parser.parse_args()

    SimulateGenome(args)
