#!/usr/bin/python
# coding=utf-8

'''
Written by Isabela Almeida
Created on Aug 28, 2023
Last modified on Aug 29, 2023    
 
Description: Takes a FASTQ file, counts the bases in each
sequence and outputs sequence to file n+ or n depending
on sequence size.    

<Envinroment: VSC>                                   
<OS: Linux (Ubuntu)>                                          
'''

from pandas import DataFrame
import pandas as pd
from Bio import SeqIO
import argparse
import sys
import gzip
import io

# Get command-line args
userInput = argparse.ArgumentParser(description="Takes a FASTQ file, counts the bases in each"
                                 " sequence and outputs sequence to file n+ or n depending"
                                 " on sequence size.")
requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--fastq', action='store',
                          required=True, type=argparse.FileType('r'),
                          help='Target fastq file')
requiredNamed.add_argument('-n', '--numberbases', action='store',
                      required=True, type=int, default='20000',
                      help='number of bases to separate sequences into different files.'
                      'E.g. 20k bp should be 20000')
requiredNamed.add_argument('-b', '--biggerthann', action='store',
                           required=True,
                           type=str, help='Path to file to save sequences with a size'
                           ' bigger than n ')
requiredNamed.add_argument('-s', '--smallerthann', action='store',
                           required=True,
                           type=str, help='Path to file to save sequences with a size'
                           ' smaller or equal to n ')

## Import user-specified command line values.
args = userInput.parse_args()
fastq_file = args.fastq
n_value = args.numberbases
out_n_plus = args.biggerthann
out_up_to_n = args.smallerthann

## Create output file
file_n_plus = open(out_n_plus, "w+" )
file_up_to_n= open(out_up_to_n, "w+" )

## Get Target fastq sequences
with gzip.open(fastq_file, "rb") as handle:
    text_wrapper = io.TextIOWrapper(handle, encoding="utf-8")

target_sequences = SeqIO.parse(text_wrapper,"fastq")

## Loop through each fastq sequence
n_plus = up_to_n = num = 0
seq0_1k = seq1k = seq1k_2k = seq2k_5k = seq5k_10k = seq10k_20k = seq20k_30k = seq30k_40k = seq40k_50k = seq50kplus = 0

for fastq in target_sequences:
    num += 1
    seqName, sequence = fastq.description, str(fastq.seq)
    thisfastq = []
    seqLen = len(sequence)
    if seqLen <1000: seq0_1k += 1
    elif seqLen == 1000: seq1k +=1
    elif (seqLen > 1000) and (seqLen <= 2000): seq1k_2k += 1
    elif (seqLen > 2000) and (seqLen <= 5000): seq2k_5k += 1
    elif (seqLen > 5000) and (seqLen <= 10000): seq5k_10k += 1
    elif (seqLen > 10000) and (seqLen <= 20000): seq10k_20k += 1
    elif (seqLen > 20000) and (seqLen <= 30000): seq20k_30k += 1
    elif (seqLen > 30000) and (seqLen <= 40000): seq30k_40k += 1
    elif (seqLen > 40000) and (seqLen <= 50000): seq40k_50k += 1
    elif (seqLen > 50000): seq50kplus += 1
    
    if seqLen > n_value:
        n_plus += 1
        file_n_plus.write("> {0}\n".format(seqName))
        file_n_plus.write("{0}\n".format(sequence))
    else:
        up_to_n += 1
        file_up_to_n.write("> {0}\n".format(seqName))
        file_up_to_n.write("{0}\n".format(sequence))

# Print info about the results to terminal
print('Sequences\tValue')
print('Total\t', num)
print('n+\t', n_plus)
print('<=n\t', up_to_n)
print('<1k\t',seq0_1k)
print('==1k\t',seq1k)
print('>1k, <=2k\t',seq1k_2k)
print('>2k, <=5k\t',seq2k_5k)
print('>5k, <=10k\t',seq5k_10k)
print('>10k, <=20k\t',seq10k_20k)
print('>20k, <=30k\t',seq20k_30k)
print('>30k, <=40k\t',seq30k_40k)
print('>40k, <=50k\t',seq40k_50k)
print('>50k\t',seq50kplus)
