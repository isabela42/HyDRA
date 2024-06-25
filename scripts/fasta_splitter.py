#!/usr/bin/python
# coding=utf-8

'''
Written by Isabela Almeida
Created on Aug 28, 2023
Last modified on Nov 09, 2023
Version: 1.0.6
 
Description: Takes a FASTA file, counts the bases in each
sequence, and outputs splits sequences into different files
according to their nucleotide lenght.

User inputs:
-f --Fasta      FASTA file
-n -numberbases n value to split sequences into bigger than or up to n files
-p --pathtostem Path/to/stem of file to split sequences into

Outputs file | description 
{pathtostem}_{numberbases}le.fasta    | sequences with lenght less than or equal to -n
{pathtostem}_{numberbases}plus.fasta  | sequences with lenght -n plus
{pathtostem}_20-30k.fasta             | sequences with lenght >20k & <= 30k
{pathtostem}_30-40k.fasta             | sequences with lenght >30k & <= 40k
{pathtostem}_40-50k.fasta             | sequences with lenght >40k & <= 50k
{pathtostem}_50kplus.fasta            | sequences with lenght >50k
{pathtostem}_counts.txt               | sequence ID and lenght for each input sequence  

<Envinroment: VSC>                                   
<OS: Linux (Ubuntu)>                                          
'''

from pandas import DataFrame
import pandas as pd
from Bio import SeqIO
import argparse
import sys

# Get command-line args
userInput = argparse.ArgumentParser(description='Takes a FASTA file, counts the bases in each '
                                                'sequence, and outputs splits sequences into '
                                                'different files according to their nucleotide '
                                                'lenght.')
requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--Fasta', action='store',
                          required=True, type=argparse.FileType('r'),
                          help='Target FASTA file')
requiredNamed.add_argument('-n', '--numberbases', action='store',
                      required=True, type=int, default='20000',
                      help='number of bases to separate sequences into different files.'
                      'E.g. 20k bp should be 20000')
requiredNamed.add_argument('-p', '--pathtostem', action='store',
                           required=True,
                           type=str, help='Path/to/stem of file to split sequences into')

## Import user-specified command line values.
args = userInput.parse_args()
fasta_file = args.Fasta
n_value = args.numberbases
out_path_stem = args.pathtostem

## Create output file
out_20kplus = "{0}_20-30k.fasta".format(out_path_stem)
out_30kplus = "{0}_30-40k.fasta".format(out_path_stem)
out_40kplus = "{0}_40-50k.fasta".format(out_path_stem)
out_50kplus = "{0}_50kplus.fasta".format(out_path_stem)
out_up_to_n = "{0}_{1}le.fasta".format(out_path_stem,n_value)
out_n_plus = "{0}_{1}plus.fasta".format(out_path_stem,n_value)
out_counts = "{0}_counts.txt".format(out_path_stem)

file_20kplus = open(out_20kplus, "w+" )
file_30kplus = open(out_30kplus, "w+" )
file_40kplus = open(out_40kplus, "w+" )
file_50kplus = open(out_50kplus, "w+" )
file_up_to_n = open(out_up_to_n, "w+" )
file_nplus = open(out_n_plus, "w+" )
file_counts = open(out_counts, "w+" )


## Get Target fasta sequences
target_sequences = SeqIO.parse(fasta_file,"fasta")

## Loop through each fasta sequence
n_plus = up_to_n = num = 0
seq0_1k = seq1k = seq1k_2k = seq2k_5k = seq5k_10k = seq10k_20k = seq20k_30k = seq30k_40k = seq40k_50k = seq50kplus = 0

for fasta in target_sequences:
    num += 1
    seqName, sequence = fasta.description, str(fasta.seq)
    thisFASTA = []
    seqLen = len(sequence)

    if num == 1:
        longest = smallest = seqLen
    
    file_counts.write("> {0}\n".format(seqName))
    file_counts.write("{0}\n".format(seqLen))

    # count lenghts
    if seqLen >= longest: longest = seqLen
    if seqLen <= smallest: smallest = seqLen
    
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

    # write sequences to files
    if seqLen > n_value:
        n_plus += 1
        file_nplus.write("> {0}\n".format(seqName))
        file_nplus.write("{0}\n".format(sequence))
    else:
        up_to_n += 1
        file_up_to_n.write("> {0}\n".format(seqName))
        file_up_to_n.write("{0}\n".format(sequence))

    if seqLen > 50000:
        file_50kplus.write("> {0}\n".format(seqName))
        file_50kplus.write("{0}\n".format(sequence))
    elif seqLen > 40000 and seqLen <= 50000:
        file_40kplus.write("> {0}\n".format(seqName))
        file_40kplus.write("{0}\n".format(sequence))
    elif seqLen > 30000 and seqLen <= 40000:
        file_30kplus.write("> {0}\n".format(seqName))
        file_30kplus.write("{0}\n".format(sequence))
    elif seqLen > 20000 and seqLen <= 30000:
        file_20kplus.write("> {0}\n".format(seqName))
        file_20kplus.write("{0}\n".format(sequence))  

# Print info about the results to terminal
print('Sequences\tValue')
print('n\t', n_value)
print('Total\t', num)
print('Longest\t', longest)
print('Smallest\t', smallest)
print('n(', n_value, ')+\t', n_plus)
print('<=n(', n_value, ')\t', up_to_n)
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