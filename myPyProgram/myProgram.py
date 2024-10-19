import sys
import json
import time
import re
import numpy as np
#import pandas as pd
import scipy
import Bio
import math
import random
import statistics
#import suffix_trees
from Bio import SeqIO
from Bio.Data import CodonTable
from scipy import stats
#from scipy.stats import chi2
#from scipy import optimize
#import scikit_posthocs as sp # this library is required for post hoc process after Kruskal Wallis tests (Dunn)
#import matplotlib as mpl
#import matplotlib.pyplot as plt
# The path for the datafiles and input files are respectively :
# /home/users/m/j/mjoiret/myData
# /home/users/m/j/mjoiret/myInput
# These paths are defined in the bash script to parse the right filename with the right argument.
# check the number of arguments that were parsed:
if len(sys.argv) < 6:
    print('Usage: python' + sys.argv[0] + 'filenames')
    exit()
if len(sys.argv) > 6:
    print('Error: too many arguments added after the python program')
    exit()
# Get the name of the 3 data files, the 1 variable parsed input file and the name of the outputfile:
datafile1 = str(sys.argv[1]) # this is the JSON file with Yeasts elongation rates for all sense codons.
datafile2 = str(sys.argv[2]) # this is the tunnel electrostatics profile for the axial forces.
datafile3 = sys.argv[3] # this is the fasta file with the transcripts sequences.
inputfile = sys.argv[4] # this is the reads counts and initiation rates for all transcripts.
outputfile = sys.argv[5] # this provides the name of the output file
                         # produced by each job(or task) run by the SLURM manager.
                         # Each input file produces a single output file (one to one mapping)

# read the JSON file with the rates for elongation cycle for Yeasts:
dataJSON = open(datafile1, 'r')
# Note you do not want to print in the general log file of the bash script !
# Yous should print in a dedicated file your specifics milestones of interest !
# You need do define a logical flag for a legit print...to write in the right file
print('--------------------------------------------------------------------')
print('Dictionnary of elongation rates parameters for all codons in Yeasts:')
print('--------------------------------------------------------------------')
dataRatesDict = json.load(dataJSON)
print(dataRatesDict)
dataJSON.close()

#read the JSON file with the electrostatic profile for the axial forces in the ribosome exit tunnel:
dataJSON = open(datafile2, 'r')
print('--------------------------------------------------------------------')
print('Dictionnary of axial forces positions:')
print('--------------------------------------------------------------------')
electrostaticAFDict = json.load(dataJSON)
print(electrostaticAFDict)
dataJSON.close()

#read the fasta file with the transcript sequences:
transcriptsSeq = open(datafile3, 'r')
#close the file:
transcriptsSeq.close()

#read the reads counts and the initiation rates for all transcripts in the correct order
readsStat = open(inputfile, 'r')
transcrID = list()
readCount = list()
initRate = list()
line = readsStat.readline()
lineString = line.rstrip('\n')
lineCount = 0
while line != '':
    col = lineString.split()
    transcrID.append(str(col[0]))
    readCount.append(int(col[1]))
    initRate.append(float(col[2]))
    lineCount += 1
    line = readsStat.readline()
    lineString = line.rstrip('\n')
#close the input file:
readsStat.close()

#find the init rate for '(Bgn)':
indexFound = transcrID.index('(Bgn)')

# produce the result file:
headerLine = 'mRNA nb: ' + '\t' +  str(int(lineCount)) + '\n'
lineResult = '(Bgn)' + '\t' + 'init rate = ' + '\t' + format(initRate[indexFound], '6.4f') + '\n'
resultFile = open(outputfile, 'w')
resultFile.write(headerLine)
resultFile.write(lineResult)
resultFile.close()

finalFile = open('result.txt', 'a') # append to the end of the file
headerLine = 'next result from next input file:' + '\n'
finalFile.write(headerLine)
finalFile.write(lineResult)
finalFile.close()
#sys.exit()
