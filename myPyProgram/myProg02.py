import sys
import json
import time
#import timer
import re
import numpy as np
#import pandas as pd
import scipy
import Bio
import math
import random
import statistics
import suffix_trees
#from timer import Timer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from scipy import stats
from suffix_trees import STree
#import dependencies for classes, methods and functions (OOP) to be used in this main program:
sys.path.append('/home/users/m/j/mjoiret/myPyPrograms')
import myClasses
from myClasses import time_millis
from myClasses import decode
from myClasses import indexOfRead
from myClasses import pdf_HYPOEXP
from myClasses import expRandomTimeSample
from myClasses import hypoExpRandomTime
from myClasses import hypoExpRandomSample
from myClasses import ADATsilencedDict # this could be commented out.

# The path for the datafiles, input files, program files and output files are respectively :
# /home/users/m/j/mjoiret/myData
# /home/users/m/j/mjoiret/myInput
# /home/users/m/j/mjoiret/myPyPrograms
# /home/users/m/j/mjoiret/myoutput
# These paths are defined in the bash script to parse the right filename with the right argument.
# check the number of arguments that were parsed:
if len(sys.argv) < 8:
    print('Usage: python' + sys.argv[0] + 'filenames')
    exit()
if len(sys.argv) > 8:
    print('Error: too many arguments added after the python program')
    exit()
# Get the name of the 3 data files, the 1 variable parsed input file, the name of the outputfile
# and extra-argument used for the simulation:
datafile1 = str(sys.argv[1]) # this is the JSON file with Yeasts elongation rates for all sense codons.
datafile2 = str(sys.argv[2]) # this is the tunnel electrostatics profile for the axial forces.
datafile3 = sys.argv[3] # this is the fasta file with the transcripts sequences.
inputfile = sys.argv[4] # this is the reads counts and initiation rates fold change wrt reference for all transcripts.
outputfile = sys.argv[5] # this provides the name of the output file
                         # produced by each job(or task) run by the SLURM manager.
                         # Each input file produces a single output file (one to one mapping)
riboRatio = float(sys.argv[6]) # this is 6th argument parsed that will determine the ribosome pool ratio used in the simulation
init_rate_SLURM = float(sys.argv[7]) # this is the 7th argument parsed that will determine the reference initiation rate of
# any ribosome on any standard transcript (if the pick up rate of any transcript was the same across all transcripts).
#
# How do I get the SLURM array task id or the task number in this instance of the python program?
# slicing or regex capture from outputfile='ratio${i}res$SLURM_ARRAY_TASK_ID.txt':
# or ---> outputfile='ratio${i}init${j}res$SLURM_ARRAY_TASK_ID.txt':
task_id_match = re.search('ratio(\d+)init\d+res.*', outputfile)
task_id = int(task_id_match.group(1))
job_array_index_match = re.search('\w+\d+\w+\d+res(\d+)\.txt', outputfile)
job_array_index = int(job_array_index_match.group(1))
# to comment out (4 lines below):
print('ribosome ratio ', format(riboRatio, '.2f'))
print('initiation rate ', format(init_rate_SLURM, '.4e'))
print('task_id_match ', task_id_match)
print('task_id ', task_id)
print('job_array_index_match ', job_array_index_match)
print('job_array_index ', job_array_index)

# read the JSON file with the rates for elongation cycle for Yeasts:
dataJSON = open(datafile1, 'r')
dataRatesDict = json.load(dataJSON)
print('dataRatesDict: is below')
print(dataRatesDict)
dataJSON.close()

#read the JSON file with the electrostatic profile for the axial forces in the ribosome exit tunnel:
dataJSON = open(datafile2, 'r')
electrostaticAFDict = json.load(dataJSON)
dataJSON.close()

#read the fasta file with the transcript sequences:
transcriptsSeq = open(datafile3, 'r')
# scan sequentielly all the transcripts and retrieve the cds sequence, the mRNA sequences
# from start codon to sto codon, count and check that the mRNA length is a integer multiple
# of 3 nucleotides (codon), count the codons and get a codon list per transcript,
# for all transcripts

#close the file:
transcriptsSeq.close()
flag4allInputFiles = 'true'
if flag4allInputFiles:
    numFasta = 0
    IDtagList = []
    descripList = []
    orfFlagList = []
    orfList = []
    lastCodonList = []
    mRNAlengthList = []

    for cds in SeqIO.parse(datafile3, "fasta"):
        numFasta += 1 # count the transcripts in the fasta sequence
        GeneMatch = re.search('(\(.*\))', cds.description)
        GeneDescript = re.search('(.*\(.*\))', cds.description) #'Description:.*\(.*\)'
        if GeneMatch is None:
            Gid = 'NA'
            GidDescript ='NA'
        else:
            #print('GeneMatch.group(0)=', GeneMatch.group(0))
            Gstr = GeneMatch.group(0)
            Gid = Gstr[1:len(Gstr)-1]
            if GeneDescript is None:
                GidDescript = 'description not read'
            else:
                Descript = GeneDescript.group(0)
                GidDescript = Descript[0:len(Descript)-len(Gstr)]
        IDtagList.append(Gid)
        descripList.append(GidDescript)
        mRNA = str(cds.seq.transcribe())
        mRNAlength = len(mRNA)
        st = STree.STree(mRNA)
        start_codon = 'AUG'
        stop_codon = ['UAA', 'UAG', 'UGA']
        STARTnucleotideIndex = st.find(start_codon)
        # Translated sequence starts at START codon and stops at STOP codon in a cds sequence:
        last_codon = mRNA[(mRNAlength-3):mRNAlength]
        lastCodonList.append(last_codon)
        TranslatedSeq = mRNA[0:mRNAlength]
        # The Open Reading Frame ORF is thus mRNAorf
        mRNAorf = TranslatedSeq[0:mRNAlength]
        orfLength = len(mRNAorf)
        mRNAlengthList.append(orfLength)
        # check for 3 tests : whole number of triplets, start and end codons correct, length > 90 nucleotides
        # if the three conditions are correct: flag = 1, else: flag = 0
        flag = 1
        if orfLength%3 == 0:
            CodonsNumb = orfLength//3
        else:
            flag = 0
            FailedTripletsCount += 1
            CodonsNumb = orfLength//3

        if (last_codon in stop_codon) & (STARTnucleotideIndex==0):
            #print('Start and stop codons are correct: last codon is: ', last_codon)
            flag4allInputFiles = 'true'
        else:
            flag = 0
            FailedStartStopCount += 1
        if mRNAlength > 90:
            #print('mRNA length larger than 90: ', mRNAlength)
            flag4allInputFiles = 'true'
        else:
            flag = 0
            FailedLengthcount += 1
        if flag == 1:
            orfFlagList.append('yes')
        else:
            orfFlagList.append('no')
        ORF = cds.seq.transcribe()[STARTnucleotideIndex: (STARTnucleotideIndex + orfLength)]
        orfList.append(ORF)

# Forcing a time lag of 3 to 5 seconds between the parallel multiprocessing of jobs or tasks
time.sleep(3*task_id+job_array_index*2) #pauses the current thread

# time in seconds elapsed since epoch reference at the start of this main:
started_time = time.time()
print('started time of this thread in ms =', format(started_time*1000, '.3f'))

# only print once in the job log file:
# print only once in the log file the data files info the list of the transcripts fasta file:
print('outputfile name=', outputfile)
#print('sys.argv[5]= ', sys.argv[5]) # /home/users/m/j/mjoiret/myOutput
if outputfile == '/home/users/m/j/mjoiret/myOutput/ratio0init0res0.txt':
    print('--------------------------------------------------------------------')
    print('transcript file name = ', datafile3)
    print('elongation rates per codon file name = ', datafile1)
    print('electrostatic profile of axial forces in exit tunnel file name = ', datafile2)
    print('reads counts and initiation rate first input file name =', inputfile)
    print('--------------------------------------------------------------------')
    print('name of the first output file = ', outputfile)
    print('--------------------------------------------------------------------')
    print('Dictionnary of elongation rates parameters for all codons in Yeasts:')
    print('--------------------------------------------------------------------')
    print(dataRatesDict)
    print('--------------------------------------------------------------------')
    print('Examples of retrievieng the three rates for codon GCG:')
    print('GCG step 1 rate = ', dataRatesDict['GCG'][0])
    print('GCG step 2 rate = ', dataRatesDict['GCG'][1])
    print('GCG step 3 rate = ', dataRatesDict['GCG'][2])
    print('--------------------------------------------------------------------')
    print('Dictionnary of axial forces positions:')
    print('--------------------------------------------------------------------')
    print(electrostaticAFDict)
    print('--------------------------------------------------------------------')
    print('list of the transcripts in the fasta format data file:')
    for i in range(numFasta):
        print(IDtagList[i] + '\t'+ descripList[i] +'\t' + orfFlagList[i] +'\t' + lastCodonList[i] +'\t' + str(mRNAlengthList[i]) +'\n')
        print(orfList[i])

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
print('--------------------------------------------------------------------')
print('input file for RNA-Seq reads counts and initiation rates: ', inputfile)
for i in range(numFasta):

    print(format(i+1, '3d') + ' ' + '{: <16}'.format(IDtagList[i]) + '==' + transcrID[i] + '\t' +\
    format(readCount[i], '4.2f') + '\t' + format(initRate[i], '4.2f'))
print('--------------------------------------------------------------------')

# produce a global file collecting meta-results for all jobs (multiprocessing the same python program)
finalFile = open('result.txt', 'a') # append to the end of the file
headerLine = '---------------------------------------------------------------' + '\n'
finalFile.write(headerLine)
headerLine = 'input file processed:' + '\t' + str(inputfile) +'\n'
finalFile.write(headerLine)
headerLine = '---------------------------------------------------------------' + '\n'
finalFile.write(headerLine)
headerLine = 'mRNA read counts and initiation rates for this input file:' +'\n'
finalFile.write(headerLine)
headerLine = 'List of transcripts variants' +'\n'
finalFile.write(headerLine)
headerLine = '---------------------------------------------------------------' + '\n'
finalFile.write(headerLine)
headerLine = '#' + '\t' + 'IDtags' + '\t' +'transcripts variants' + '\t' + 'ORF' +'\t'+ \
'last codon' +'\t' + 'peptide length' + '\t' + 'mRNA length' + '\t' + 'read copy' + '\t' + 'init. rate'+'\n'
finalFile.write(headerLine)

# number of transcripts that were read in fasta file:
for i in range(numFasta):
    lineResult = format(i+1, '3.0f') + ' ' + '{: <6}'.format(IDtagList[i]) +'\t' +\
    '{: <85}'.format(descripList[i]) + '\t' + str(orfFlagList[i]) +\
    '\t'+ str(lastCodonList[i]) + '\t' + str(int((mRNAlengthList[i]-3)//3)) + '\t' +\
    str(int(mRNAlengthList[i])) +\
    '\t' + format(readCount[i], '4.2f') + '\t' + format(initRate[i], '4.2f') + '\n'
    finalFile.write(lineResult)
# last line (check)
lastLine = 'The number of transcripts in fasta file was ' + format(numFasta, '3d') +\
' and number of RNA-seq reads type was '+ format(len(readCount), '3d') +'.' +'\n'
finalFile.write(lastLine)
headerLine = '---------------------------------------------------------------' + '\n'
finalFile.write(headerLine)

#finalFile.close()

# -----------------------------------------------------------------------------
# Functions: set-up
# -----------------------------------------------------------------------------
def riboCountOnTranscript(thisReadTranscriptID):
    # This function returns the number of ribosomes that are currently elongating on the
    # peculiar read of the transcript copy object (i.e., a lattice instance) given in argument.
    riboCountOnTranscript = 0
    for i_rib in range(len(ribosome)):
        if ribosome[i_rib].uniqueReadID == thisReadTranscriptID.get_uniqueID():
            riboCountOnTranscript += 1
    return riboCountOnTranscript

def riboTypeCounts(riboList):
    # counts the ribosome types of ribosomes. The argument riboList is the list of
    # the ribosoes objects that are currently instantiated and are either 'free' or 'initiated' or 'translating'.
    # This function returns as tuple the tally of the free ribosomes and the tally of
    # the initiated ribosomes:
    freeRibCount = 0
    initiatedRibCount = 0
    for i_rib in range(len(riboList)):
        if riboList[i_rib].type == 'free':
            freeRibCount += 1
        if ((riboList[i_rib].type == 'initiated') or (riboList[i_rib].type == 'translating')):
            initiatedRibCount += 1
    return [freeRibCount, initiatedRibCount]

def maxNbRibo(latticeList):
    # This function scans all the transcripts reads and retrieve the current number of
    # elongating ribosomes on each. It returns the largest number obtained across all transcriptome:
    maxNriboPerRead = 0
    for i_read in range(len(latticeList)):
        if latticeList[i_read].get_rpfCount() > maxNriboPerRead:
            maxNriboPerRead = latticeList[i_read].get_rpfCount()
    return maxNriboPerRead

def constrain(val, min_val, max_val):
    # This function is a python way to constrain a value as it
    # would be done in javascript p5.js processing
    if val < min_val:
        val = min_val
    elif val > max_val:
        val = max_val
    return val

def boolToDigit(status):
    # this function receives in argument aboolean value (True or False) and
    # returns 1 if True or 0 if False:
    if status:
        return 1
    else:
        return 0

# ----------------------------------------------------------------------------------
# Factors affecting the elongation kinetics (each factor is a boolean entry):
# ----------------------------------------------------------------------------------
# Proline factor and EF-P or eIF5A elongation factor depletion affects
# the peptide bond formation rate (rate 2):
prolineStatus = True
EFPdepletedStatus = False
# Ribosome exit tunnel electrostatic interaction also affects rate 2:
tunnelStatus = True
# ADAT2 silencing mainly affects tRNA accommodation and proofreading (rate 1)
# for 37 codons / anticodons out of the 61 (8 amino acid residues TAPS LIVR):
ADAT2silencedStatus = False
# ELP3-URM1 silencing affects tRNA accommodation and proofreading (rate 1)
# for 6 codons / anticodons out of the 61 (3 amino acid residues K, Q, E):
U34hypoModificationStatus = False # experiment 1
# mRNA secondary structure:
mRNAsecondaryStrucStatus = False # not implemented yet

# ----------------------------------------------------------------------------------
# setup and instantiations of objects:
# ----------------------------------------------------------------------------------
# Instantiation of transcripts: assigning the transcripts to their objects'attributes
# -----------------------------------------------------------------------------------
transcriptCardinality = len(orfList)
# Define an array of lattices = array of transcripts objects
lattice = []
# We initiate 'riboRecruitScore' to generate a multinomial random sampling of the transcripts
# for recruitment of the ribosomes for transcript initiation relative probabilities (not a uniform
# probability):
indicesListToSampleFrom = []
# is the list of transcript id number to sample from. In this list, the same transcript id
# number will appear proportionally to the number of recruitment score for initiation.
riboRecruitScore = np.zeros(transcriptCardinality)
# is the array of the readInitRate for each transcript but multiplied by 10 for 10% relative difference.
for i_tr in range(transcriptCardinality):
    # retrieve the initiation rate 'initRate[i_tr]' per transcript:
    riboRecruitScore[i_tr] = int(10*initRate[i_tr])
    # riboRecruitScore.append(int(10*initRate[i_tr]))
    # scaled by a factor of 10 to take into account only the first decimal digit of initiation rates.
    # retrieve the number of reads copy per transcript 'readCount[i_tr]':
    for i_read in range(readCount[i_tr]):
        read_str = str(IDtagList[i_tr])+'#'+str(i_read)
        for i_score in range(int(riboRecruitScore[i_tr])):
            indicesListToSampleFrom.append(indexOfRead(IDtagList, readCount, read_str))
        # lattice objects initiations:
        # Below, we add the transcript copies if there are copies:
        # transcLength, transcSequence, transcCopy, transcName, transcNumID)
        lattice.append(myClasses.Lattice(mRNAlengthList[i_tr]/3, orfList[i_tr], readCount[i_tr], IDtagList[i_tr], i_read))
    # end of lattice objects initiations.
# end of outer for loop
# calculate the full transcripts cardinality:
transcriptFullCardinality = len(lattice)
# Initialize a 'pseudo-matrix that will record the (cumulated) queueing time set points for each codon of each transcript read copy'
# Initialize a 'pseudo-matrix that will record the (cumulated) real residence time for each codon of each transcript read copy'
# Initialize a 'pseudo-matrix that will record the (cumulated) ribosome footprint occurence for each codon of each transcript read copy'
RibosSeqSPtruth = []
RibosSeqRRTtruth = []
RibosSeqFCsampled = []
# compute the total length of the instantiated transcriptome:
transcriptomeLength_nts = 0
for i in range(transcriptFullCardinality):
    transcriptomeLength_nts += lattice[i].get_mRNAlength()
    RibosSeqSPtruth.append(np.zeros(lattice[i].get_length())) # initialize second dimension of array
    RibosSeqRRTtruth.append(np.zeros(lattice[i].get_length())) # initialize second dimension of array
    RibosSeqFCsampled.append(np.zeros(lattice[i].get_length())) # initialize second dimension of array

transcriptomeLength_codons = transcriptomeLength_nts/3
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# Instantiation of ribosomes: assigning ribosomes' objects'attributes:
# -----------------------------------------------------------------------------------
# ribosome footprint length (nt scale):
riboFootprintLength = 33 # 33 nucleotides = 33 pixels = 11 Codons
# Initial number (at program start) of free ribosomes
# (results from ribosome pool ratio times transcriptFullCardinality)
freeRibosomeInitialNumber = int(transcriptFullCardinality * riboRatio)

# Define an array of ribosomes = array of ribosomes objects
ribosome = []
# At the very start, all the ribosomes are free ribosomes.
# How many ribosomes are there ? The answer is the number of transcript ID times the ribosome pool factor
# that was calculated just above, i.e., freeRibosomeInitialNumber:
# when instantiated, all ribosomes should be just static. But, when initiated on a given transcript,
# it receives a peculiar initiation rate fold change associated to this transcript (multinomial prob. distr.).
# The reference initiation rate comes from a SLURM parsed argument, called init_rate_SLURM.
for i_rib in range(freeRibosomeInitialNumber):
    ribosome.append(myClasses.Ribosome(started_time))
    ribosome[i_rib].set_initiationRate(init_rate_SLURM)

# ----------------------------------------------------------------------------------
# START MAIN LOOP : time is running...(similar to draw function in p5.js processing)
# ----------------------------------------------------------------------------------
#time_out = 2.5 # 2.5 hours
time_out = 1.9 # (this is 114 minutes )#0.05  (this is 3 minutes)
time_flag = False

replicateCount = 0
# find the index for '(EIF5A)':
indexFound = IDtagList.index('EIF5A')
indexOfReadFound = indexOfRead(IDtagList, readCount, 'EIF5A#0')
# Number of time points required:
time_point_nb = 24
time_flag_arr = [False for i in range(time_point_nb)] # naive list comprehension
time_point_item = np.linspace(0, time_point_nb, time_point_nb, endpoint=False)
spare_time = 20 * 1000 # 20 seconds = 20,000 ms
start_sample_time_ms = spare_time
stop_sample_time_ms = (1.8 * 60 * 60 * 1000 + spare_time) - spare_time
time_point_t = np.linspace(start_sample_time_ms, stop_sample_time_ms, time_point_nb, endpoint=True)
time_array_min = np.zeros((time_point_nb, 2))
# Number of snapshots required to build a ribosome footprint density map:
time_lap_between_snaps = 10 * 1000 # 10 seconds = 10,000 ms = 10 * 1000
snapshots_nb = math.floor(stop_sample_time_ms/time_lap_between_snaps) # if 1.8 hours with 10 s lap, this is 648 snaps
eventSnaps_arr = [False for i in range(snapshots_nb)]
time_snap_item = np.linspace(0, snapshots_nb, snapshots_nb, endpoint=False)
time_snap_t = np.linspace(spare_time/2, stop_sample_time_ms + spare_time/2, snapshots_nb, endpoint=True)
footprintedFragmentsCount = 0
# define and initiate important array data types:
# chronology of free and initited ribosomes:
chronoRiboFreeInitiated = np.zeros((time_point_nb, 2))
# chronology of avg distance between two ribosomes across all transcriptome:
chronoAvgDistBetweenRibinNTS = np.zeros(time_point_nb)
# chronology of protein abundance per transcript:
chronoProtCountForTr = np.zeros((time_point_nb, transcriptCardinality))
# chronology of RPF (number of ribosome footprint fragment per transcript):
chronoRPFcountForTr = np.zeros((time_point_nb, transcriptCardinality))

# define and initiate codon_QueueTime_Mat:
# numpy zeros: (# nb of codons in the transcript in row dimension and # nb of time points in column dimension)
col_nb = time_point_nb
#row_nb = lattice[indexOfReadFound].get_length()
#codon_QueueTime_Mat = np.zeros((row_nb, col_nb,))
codon = []
i_time = 0
i_snap = 0
maxPolysomeNb = 0
maxNbPoly = 20
# chronology of polysome fragment profiling across all transcripts:
chronoPolysomeFragProfile = np.zeros((time_point_nb, maxNbPoly+1)) # warning this will be a list of arrays with a variable nb of col
# snapshots chronology of polysome fragmentation profiling per transcript:
chronoPolyFragProfileForTr = np.zeros((snapshots_nb, transcriptCardinality, maxNbPoly + 1))
polyFragProfileForTr = np.zeros((transcriptCardinality, maxNbPoly + 1)) # will be aggregated over all snapshots

while time_millis(started_time)/(1000.0 * 60 * 60) <= time_out:
    # identation while until time-out mainloop (loop starts here):
    # -- start Monte Carlo step recapitulates like this: ***
    # -- scan all ribosomes on all transcripts: (2 for loops: outer + inner nested loop)
    # -- outer loop:
    for i_rib in range(freeRibosomeInitialNumber):
    # .. pick a transcript randomly following the multinomial distribution in the transcript ID category
        i_multinomial_random = random.randint(0, len(indicesListToSampleFrom) - 1)
        picked_transcript = int(indicesListToSampleFrom[i_multinomial_random])
        # .. initiate a ribosome on this randomly picked transcript:
        #ribosome[i_rib].initiates(lattice[picked_transcript], started_time)
    # .... inner for loop:
        for i_tr in range(transcriptFullCardinality):
            currentTime = time_millis(started_time)

    # .. pick a transcript randomly following the multinomial distribution in the transcript ID category
            #i_multinomial_random = random.randint(0, len(indicesListToSampleFrom) - 1)
            #picked_transcript = int(indicesListToSampleFrom[i_multinomial_random])

    # .. initiates a ribosome on this randomly picked transcript:
            if ribosome[i_rib].type == 'free':
                ribosome[i_rib].initiates(lattice[picked_transcript], started_time)

    # .. attempt new translocations:
            if ribosome[i_rib].type != 'free':
                # this ribosome instance can only translocates if the transcript it is on
                # is the correct transcript currently scanned (this condition is checked in the ribosome
                # translocates method). For sanity, the condition is checked again here:
                if ribosome[i_rib].uniqueReadID == lattice[i_tr].get_uniqueID():
                    ribosome[i_rib].translocates(lattice[i_tr], dataRatesDict, electrostaticAFDict, started_time, IDtagList, readCount, RibosSeqSPtruth, RibosSeqRRTtruth, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus)

    # .... end inner for loop
    # -- end outer for loop
    # -- end Monte Carlo recapitulating steps  ***
    # -- update the number of produced proteins, i.e., aggregate the read copies translated for each
    # -- transcript ID (name):
    translated_count_list = []
    for i_tr in range(transcriptCardinality):
        cumul_translated = 0
        for i_read in range(transcriptFullCardinality):
            if lattice[i_read].get_name() == IDtagList[i_tr]:
                cumul_translated += lattice[i_read].get_readTranslated()
        translated_count_list.append(cumul_translated)
    # -- update the number of proteins produced for each transcript name in each transcript read:
    for i_tr in range(transcriptCardinality):
        for i_read in range(transcriptFullCardinality):
            if lattice[i_read].get_name() == IDtagList[i_tr]:
                lattice[i_read].set_countTranslated(translated_count_list[i_tr])
    # -- extract the polysome fragment current profile:
    # -- i.e., the frequency distribution of the number of reads having a number of
    # -- ribosomes equal to n, with i_n_some=0, 1, 2, 3, ... maxNbRibo.
    currentMaxNbRibInRead = maxNbRibo(lattice)
    if currentMaxNbRibInRead > maxPolysomeNb:
        maxPolysomeNb = currentMaxNbRibInRead
    #maxNbPoly = 20
    polysomeFragProfile = np.zeros(maxNbPoly + 1) # maxNbPoly + 1
    polysomeFragProfileOnRead = []
    # -- update the number of ribosome protected fragments on each transcript
    for i_read in range(transcriptFullCardinality):
        polysomeFragProfileOnRead.append(np.zeros(maxNbPoly+1))
        lattice[i_read].set_rpfCount(riboCountOnTranscript(lattice[i_read]))
        for i_n_some in range(maxNbPoly + 1):
            if lattice[i_read].get_rpfCount() == i_n_some:
                polysomeFragProfile[i_n_some] += 1
                polysomeFragProfileOnRead[i_read][i_n_some] += 1
    # determine the current number of free ribosomes and actively translating ribosomes:
    n_ribo_freetype = riboTypeCounts(ribosome)[0]
    n_ribo_translating = riboTypeCounts(ribosome)[1]
    # determine average distance between two ribosomes transcriptome wide:
    if n_ribo_translating > 0:
        avg_dist_between_ribo_inCodons = float(format(transcriptomeLength_codons/n_ribo_translating, '.2f'))
        avg_dist_between_ribo_inNts = float(format(transcriptomeLength_nts/n_ribo_translating, '.2f'))
    else:
        avg_dist_between_ribo_inCodons = float(format(transcriptomeLength_codons, '.2f'))
        avg_dist_between_ribo_inNts = float(format(transcriptomeLength_nts, '.2f'))
    # format the current value of the clock (elapsed time since simulation start time): min & sec:
    current_simul_time = time_millis(started_time)
    current_simul_t_sec = current_simul_time/1000
    elapsed_simul_min = math.floor(current_simul_time/(1000 * 60))
    elapsed_simul_sec = int((current_simul_time/1000) % 60)
    elapsed_time_txt = format(elapsed_simul_min, '3d') + ' min ' + format(elapsed_simul_sec, '2d') + ' sec.'

    if ((not time_flag) and (time_millis(started_time)/(1000.0 * 60 * 60)>= time_out/2.0)):
        time_flag = True
        print('Half-time of timeout was just reached at ', \
        format(time_millis(started_time)/(1000.0 * 60), '.3f'), ' minutes.')
    #
    if ((not time_flag_arr[i_time]) and (time_millis(started_time) >= time_point_t[i_time])):
        time_flag_arr[i_time] = True
        # track and save result data ... here if you want more
        i_time += 1
        #replicateCount += 1
        if i_time >= time_point_nb:
            break
        # track and save result data at time point:
        time_array_min[i_time - 1] = [elapsed_simul_min, elapsed_simul_sec]
        #++++++++++++++++++++++++
        # save the chronology of:
        # +++++++++++++++++++++++
        # - ribosomes free and initiated:
        chronoRiboFreeInitiated[i_time - 1][0]= n_ribo_freetype
        chronoRiboFreeInitiated[i_time - 1][1]= n_ribo_translating
        # - average distance in nts between two ribosome across all transcriptome:
        chronoAvgDistBetweenRibinNTS[i_time -1] = avg_dist_between_ribo_inNts
        # - protein abundance per transcript ID (aggregated over all its read copies):
        # - number of ribosome protected fragment (per transcript ID):
        for i_tr in range(transcriptCardinality):
            chronoProtCountForTr[i_time -1][i_tr]= translated_count_list[i_tr]
            cumulRPFacrossReadCopies = 0
            for i_read in range(transcriptFullCardinality):
                if lattice[i_read].get_name() == IDtagList[i_tr]:
                    cumulRPFacrossReadCopies += lattice[i_read].get_rpfCount()
            # (divide by the number of read copies for this transcript,i.e., by readCount[i_tr])
            chronoRPFcountForTr[i_time -1][i_tr] = float(format(cumulRPFacrossReadCopies/readCount[i_tr], '5.2f'))
        # - polysome fragment profiles:
        for i_some in range(len(polysomeFragProfile)):
            chronoPolysomeFragProfile[i_time - 1][i_some] = polysomeFragProfile[i_some]

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Accumulating ribosome footprint count at codon resolution from systematic Ribosome
    # occupancy snapshots to build Ribo-seq. Each and every time_lap, update and sum the
    # ribosome count and footprint count on the transcripts reads. Do this till the end of
    # the simulation.
    # Aggregate per transcript id over the read copies
    if ((not eventSnaps_arr[i_snap]) and (time_millis(started_time) >= time_snap_t[i_snap])):
        eventSnaps_arr[i_snap] = True
        i_snap += 1
        if i_snap >= snapshots_nb:
            break
        # take a snapshot of the ribosomes footprint, cumulated count of the RPFs,
        # update and aggregate RPF for all transcripts id:
        for i_read in range(transcriptFullCardinality):
            #if lattice[i_read].get_rpfCount() > 0:
            if ((lattice[i_read].get_currentFootprint()[0]>=0) and (len(lattice[i_read].get_currentFootprint())>0)):
                #footprintedFragmentsCount += lattice[i_read].get_rpfCount()
                footprintedFragmentsCount += len(lattice[i_read].get_currentFootprint())
                # update the footprint count on this read at this snapshot event:
                #for i_rpf in range(lattice[i_read].get_rpfCount()):
                for i_rpf in range(len(lattice[i_read].get_currentFootprint())):
                    # cover all codons of this footprint:
                    lowerLeftCodon = int(lattice[i_read].get_currentFootprint()[i_rpf]/3-4)
                    upperRightCodon = int(lattice[i_read].get_currentFootprint()[i_rpf]/3+6)
                    # constrained by length of the transcript read in codon number:
                    lowerBound = int(constrain(lowerLeftCodon, 0, lattice[i_read].get_length()))
                    upperBound = int(constrain(upperRightCodon, 0, lattice[i_read].get_length()))
                    for i_codon in range(lowerBound, upperBound):
                        RibosSeqFCsampled[i_read][i_codon] += 1
        # take a snapshot of the polysome fragment profile per transcript (aggregate over the reads):
        for i_tr in range(transcriptCardinality):
            first_read_string = str(IDtagList[i_tr]) + '#' + str(0)
            indOfTranscriptFirstRead = indexOfRead(IDtagList, readCount, first_read_string)
            read_indices = [i_ind for i_ind in range(indOfTranscriptFirstRead, indOfTranscriptFirstRead + readCount[i_tr])]
            for i_read in read_indices:
                for i_some in range(maxNbPoly + 1):
                    chronoPolyFragProfileForTr[i_snap - 1][i_tr][i_some] += polysomeFragProfileOnRead[i_read][i_some]

    # identation while until time-out mainloop (loop else branch ends here).

# Aggregate per transcript id over the read copies:
# - Retrieve the indices of the reads of a given transcript index
# - (the list of all indices of the reads corresponding to this given transcript index)
# - Keep in mind the number of reads of this transcript for later normalization across read copies:
# - --> readCount[i_tr]
# Aggregate per transcript id over the read copies:
# - the polysome fragment profile per transcript id (over all snapshots ?)... to be continued...
ribosSeqFCaggregated = []
FC = []
for i_tr in range(transcriptCardinality):
    first_read_string = str(IDtagList[i_tr]) + '#' + str(0)
    indOfTranscriptFirstRead = indexOfRead(IDtagList, readCount, first_read_string)
    ribosSeqFCaggregated.append(np.zeros(lattice[indOfTranscriptFirstRead].get_length())) # instantiation of second dimension of this array (=length in codons)
    FC.append(np.zeros(lattice[indOfTranscriptFirstRead].get_length()))
    read_indices = [i_ind for i_ind in range(indOfTranscriptFirstRead, indOfTranscriptFirstRead + readCount[i_tr])]
    for i_codon in range(lattice[indOfTranscriptFirstRead].get_length()):
        for i_read in read_indices:
            ribosSeqFCaggregated[i_tr][i_codon] += RibosSeqFCsampled[i_read][i_codon]
        # normalize the rpf per # of trancript reads:
        FC[i_tr][i_codon] = ribosSeqFCaggregated[i_tr][i_codon]/readCount[i_tr]
# Calculate the translation efficiency (TE) = protein count per nb of transcripts per transcript kb per unit time
# Calculate TE (also not per kb)
TE = []
elapsed_t_in_hour = time_point_t[time_point_nb - 2]/(1000 * 60 * 60)
for i_tr in range(transcriptCardinality):
    first_read_string = str(IDtagList[i_tr]) + '#' + str(0)
    indOfTranscriptFirstRead = indexOfRead(IDtagList, readCount, first_read_string)
    this_kb = float(lattice[indOfTranscriptFirstRead].get_mRNAlength()/1000)
    TE.append(chronoProtCountForTr[time_point_nb - 2][i_tr]/(readCount[i_tr] * this_kb * elapsed_t_in_hour))

# Aggregate the polysome fragmentation profile over all snapshots for each transcript:
polyFragProfileForTr = chronoPolyFragProfileForTr.sum(axis=0)

print('task ', task_id, ' was processed.')
print('job array index ', job_array_index, ' was processed.')
print('ribo ratio = ', riboRatio, 'free ribosomes = ', n_ribo_freetype,'translating ribosomes = ', n_ribo_translating)
print('initiation rate = ', init_rate_SLURM, 'avg dist between ribo (codons) = ', avg_dist_between_ribo_inCodons, 'avg dist (nts) = ', avg_dist_between_ribo_inNts)
print('This instance of the main program was finished after ', \
    format(time_millis(started_time)/(1000.0 * 60), '.3f'), ' minutes.')

# produce the result file:
resultFile = open(outputfile, 'w')
headerLine = '+++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'Parameter space'+'\n'
resultFile.write(headerLine)
headerLine = '+++++++++++++++'+'\n'
resultFile.write(headerLine)
for i in range(numFasta):
    lineResult = format(i+1, '3d') + ' ' + '{: <16}'.format(IDtagList[i]) + '==' + transcrID[i] + '\t' +\
    format(readCount[i], '4.2f') + '\t' + format(initRate[i], '4.2f') +'\n'
    resultFile.write(lineResult)
lineResult = 'transcriptome total size (nts, codons) = ' + format(int(transcriptomeLength_nts), 'd') +'\t'+ format(int(transcriptomeLength_codons), 'd') +'\n'
resultFile.write(lineResult)
lineResult = 'average initiation rate (reference) = ' + format(init_rate_SLURM, '.5e') +'\n'
resultFile.write(lineResult)
lineResult = 'ribosome pool ratio = ' + format(riboRatio, '6.2f') +'\n'
resultFile.write(lineResult)
lineResult = 'transcripts ID number = ' + format(transcriptCardinality, 'd') +'\n'
resultFile.write(lineResult)
lineResult = 'transcript copies number = ' + format(transcriptFullCardinality, 'd') +'\n'
resultFile.write(lineResult)
lineResult = 'ribosome pool size = ' + format(freeRibosomeInitialNumber, 'd') +'\n'
resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
headerLine = '+++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'Activated factors for this simulation run'+'\n'
resultFile.write(headerLine)
headerLine = '+++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
lineResult = 'proline' + '\t' + 'EFP depl' + '\t' + 'tunnel' + '\t' + 'ADAT2 depl' + '\t' +'U34hypo' + '\t' + 'mRNAsecond' +'\n'
resultFile.write(lineResult)
lineResult = str(prolineStatus) +'\t' + str(EFPdepletedStatus) +'\t' + str(tunnelStatus) +'\t' + \
str(ADAT2silencedStatus) + '\t' + str(U34hypoModificationStatus) +'\t' +str(mRNAsecondaryStrucStatus) +'\n'
resultFile.write(lineResult)
# binary boolean format (0 or 1)
lineResult = format(boolToDigit(prolineStatus), 'd') + '\t' + format(boolToDigit(EFPdepletedStatus), 'd') + '\t'+\
format(boolToDigit(tunnelStatus), 'd') + '\t' + format(boolToDigit(ADAT2silencedStatus), 'd') + '\t' +\
format(boolToDigit(U34hypoModificationStatus), 'd') + '\t' + format(boolToDigit(mRNAsecondaryStrucStatus), 'd') + '\n'
resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
# proceed from here with the simulation results...++++***
# record the chronological monitoring of variables:
headerLine = '+++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'History through time points monitoring:  '+'\n'
resultFile.write(headerLine)
headerLine = '+++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
lineResult = 'number of time points = ' + format(time_point_nb-1, 'd') +'\n'
resultFile.write(lineResult)
lineResult = 'time point list= ' + str(time_point_t) +'\n'
resultFile.write(lineResult)
lineResult = '# [free, initiated] ribosomes at time points:'+'\n'
resultFile.write(lineResult)
for line in range(time_point_nb-1):
    lineResult = format(line, 'd')+'\t'+ format(int(chronoRiboFreeInitiated[line][0]), 'd') + '\t' + format(int(chronoRiboFreeInitiated[line][1]), 'd') +'\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
lineResult = '# min  sec  avg_dist_betw_ribosomes (nts):'+'\n'
resultFile.write(lineResult)
for line in range(time_point_nb-1):
    lineResult = format(line, 'd')+'\t'+ format(int(time_array_min[line][0]), 'd') + '\t' + format(int(time_array_min[line][1]), 'd') + '\t' + format(chronoAvgDistBetweenRibinNTS[line], '6.2f') + '\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
lineResult = '# min  sec  protein counts array for each transcript:'+'\n'
resultFile.write(lineResult)
for line in range(time_point_nb-1):
    lineResult = format(line, 'd')+'\t'+ format(int(time_array_min[line][0]), 'd') + '\t' + format(int(time_array_min[line][1]), 'd') + '\t' + str(chronoProtCountForTr[line]) + '\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
lineResult = '# min  sec  RPF count for each transcript:'+'\n'
resultFile.write(lineResult)
for line in range(time_point_nb-1):
    lineResult = format(line, 'd')+'\t'+ format(int(time_array_min[line][0]), 'd') + '\t' + format(int(time_array_min[line][1]), 'd') + '\t' + str(chronoRPFcountForTr[line]) + '\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
lineResult = '# min  sec  nbOfBins(n-some)  polysome fragments profiling:'+'\n'
resultFile.write(lineResult)
for line in range(time_point_nb-1):
    lineResult = format(line, 'd')+'\t'+ format(int(time_array_min[line][0]), 'd') + '\t' + format(int(time_array_min[line][1]), 'd') + '\t'+ format(len(chronoPolysomeFragProfile[line]), 'd') +'\t' + str(chronoPolysomeFragProfile[line]) + '\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
headerLine = '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'Protein synthesis relative abundance and translational efficiency results:  '+'\n'
resultFile.write(headerLine)
headerLine = '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerline = 'IDtag'+'\t'+'readCount' +'\t' + 'initRate' +'\t' + 'protCount' +'\t' + 'TE' +'\n'
resultFile.write(headerLine)
for i_tr in range(transcriptCardinality):
    lineResult = str(IDtagList[i_tr]) +'\t' + format(readCount[i_tr], 'd') +'\t' + format(initRate[i_tr], '.4f') +'\t'+\
    format(translated_count_list[i_tr], 'd') + '\t' + format(TE[i_tr], '6.2f') + '\n'
    resultFile.write(lineResult)
lineResult = 'simulation run elapsed time (hour) = ' + '\t' +  format(elapsed_t_in_hour, '6.4f') +'\n'
resultFile.write(lineResult)

headerLine = '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'Ribo-seq - Ribosome occupancy maps - RPF density maps results:  '+'\n'
resultFile.write(headerLine)
headerLine = '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
lineResult = 'IDtag' +'\t' + 'length(codons)' + '\t' + 'FC[tr][:]' +'\n'
resultFile.write(lineResult)
for i_tr in range(transcriptCardinality):
    # FC[i_tr][i_codon] per # read copy for this transcript...
    #lineResult = str(IDtagList[i_tr]) +'\t' + format(len(FC[i_tr]), 'd') + '\t' + str(FC[i_tr]) +'\n'
    #resultFile.write(lineResult)

    i_codon = 0
    lineResult = str(IDtagList[i_tr]) +'\t' + format(len(FC[i_tr]), 'd') + '\t' + '['
    sub_lineToAdd = str()
    while i_codon < len(FC[i_tr]):
        if i_codon == len(FC[i_tr])-1:
            sub_lineToAdd += format(FC[i_tr][i_codon], '.2f') +']'+'\n'
            lineResult += sub_lineToAdd
            resultFile.write(lineResult)
            lineResult = str()
            sub_lineToAdd = str()
        else:
            if ((i_codon % 12 != 0) or (i_codon == 0)):
                sub_lineToAdd += format(FC[i_tr][i_codon], '.2f') + '  '
            else:
                sub_lineToAdd += format(FC[i_tr][i_codon], '.2f') +'\n'
                lineResult += sub_lineToAdd
                resultFile.write(lineResult)
                lineResult = str()
                sub_lineToAdd = str()
        i_codon += 1

headerLine = '+++++++++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
headerLine = 'Polysome fragments profiling results:  '+'\n'
resultFile.write(headerLine)
headerLine = '+GLOBAL++++++++++++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
bin_some_nb = len(chronoPolysomeFragProfile[-2])
lineResult = ''
for i_some in range(bin_some_nb):
    lineResult += format(i_some, 'd') + '\t'
lineResult += format(bin_some_nb, 'd') + '\n'
resultFile.write(lineResult)
for i_some in range(bin_some_nb - 1):
    lineResult += format(int(chronoPolysomeFragProfile[-2][i_some]), 'd') + '\t'
lineResult += format(int(chronoPolysomeFragProfile[-2][bin_some_nb - 1]), 'd') + '\n'
resultFile.write(lineResult)
lineResult = str(chronoPolysomeFragProfile[-2]) +'\n'
resultFile.write(lineResult)
# proceed from Here with polysome fragment profiling per transcript ID ... to be continued...
headerLine = '+per TRANSCRIPT +++++++++++++++++++++++'+'\n'
resultFile.write(headerLine)
lineResult = 'i_tr' + '\t' + 'IDtag' + '\t'  + 'polyFragProfileForTr[tr][:]' +'\n'
resultFile.write(lineResult)
for i_tr in range(transcriptCardinality):
    lineResult = format(i_tr, 'd') + '\t' + str(IDtagList[i_tr]) + '\t' + str(polyFragProfileForTr[i_tr]) + '\n'
    resultFile.write(lineResult)
tailerLine = '----------------------------------------------------------------------------' + '\n'
resultFile.write(tailerLine)
resultFile.close()
sys.exit()
