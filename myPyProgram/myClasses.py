# Classes and methods to be used as Object Oriented Programming (OOP)
# to implement the Agent-based model (ABM) of protein synthesis by
# mRNA translation of a pool of ribosomes.
#--------------------------------------------------------------------
#
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
from Bio.Seq import Seq ### comment out
from Bio.Data import CodonTable
from scipy import stats
from suffix_trees import STree
#from myProg02 import time_millis
#from myProg01 import dataRatesDict
#import dataRatesDict
# Define the dictionary of the factors affecting step 1 rate (accommodation and proofreading) upon
# ADAT enzyme silencing. This a dictionary of the factors by witch the rate of step 1 (accomodation)
# will be multiplied because of the inactivation of the A to I editing ADAT enzyme. The NNC codons will
# be slower while the NNU codons will be faster (except for isoleucine for witch AUU is slower)
# upon ADAT silencing.
ADATsilencedDict = {'AUC': 0.57, 'AUU': 0.4, 'AUA': 2.5, 'ACC': 0.85, 'ACU': 1.5, 'ACG': 1.0, 'ACA': 0.67,
'GCC': 0.57, 'GCU': 2.0, 'GCG': 0.71, 'GCA': 0.9, 'CCC': 0.57, 'CCU': 2.0, 'CCG': 0.71, 'CCA': 0.9,
'GUC': 0.85, 'GUU': 1.5, 'GUA': 0.67, 'GUG': 1.05, 'UCC': 0.5, 'UCU': 1.05, 'UCG': 1.1, 'UCA': 0.8,
'AGC': 1.15, 'AGU': 1.05, 'CUC': 0.9, 'CUU': 1.35, 'CUG': 1.4, 'CUA': 0.5, 'UUG': 0.8, 'UUA': 1.2,
'CGC': 0.4, 'CGU': 2.0, 'CGG': 1.0, 'CGA': 1.2, 'AGG': 0.85, 'AGA': 0.95}

# Define the dictionary of the factors affecting rate 1 (accommodation and proofreading) upon
# ELP3-URM1 enzyme silencing. The lack of U34 modification or U34 hypomodification (mcm5s2) or
# lack of s2 (thiolation) slow down rate 1 for codons ending with A in K, E, Q amino acid and affects all
# rate 1 kinetics for 6 codons (3 x 2) for the 3 amino acids K, E, Q (lysine, glutamate, glutamine).
# The G ending codons of these amino acids are read slightly faster in the case of yeast species
# see Nedialkova and Leidel, Cell, (2015) and Ranjan and Rodnina, JACS (2017).
# The factors in the dictionary below are from figure 1 in Nedialkova and Leidel, Cell (2015).
uridine34SilencedDict = {'AAA': 0.33, 'AAG': 1.10, 'CAA': 0.55, 'CAG': 1.0, 'GAA': 0.62, 'GAG': 1.18}
#------------------------------------------------------------------------------
# ribosome footprint length (nt scale):
riboFootprintLength = 33 # 33 nucleotides = 33 pixels = 11 Codons
# -----------------------------------------------------------------------------
# Time clock and Timer set up
# -----------------------------------------------------------------------------
def time_millis(started_t):
    # this function returns the elapsed time in milliseconds since this instance
    # of the main process
    # started when current_time was first recorded at the start of this program
    # instance on this CPU/thread:
    current_time = time.time()
    elapsed_ms = float(format((current_time - started_t) * 1000, '.3f'))
    return elapsed_ms

def decode(codon):
    # This function translates a codon to its cognate amino acid.
    # The function uses the Seq object and its methods fom biopython module.
    # The .translate() function uses the standard genetic code by default
    rna_codon = Seq(str(codon))
    if str(rna_codon.translate()) == '*':
        residue = 'stop'
    else:
        residue = str(rna_codon.translate())
    return residue

def indexOfRead(geneIDlist, readCopyList, readIDstring):
    # This function builds a list of all reads copy.
    # This function returns the index (integer) of the transcripts read list (tr index)
    # corresponding to the geneIDlist given in first argument, the readCopyList
    # given in second argument and the unique ID string given as the third argument
    # (unique ID to query for)
    itemList = []
    for i in range(len(geneIDlist)):
        for j in range(readCopyList[i]):
            itemList.append(geneIDlist[i]+"#"+str(j))
    return itemList.index(readIDstring)

def MaxwellBoltzmannTunnelFactor(upstreamSeq, electrostaticAFDict):
    # This function uses the upstream window sequence of a codon of a given transcript
    # to calculate the Maxwell-Boltzmann factor contributed by the tunnel electrostatic
    # interaction with the mobile window of the last max50 aa upstream of the P-site.
    # The unique input required to compute this factor is the codon.upstreamWindow
    # itself and the model published in Physical Review E and CSBJ by Joiret et al.
    # 1) We convert the upstreamSeq of nucleotides into an ordered sequence of 0, +1 and -1
    # for neutral, + and - amino acids decoded from their triplets in the direction 3' to 5'.
    # Reverse direction, starting from from carboxy-terminal P-site aa to N-terminal aa at
    # the distal end of the tunnel. When you pop the list, you first get the C-terminal
    # end of the nascent chain (ListSign is a LIFO stack list).
    listSign = []
    seqLength = len(upstreamSeq)
    codNumb = int(seqLength/3)
    for i in range(codNumb):
        this_codon = upstreamSeq[i*3: i*3 + 3]
        if ((decode(this_codon) == 'D') or (decode(this_codon) == 'E')):
            listSign.append(-1.0) # negatively charged aa
        else:
            if ((decode(this_codon) == 'K') or (decode(this_codon) == 'R')):
                listSign.append(+1.0) # positively charged aa
            else:
                listSign.append(0.0) # neutral aa
            # end if
        # end if
    # end for
    # 2) The list of axial forces mapping the positions of all aa
    # in window of max 50 codons was retrieved from JSON file to the
    # preloaded dictionary called 'electrostaticAFDict'.
    # If 'K' or 'R' (resp 'D' or 'E') are at position 1 (P-site), the
    # AF is due to PTC electrostatics and is much stronger (=21.2 pN)
    # than the contribution of tunnel: reference: Joiret et al. in CSBJ, 2023.
    # 3) The elementwise dot product of the two lists gives the
    # elementwise list of axial forces acting on the peptidy-tRNA at P-site.
    # 4) The sum of these elementwise forces is the total force. You multiply
    # with axial displacement to get the mechanical work and for a given temperature,
    # you get the MB factor affecting the rate of step 2 (peptide bond formation rate).
    windowLength = len(listSign)
    cumulatedAF = 0.0
    aaSign = 0.0
    afPTC = -21.2
    for i in range(windowLength):
        aaSign = listSign.pop()
        if ((aaSign != 0.0) and (i==0)):
            cumulatedAF += aaSign * afPTC
        cumulatedAF += aaSign * electrostaticAFDict[str(i+1)]
    # end for
    # displacement parallel to axial forces: positive work!
    mechaWork = -1.0 * cumulatedAF * 0.25 # pN.nm
    kBT = 4.282 # pN.nm (at 310.25 K or 37°C)
    maxBoltzFactor = np.exp(mechaWork/kBT)
    return maxBoltzFactor

def ADATsilencingFactor(A_codon):
    # This function returns the factor by witch the rate of step 1 (accommodation) will be multiplied
    # because of the inactivation of the A to I editing ADAT enzyme. Upon ADAT silencing, the NNC codons will
    # be slower while the NNU codons will be faster (except for isoleucine for witch AUU is slower).
    if decode(A_codon) in ['T', 'A', 'P', 'S', 'L', 'I', 'V', 'R']:
        return ADATsilencedDict[str(A_codon)]

def U34hypoModificationFactor(A_codon):
    # This function returns the factor by witch the rate of step 1 (accommodation) will be multiplied
    # because of the inactivation of the U34 modiying enzyme (ELP3, URM1, ...). Upon U34 modifcation
    # silencing, the NNA codons will
    # be slower while the NNG codons will be faster.
    # Note that A_codon and Acodon would both be the same here because there are no U at all in
    # these 6 targeted codons. So there was no need to back_transcribe the A_codon (rna_codon) into a dna_codon.
    if decode(A_codon) in ['K', 'E', 'Q']:
        return uridine34SilencedDict[str(A_codon)]

def pdf_HYPOEXP(x, lamb1, lamb2, lamb3):
    # This function returns the probability density function of an hypo-exponential distribution
    # with three exponential rates parameters.
    # lambda_list is the list of the n different parameters of the n exponential distribution
    # for which we want the distribution of their sum or the convolution product of the n
    # exponential distributions
    n = 3
    lambda_list = [lamb1, lamb2, lamb3]
    productsOfLambdas = 1.0
    frac_differences = []
    for param_idx in range(n):
        productsOfLambdas *= lambda_list[param_idx]
        lambda_ref = lambda_list[param_idx]
        lambda_denum = []
        productOfDifference = 1.0
        for i in range(n):
            if lambda_list[i] != lambda_ref:
                lambda_denum.append(lambda_list[i])
                productOfDifference *= (lambda_list[i] - lambda_ref)
        frac_differences.append(productOfDifference)
    sum = 0
    for param_idx in range(n):
        sum += np.exp(-1.0*lambda_list[param_idx]*x)/frac_differences[param_idx]
    sum *= productsOfLambdas
    return sum

# Function to sample a value drawn from an exponential distribution specified by its rate parameter:
def expRandomTimeSample(rate):
    # rate (lambda) is the inverse of scale parameter (beta).
    # beta, the scale parameter, is supposed to be in milliseconds units. So, the rate is [1/ms] units.
    # The function returns a random time (in milliseconds) sampled from the exponential pdf
    #return np.random.Generator.exponential(scale=1.0/rate, size=None)
    return np.random.exponential(scale=1.0/rate, size=None)

# Function to sample a value from a hypo-exponential distribution
# specified with its 3 rates (1/time) parameters:
def hypoExpRandomSample(rate1, rate2, rate3):
    # The function returns a random time (in ms units) sampled fom the hypoexponential pdf.
    # The hypoexponential distribution is the convolution product of 3 exponetial pdfs.
    # Equivalently, the hypoexponential random variable is the sum of three iid hypoexponential pdfs:
    expTime1 = expRandomTimeSample(rate1)
    expTime2 = expRandomTimeSample(rate2)
    expTime3 = expRandomTimeSample(rate3)
    hypoExpTime = expTime1 + expTime2 + expTime3
    return hypoExpTime

# Function assigning the queueing time of a ribosome on a codon at the P-site
# depending on codon at P and A-sites, upstream sequence context, proline status factor,
# ADAT silencing factor and possibly other factors (downstream secondary structures, ...):
def hypoExpRandomTime(P_codon, A_codon, upSeq, ratesDict, electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus):
    # The function returns the queueing time a ribosome has to wait when on the P-site codon.
    # This function needs the data of the dictionary 'DataRatesDict' that was downloaded
    # from dataFile1 = dataJSONyeast.json. This dictionary is built in the main program.
    # The dictionary is passed here in this function argument.
    # The function is imported from this class module in the main
    # program. P_codon and A_codon are in mRNA nucleotides (U instead of T) but the rates dictionary
    # has cds (dna) nucleotides and so U has to be replaced with T with biopython.back_transcribe()
    # to properly access the keys of the dictionary:
    Pcodon = P_codon.back_transcribe()
    Acodon = A_codon.back_transcribe()
    if ((Acodon == 'TAA') or (Acodon == 'TAG') or (Acodon == 'TGA')):
        Acodon = Pcodon
    rate1 = ratesDict[Acodon][0]
    rate2 = ratesDict[Pcodon][1]
    rate3 = ratesDict[Pcodon][2]
    # Depending on context, these rates can be modulated by different factors:
    # - Proline factor and EF-P or eIF5A elongation factor depletion affects rate2:
    # ...if species == 'yeast', the rate is twofold slower than for E.coli.
    # ... and the peptide bond formation rate for E.coli corresponds at least to t_1/2 = 5102 ms
    if ((prolineStatus) and (decode(P_codon) == 'P')):
        rate2 = 0.5 * 0.6931 /5102 # ln2 / t_1/2
    if ((EFPdepletedStatus) and (decode(A_codon) == 'P')):
        # there is a fold change of 90 (90 times slower) in the elongation rate
        # when proline is at the A site and EF-P or eIF5A is depleted:
        rate2 = rate2/90.0
    # - Exit tunnel electrostatics interaction mainly affects rate2:
    if (tunnelStatus):
        rate2 = rate2 * MaxwellBoltzmannTunnelFactor(upSeq, electrostaticAFDict)
    # - ADAT2 silencing affects rate1 for 8 amino acids and 37 ADAT2 sensitive codons:
    if ((ADAT2silencedStatus) and (decode(A_codon) in ['T', 'A', 'P', 'S', 'L', 'I', 'V', 'R'])):
        rate1 = rate1 * ADATsilencingFactor(A_codon)
    # - U34 hypomodification affects rate 1 for 3 amino acids and 6 ELP3-URM1 sensitive codons:
    if ((U34hypoModificationStatus) and (decode(A_codon) in ['K', 'E', 'Q'])):
        rate1 = rate1 * U34hypoModificationFactor(A_codon)
    # ...
    # Finally, the queueing time is sampled from the hypoexponetial distribution:
    return hypoExpRandomSample(rate1, rate2, rate3)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CLASSES for the agents used in the Agent-based model of protein synthesis by ribosome:
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# -----------------------------------------
# Transcript class, attributes and methods:
# -----------------------------------------
#
# The Transcript class is called the 'Lattice' class. The name is borrowed from the TASEP model.
#
class Lattice:
    # constructor:
    # an instance of a lattice, aka transcript, keeps the following (hidden) data:
    # constructor attributes assigned at first instantiation
    def __init__(self, transcLength, transcSequence, transcCopy, transcName, transcNumID):
        self.__length = transcLength  # in codons scale
        self.__sequence = transcSequence
        self.__mRNAlength = len(transcSequence)
        self.__copy = transcCopy
        self.__name = transcName
        self.__numID = transcNumID
        self.__uniqueID = str(transcName)+'#'+str(transcNumID)
        #self.__ProStatus = transcProStatus # maybe unnecessary
        #self.__TunnelStatus = transcTunnelStatus # maybe unnecessary
    # attributes (other housekeeping attributes that can be affected within an instantiation
    # and updated according to interactions with other objects such as ribosomes object and their methods)
        self.__currentFootPrint = [-1] # array of ribosome positions when at P-sites (nucleotide scale)
        self.__translatingTime = [] # list of all translating times that each translating ribosome
        #has spent on this particular transcript.
        self.__readTranslated = 0 # only this (read) copy of the transcript is iterated upon translation
        self.__countTranslated = 0 # CAUTION: aggregated over all the copies of this transcript name at a
        # given time point.
        self.__rpfCount = 0 # number of RPF (ribosome protected fragment) on the transcript.
    # methods: getters (accessors):
    def get_length(self):
        return int(self.__length)
    def get_sequence(self):
        return self.__sequence
    def get_mRNAlength(self):
        return int(self.__mRNAlength)
    def get_copy(self):
        return int(self.__copy)
    def get_name(self):
        return self.__name
    def get_numID(self):
        return int(self.__numID)
    def get_uniqueID(self):
        return self.__uniqueID
    def get_rpfCount(self):
        return self.__rpfCount
    def get_currentFootprint(self):
        return self.__currentFootPrint
    def get_readTranslated(self):
        return self.__readTranslated
    def get_countTranslated(self):
        return self.__countTranslated
    # methods setters (mutators):
    def set_currentFootprint_initiation(self, nt_pos):
        if ((nt_pos == 0) and (self.__currentFootPrint[0] == -1)):
            # this is the first initiation for this transcript:
            self.__currentFootPrint.pop()
            self.__currentFootPrint.append(nt_pos) # push method does not exist in python lists
        #else:
        if ((nt_pos == 0) and (self.__currentFootPrint[0] > -1)):
            # this is a polysome initiation attempt:
            freeForInitiation = True
            for i_pos in range(len(self.__currentFootPrint)):
                if self.__currentFootPrint[i_pos] <= nt_pos + riboFootprintLength:
                    freeForInitiation = False
            if freeForInitiation:
                self.__currentFootPrint.append(nt_pos)

    def set_currentFootprint_translocation(self, nt_pos):
        if nt_pos == (self.__mRNAlength - 3):
            # this is a termination and the footprint list should shrink by one element:
            # slicing assignment or assigning empty list to a slice or removing the element:
            #indexToRemove = self.__currentFootPrint.index(nt_pos)
            # or ... del self__currentFootPrint[indexToRemove]
            # self.__currentFootPrint.remove(nt_pos)
            #print('len of self.__currentFootPrint=', len(self.__currentFootPrint))
            if (len(self.__currentFootPrint)>1) and (self.__currentFootPrint[0]==nt_pos):
                del self.__currentFootPrint[0]
            if (len(self.__currentFootPrint)==1) and (self.__currentFootPrint[0]==nt_pos):
                self.__currentFootPrint = [-1]

        else:
            # this is a standard translocation requiring to update the positions list:
            for i_pos in range(len(self.__currentFootPrint)):
                if self.__currentFootPrint[i_pos] == nt_pos:
                    self.__currentFootPrint[i_pos] = nt_pos + 3
                    break
                # end if
            # end for
        # end else

    def set_readTranslated(self):
        # This mutator increases the number of protein counts that were produced from
        # translating this read up to and including the stop codon:
        self.__readTranslated += 1
        return self.__readTranslated
    def set_countTranslated(self, addCount):
        # This mutator increases the number of protein counts that were produced from
        # translating this transcript ID (aggregating all its read copies):
        self.__countTranslated += addCount
        return self.__countTranslated
    def set_rpfCount(self, riboCountOnThisRead):
        self.__rpfCount = riboCountOnThisRead
        return self.__rpfCount

# end Lattice class
# -----------------------------------------------------------------------------------------------------------------------------
#
#  The codon class defines the codon object type. A codon instantiation has a type (nucleotides triplet), belongs to a given transcript
#  and has a precise position on this transcript (nucleotide scale). A codon can be a P-site or a A-site codon.
#  Each codon is assigned a queuing time for a ribosome to proceed an elongation cycle.
#
class Codon:
    # constructor:
    # an instance of a codon keeps the following (hidden) data:
    def __init__(self, tempTranscript, tempPosition, dataRatesDict):
        self.type = tempTranscript.get_sequence()[tempPosition:tempPosition + 3] # slicing in python
        self.belongs = tempTranscript
        self.position = tempPosition # index position of codon in the transcript it comes from (in nucleotide scale
        if tempPosition <= tempTranscript.get_mRNAlength() - 6:
            self.typeAsite = tempTranscript.get_sequence()[tempPosition + 3:tempPosition + 6] # A-site codon
        if tempPosition <= 147:
            self.upstreamWindow = tempTranscript.get_sequence()[0:tempPosition]
        if tempPosition >= 148:
            self.upstreamWindow = tempTranscript.get_sequence()[tempPosition - 49*3:tempPosition + 3]
        self.ratesDict = dataRatesDict
    # end constructor

    # methods:
    def queueTimeAtPsite(self, electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus): # sp = ''yeast' species("yeast" or "coli" or ...)
        if self.position == self.belongs.get_mRNAlength() - 3:
            # assign termination rate:
            terminationRate = 0.5 # mean time for termination is 1/0.5 millisecond (2 ms)
            return expRandomTimeSample(terminationRate)
        else:
            # 'still elongating...', assign elongation rate (queueing time on the codon) :
            return hypoExpRandomTime(self.type, self.typeAsite, self.upstreamWindow, self.ratesDict, electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus)
    # end methods
# end Codon class
# -----------------------------------------------------------------------------------------------------------------------------
#
#  The Ribosome class defines the ribosome object type. A ribosome instantiation has a type 'free' at
#  instantiation (but can change to 'initiated', 'translating'), a traffic jam
#  status (ribosome stalling), time events records: -clockTimer, -startedInitiation, -startedTranslocation
#  - endTermination, -
#
class Ribosome:
    # constructor:
    def __init__(self, started_time):
        self.type ='free'
        self.clockTimer = time_millis(started_time) # time()
        #self.clockTimer = time.time()*1000
        self.startedInitiation = 0
        self.timeSinceInitiation = 0
        self.transcriptName = ''
        self.uniqueReadID = ''
        self.goStatus = True
        self.indexCodon = -3 # indexCodon is the ribosome's position (nucleotide scale) on the codon
                             # in the mRNA sequence (multiple of 3)
        self.startedTranslocation = 0
        self.translocqueueingTimer = 0
        self.dwellTimePrevious = 0
        self.endedTermination = 0
        self.__initiationRate = 1.0/60000 # t1/2 = 41.6 seconds. Half of the free ribosomes will be engaged after 41 sec.
        self.InitiationTimerSP = expRandomTimeSample(self.__initiationRate) # sample a time from an exponential distribution
        self.TranslocationTimerSP = 0
        # This will be the list of clocktime events records upon each translocation:
        self.cumulatedTimeAtThisCodon = []
        self.cumulatedTranslocationSP = 0
    # methods:
    # methods: getters (accessors):
    def get_indexCodon(self):
        return int(self.indexCodon)
    # mutators:
    def set_initiationRate(self, init_Rate):
        self.__initiationRate = init_Rate
    # method - a ribosome can take an initiation action
    # initiation method:
    def initiates(self, objTranscript, started_time):
        # A free ribosome may initiate on a given unique transcript read.
        # A free ribosome receives in argument a randomly chosen unique transcript read,
        # monitors the timer and waits until the initiation time setpoint is reached.
        # When elapsed time larger than setpoint:
        # checks if the first 10 codons (from 0 to nine included) of the given transcript read
        # are free of any previous P-site occupying ribosomes.
        # If yes:
        # 1) This 'free' ribosome becomes an 'initiated' ribosome and this free ribosome disappears:
        # it has switched the population category.
        # 2) The new position is the one of the start codon (the P site is at 0 of the transcript read).
        # If no:
        # 1) The ribosome stays 'free'
        # 2) Its timer is reset to zero.
        gotFootprint = objTranscript.get_currentFootprint()
        footprinted = False
        listLength = len(gotFootprint)
        for psite in range(listLength):
            if ((gotFootprint[psite] >= 0) and (gotFootprint[psite] < riboFootprintLength)):
                footprinted = True
                break # you can get out of the for loop if a ribosome footprint was detected
        # end for LOOP
        if (((time_millis(started_time) - self.clockTimer) >= self.InitiationTimerSP) \
        and (self.type == 'free') and (not footprinted)):
            # turn the ribosome into an initiated ribosome:
            self.type = 'initiated'
            # record time event when initiation started for this ribosome:
            self.startedInitiation = time_millis(started_time)
            self.startedTranslocation = time_millis(started_time)
            # keep track of the name of the current transcript on which this ribosome is on:
            self.transcriptName = objTranscript.get_name()
            self.uniqueReadID = objTranscript.get_uniqueID()
            self.indexCodon = 0 # upon initiation, the ribosome P site is on codon 0 = 'AUG'
            # update transcript data wrt ribosomal occupancy:
            objTranscript.set_currentFootprint_initiation(0) # indexCodon=0 upon initiation
            self.cumulatedTimeAtThisCodon.append(time_millis(started_time))
            # redraw a new random initiation time setpoint from the exponential distribution:
            self.InitiationTimerSP = expRandomTimeSample(self.__initiationRate)
    # end initiate method

    # method - a ribosome can take an action to check occupancy for other ribosomes on a given
    # transcript read.
    # checkOccupancy method:
    def checkOccupancy(self, forThisTranscript):
        # This method checks if the next 6 codons from the current P-site position of this ribosome
        # are free for the next translocation. It returns a boolean (called freeToGo) True if the ribosome
        # is free to translocate of false if not. If freeToGo is False: a traffic jam (ribosome congestion)
        # situation arise.
        freeToGo = True
        # check for the right transcript read:
        if self.uniqueReadID == forThisTranscript.get_uniqueID():
            # compare current index position of the ribosome with all occupation sites
            # for this current transcript read:
            for occupiedByOtherRibo in range(len(forThisTranscript.get_currentFootprint())):
                if ((self.indexCodon < forThisTranscript.get_currentFootprint()[occupiedByOtherRibo]) \
                and (forThisTranscript.get_currentFootprint()[occupiedByOtherRibo] - self.indexCodon <= riboFootprintLength)):
                    freeToGo = False
                    break
                # end if
            # end for
        # end if
        return freeToGo
    # end checkOccupancy method

    # method - a ribosome can take an action to translocate on a given
    # transcript read.
    # translocates method:
    def translocates(self, objTranscript, dataRatesDict, electrostaticAFDict, started_time, IDtagList, readCount, RibosSeqSPtruth, RibosSeqRRTtruth, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus):
        # This ribosome method receives in argument a transcript read object and the rates dictionary..
        # Monitors the timer and wait until the setpoint for translocation is reached.
        # When elapsed time larger than setpoint:
        # Checks if the ribosome is free to go by using the .checkOccupancy method.
        # if yes:
        # 1) This initiated ribosome or translating ribosome can translocate.
        # 2) if the ribosome is currently on the last 'stop' codon, it will terminate in
        # an elapsed time span determined by the termination rate and turn into a 'free'
        # ribosome type.
        # Inquire about the codon type that is currently in the A-site. The A-site is one codon
        # downstream the P-site which is known from the current position of the ribosome. We
        # kept track of the name and uniqueID of the current transcript on which this ribosome is
        # on and the translocates method only applies to this transcript read.
        # Check the name and ID of transcript being translated:
        right_transcript = False
        if self.uniqueReadID == objTranscript.get_uniqueID(): # are you on the right transcript?
            right_transcript = True

        # Check occupancy status downstream this codon:
        self.goStatus = self.checkOccupancy(objTranscript)
        PsiteCodon = [] # this is defined as a list but will be used as LIFO stach (pop, push=append)
        queueingTime = time_millis(started_time) - self.startedTranslocation
        # The ribosome is at least initiated, meaning that the indexCodon is >=0 and <= stop - 3 nts:
        if ((right_transcript==True) and (self.indexCodon >= 0) and (self.indexCodon <= objTranscript.get_mRNAlength() - 6)):
            # identify and instantiate the codon object at the P-site:
            PsiteCodon.append(Codon(objTranscript, self.indexCodon, dataRatesDict))
            # draw a timer setpoint for the queueing timeout (queues to translocate):
            self.TranslocationTimerSP = PsiteCodon.pop().queueTimeAtPsite(electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus)
            self.cumulatedTranslocationSP += self.TranslocationTimerSP
            #queueingTime = time_millis() - self.startedTranslocation
            if ((queueingTime >= self.TranslocationTimerSP) \
            and (self.type != 'free') \
            and (self.goStatus)):
                self.type = 'translating'
                # translocates:
                # updates transcript read data wrt ribosomal occupancy
                objTranscript.set_currentFootprint_translocation(self.indexCodon)
                # you have just translocated to the next codon, so you can record the dwell time
                # on the codon you just left and save the cumulated SP (in case traffic jam occured):
                index_tr = indexOfRead(IDtagList, readCount, objTranscript.get_uniqueID())
                RibosSeqSPtruth[index_tr][int(self.indexCodon/3)] += self.cumulatedTranslocationSP
                # calculate dwell time on the codon you just left:
                self.cumulatedTimeAtThisCodon.append(time_millis(started_time))
                self.dwellTimePrevious = self.cumulatedTimeAtThisCodon[len(self.cumulatedTimeAtThisCodon)-1] \
                                        - self.cumulatedTimeAtThisCodon[len(self.cumulatedTimeAtThisCodon)-2]
                if self.indexCodon >= 3:
                    RibosSeqRRTtruth[index_tr][int(self.indexCodon/3)-1] += self.dwellTimePrevious
                # update ribosome position:
                self.indexCodon += 3
                # reset cumulated time on this codon (in case traffic jam occurred):
                self.cumulatedTranslocationSP = 0
                # reset timer:
                self.startedTranslocation = time_millis(started_time)
            if ((not self.goStatus) and (self.type != 'free')): # there is another ribosome downstream that jams the translocation:
                PsiteCodon.append(Codon(objTranscript, self.indexCodon, dataRatesDict))
                # re-draw a timer setpoint for the queueing timeout (queues to translocate):
                self.TranslocationTimerSP = PsiteCodon.pop().queueTimeAtPsite(electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus)
                self.cumulatedTranslocationSP += self.TranslocationTimerSP
                self.startedTranslocation = time_millis(started_time)
            # end if
        # end if
        if ((right_transcript==True) and (self.indexCodon == objTranscript.get_mRNAlength() - 3)): # The ribosome is at a stop codon.
        # The next translocation is actually here a termination step:
            # queues to terminate:
            PsiteCodon.append(Codon(objTranscript, self.indexCodon, dataRatesDict))
            # draw a timer setpoint for the queueing timeout (queues to terminate):
            self.TranslocationTimerSP = PsiteCodon.pop().queueTimeAtPsite(electrostaticAFDict, prolineStatus, EFPdepletedStatus, tunnelStatus, ADAT2silencedStatus, U34hypoModificationStatus)
            self.cumulatedTranslocationSP += self.TranslocationTimerSP
            if ((queueingTime >= self.TranslocationTimerSP) \
            and (self.type != 'free')):
                # terminates:
                #print('index Codon=', self.indexCodon, 'mrna length=', objTranscript.get_mRNAlength())
                #print('current footprint of this ribosome=')
                #print('current footprint =', objTranscript.get_currentFootprint())
                # updates transcript object data wrt ribosomal occupancy:
                objTranscript.set_currentFootprint_translocation(self.indexCodon)
                #print('updated footprint=', objTranscript.get_currentFootprint())
                # calculate dwell time on the codon you just left:
                self.cumulatedTimeAtThisCodon.append(time_millis(started_time))
                self.dwellTimePrevious = self.cumulatedTimeAtThisCodon[len(self.cumulatedTimeAtThisCodon)-1] \
                                        - self.cumulatedTimeAtThisCodon[len(self.cumulatedTimeAtThisCodon)-2]
                self.startedTranslocation = time_millis(started_time)
                self.type ='free'
                self.clockTimer = time_millis(started_time)
                self.startedInitiation = time_millis(started_time)
                self.timeSinceInitiation = 0
                self.transcriptName = ''
                self.uniqueReadID = ''
                self.goStatus = True
                self.indexCodon = -3 # indexCodon is the ribosome's position (nucleotide scale) on the codon
                                     # in the mRNA sequence (multiple of 3)
                self.startedTranslocation = 0
                self.translocqueueingTimer = 0
                self.dwellTimePrevious = 0
                self.endedTermination = 0
                #self.__initiationRate = 1.0/60000 # t1/2 = 41.6 seconds. Half of the free ribosomes will be engaged after 41 sec.
                self.InitiationTimerSP = expRandomTimeSample(self.__initiationRate) # sample a time from an exponential distribution
                self.TranslocationTimerSP = 0
                # This will be the list of clocktime events records upon each translocation:
                self.cumulatedTimeAtThisCodon = []
                self.cumulatedTranslocationSP = 0
                # a protein molecule was just produced: update the count of protein
                # for this transcript unique read:
                objTranscript.set_readTranslated()
                # and update the count ofprotein for this transcript name: No don't, not here
                # objTranscript.set_countTranslated()
            # end if (for stop codon)
        # end termination
    # end translocates method
# end class Ribosome
