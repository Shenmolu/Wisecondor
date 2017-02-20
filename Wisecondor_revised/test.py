##############################################################################
#                                                                            #
#    Test a maternal DNA sample for fetal Copy Number Aberrations.           #
#    Copyright(C) 2013  TU Delft & VU University Medical Center Amsterdam    #
#    Author: Roy Straver, r.straver@vumc.nl                                  #
#                                                                            #
#    This file is part of WISECONDOR.                                        #
#                                                                            #
#    WISECONDOR is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    WISECONDOR is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with WISECONDOR.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                            #
##############################################################################



import sys
import pickle
import numpy
import math
import argparse
import warnings

numpy.seterr('ignore')

lookUpTable = None

def getZScore(freq, reference):
    average = numpy.average(reference)
    stddev  = numpy.std(reference)
    if stddev == 0:
        return 0
    Z = (freq - average) / stddev
    return Z

def getReference(readFreq,chrom,iBin,markedBins,minBins,maxBins,maxDist):
    '''Get locations of all qualified reference bins for the given target bin'''
    reference = []
    if len(lookUpTable[chrom]) <= iBin:
        return reference
    for value in lookUpTable[chrom][iBin]:
        if (value[0],value[1],) in [marked[:2] for marked in markedBins]:
            continue
        if len(readFreq[value[0]]) > value[1]:
            # Only add bin if the distance is small enough
            if value[2] <= maxDist:
                if readFreq[value[0]][value[1]] != 'NA':
                    reference.append(readFreq[value[0]][value[1]])
                else:
                    continue
            else:
                break # Stop trying, only worse bins to come

        if len(reference) >= maxBins:
            break

    # Ignore bin because of too few reference bins
    if len(reference) < minBins:
        return []
    # Ignore bins after maxBins is reached
    return reference

def markBins(readFreq,maxRounds,minBins,maxBins,maxDist,smoothRange):
    totalBins = sum([len(readFreq[chrom]) for chrom in chromList])
    prevMarks = [('',0,0)]
    markedBins = []
    rounds = 1
    zScoresDict = dict()
    zSmoothDict = dict()
    
    while ([marked[:2] for marked in prevMarks] != [marked[:2] for marked in markedBins]) and rounds <= maxRounds:
        print '\tRound: ' + str(rounds) + '\tMarks: ' + str(len(markedBins))
        rounds += 1
        prevMarks = markedBins
        markedBins = []

        for chrom in chromList:
            zScores = []

            for tBin in range(len(readFreq[chrom])):
                freq = readFreq[chrom][tBin]
                if freq == 'NA':
                    zScores.append('NA')
                    continue
                reference = getReference(readFreq,chrom,tBin,prevMarks,minBins,maxBins,maxDist)
                zValue = getZScore(freq, reference)

                if (abs(zValue) >= 3):
                    markedBins.append((chrom,tBin,zValue))

                zScores.append(zValue)

            zScoresDict[chrom] = zScores

    print 'Stopped\tMarks: ' + str(len(markedBins))

    for chrom in zScoresDict:
        zSmooth = [1] * len(zScoresDict[chrom])
        for tBin in range(len(zScoresDict[chrom])):
            temp = zScoresDict[chrom][max(0,tBin-smoothRange):min(tBin+smoothRange+1,len(zScoresDict[chrom]))]
            # Remove NA values from the set
            temp = [val for val in temp if not (val == 'NA')]
            # Get lost, frigging peak.
            temp.sort()
            temp = temp[1:-1]
            if len(temp) > 0:
                zSmooth[tBin] = numpy.sum(temp)/numpy.sqrt(len(temp))
            else:
                zSmooth[tBin] = 'NA'
        zSmoothDict[chrom] = zSmooth

    return zScoresDict, zSmoothDict

# --- MAIN ---
import argparse
parser = argparse.ArgumentParser(description='Calculate z-scores',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('sample', type=str,
                   help='sample to be tested (.correct)')
parser.add_argument('reference', type=str,
                   help='reference table used for within sample comparison (pickle)')
parser.add_argument('outfile', type=str,
                   help='output results to this file')

parser.add_argument('-female', action='store_true', default=False,
                    help='turn on if gender is female')

parser.add_argument('-maxrounds', default=5, type=int,
                   help='maximum amount of rounds used to calculate z-score')
parser.add_argument('-refminbin', default=10, type=int,
                   help='minimum number of reference bins, ignore target bin if there are less reference bins available')
parser.add_argument('-refmaxbin', default=100, type=int,
                   help='maximum number of reference bins, ignore any reference bin after')

parser.add_argument('-window', default=5, type=int,
                   help='window size for sliding window approach, number of bins is considered in each direction (i.e. using 3 results in using 3+1+3=7 bins per call)')

args = parser.parse_args()


print '# Script information:'
print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
    print '\t'.join([arg,str(argsDict[arg])])

print '\n# Processing:'
print 'Loading:\tSample:\t' + args.sample

chromList = [str(chrom) for chrom in range(1,23)]
chromList.append('X')
chromList.append('Y')

readFreq = dict()
for chrom in chromList:
    readFreq[chrom] = []
try:
    with open(args.sample,'r') as sampleFile:
        next(sampleFile)
        for line in sampleFile:
            words = line.split()
            if words[8] != 'NA':
                if (words[0][3:] == 'X' or words[0][3:] == 'Y') and not args.female:
                    readFreq[words[0][3:]].append(float(words[8])*2)
                else:
                    readFreq[words[0][3:]].append(float(words[8]))
            else:
                readFreq[words[0][3:]].append('NA')
except IOError as err:
    print 'IOError:' + str(err)
    sys.exit()

# Load reference table
print 'Loading:\tReference Table\t' + args.reference
try:
    with open(args.reference,'rb') as refFile:
        ref = pickle.load(refFile)
        lookUpTable = ref['lookUp']
        maxDist = ref['maxDist']
except pickle.PickleError as perr:
    print 'Pickle error:\t' + str(perr)
    sys.exit()


print ''
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    zScoresDict,zSmoothDict = markBins(readFreq,args.maxrounds,args.refminbin,args.refmaxbin,maxDist,args.window)

zScores = []
zSmooth = []
for chrom in chromList:
    zScores.extend(zScoresDict[chrom])
    zSmooth.extend(zSmoothDict[chrom])

try:
    with open(args.outfile,'wb') as outfile:
        with open(args.sample,'rU') as infile:
            header = next(infile)
            header = header.strip()
            outfile.write(header + '\tZ_score\tsw_Z_score')
            i = 0
            for line in infile:
                line = line.strip()
                outfile.write('\n'+line+'\t'+str(zScores[i])+'\t'+str(zSmooth[i]))
                i += 1
except IOError as err:
    print 'IOError:' + str(err)
    sys.exit()

print '\n# Finished'
