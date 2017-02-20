##############################################################################
#                                                                            #
#    Find optimal cutoff for 'good' bins for a certain reference set.        #
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



import argparse
import glob
import sys
import pickle
import numpy

parser = argparse.ArgumentParser(description='Create a new reference table from a set of reference samples, outputs table as pickle to a specified output file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('refin', type=str,
                    help='directory containing samples (.correct) to be used as reference')
parser.add_argument('refout', type=str,
                    help='reference table output, used for sample testing (pickle)')
parser.add_argument('-female', action='store_true', default=False,
                    help='turn on if gender is female')
parser.add_argument('-ignore', type=int, default=0,
                    help='ignore x highest scoring sample samples per bin distance calculation, 0 to keep all, use if you are unsure about what samples have aberrations and assume at most x have the same aberration, beware: results become less accurate')
parser.add_argument('-maxbin1', default=250, type=int,
                    help='maximum number of reference bins for each target bin before removing neighboring bins')
parser.add_argument('-maxbin2', default=100, type=int,
                    help='maximum number of reference bins for each target bin after removing neighboring bins (final)')
parser.add_argument('-refmaxval', default=1000000, type=int,
                    help='start cutoff value for determining good quality reference bins')
parser.add_argument('-refmaxrep', default=3, type=int,
                    help='amount of improval rounds for determining good quality reference bins')
args = parser.parse_args()

print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
    print '\t'.join([arg,str(argsDict[arg])])

print '\n# Processing:'

def getReferenceBins(samples,tChrom):
    '''Get qualified reference bins to bins on iChrom'''
    chromosomeDistances = []

    # For each target bin on tChrom
    tLen = min([len(samples[key][tChrom]) for key in samples.keys()])
    for tBin in range(0,tLen):
        binDistances = []
        for rChrom in chromList:
            if rChrom != tChrom:
                # For each reference bin on rChrom
                rLen = min([len(samples[key][rChrom]) for key in samples.keys()])
                for rBin in range(0,rLen):
                    # Take the difference over all samples...
                    distances = []
                    for sample in samples:
                        tVal = samples[sample][tChrom][tBin]
                        rVal = samples[sample][rChrom][rBin]
                        # Don't try to match with 0 bins later on
                        if (tVal == 0) or (rVal == 0):
                            distances = []
                            break
                        distances.append(pow((tVal - rVal),2))

                    # Append new found distance
                    if len(distances) > 0:
                        if args.ignore > 0:
                            distances.sort()
                            binDistances.append((rChrom, rBin, sum(distances[:-args.ignore])))
                        else:
                            binDistances.append((rChrom, rBin, sum(distances)))

        # Get bins with smallest distances
        binDistances = sorted(binDistances, key=lambda x : x[2])
        binDistances = binDistances[:args.maxbin1]

        # Don't take bins close to eachother, take the best one instead
        avail = dict()
        for rChrom in chromList:
            avail[rChrom] = [1] * (min([len(samples[key][rChrom]) for key in samples.keys()])+1)
        tempList = []
        for rBin in binDistances:
            if avail[rBin[0]][rBin[1]] == 1:
                tempList.append(rBin)
                avail[rBin[0]][max(0,rBin[1]-1)] = 0
                avail[rBin[0]][rBin[1]+1] = 0

        chromosomeDistances.append(tempList[:args.maxbin2])
    return chromosomeDistances


def getOptimalCutoff(refTable, repeats, optimalCutoff):
    '''Return the value of optimal cutoff'''
    for i in range(0,repeats):
        # Select the best matching reference bins
        bestMatch = []
        for chrom in refTable:
            for tbin in refTable[chrom]:
                if len(tbin) > 0:
                    if float(tbin[0][2]) < optimalCutoff:
                        bestMatch.append(float(tbin[0][2]))

        average = numpy.average(bestMatch)
        stddev  = numpy.std(bestMatch)
        optimalCutoff = average + 3 * stddev
    return optimalCutoff


# Load reference samples
print 'Loading reference samples'
referenceFiles = glob.glob(args.refin + '/*.correct')

chromList = [str(chrom) for chrom in range(1,23)]
chromList.append('X')
chromList.append('Y')

samples = dict()
for refFile in referenceFiles:
    print '\tLoading:\t' + refFile
    readFreq = dict()
    for chrom in chromList:
        readFreq[chrom] = []

    try:
        with open(refFile, 'r') as infile:
            next(infile)
            for line in infile:
                words = line.split()
                if words[8] != 'NA':
                    if (words[0][3:] == 'X' or words[0][3:] == 'Y') and not args.female:
                        readFreq[words[0][3:]].append(float(words[8])*2)
                    else:
                        readFreq[words[0][3:]].append(float(words[8]))
                else:
                    readFreq[words[0][3:]].append(0)
            samples[refFile] = readFreq
    except IOError as err:
        print 'Fail to read:\t' + refFile

# Build reference table
print 'Building reference table'
refTable = dict()
for chrom in chromList:
    refTable[chrom] = []

for tChrom in chromList:
    print '\tTargeting chromosome:\t' + tChrom
    refTable[tChrom] = getReferenceBins(samples,tChrom)


# Remove bins based on optimal cutoff
print '\nDetermining reference cutoffs'
maxDist = getOptimalCutoff(refTable,args.refmaxval,args.refmaxrep)

print '\tRemoving outliers'
lookUp = dict()
for chrom in chromList:
    lookUp[chrom] = []

for chrom in chromList:
    for tBin in range(0,len(refTable[chrom])):
        lookUp[chrom].append([])
        for rBin in refTable[chrom][tBin]:
            if rBin[2] < maxDist:
                lookUp[chrom][tBin].append(rBin)

# Write reference table and reference cutoff to file
output = dict()
output['lookUp'] = lookUp
output['maxDist'] = maxDist
try:
    with open(args.refout, 'wb') as outfile:
        print 'Writing reference to file'
        pickle.dump(output, outfile)
except pickle.PickleError as perr:
    print('Pickle error:' + str(perr))
    sys.exit()

print '\n# Finished'