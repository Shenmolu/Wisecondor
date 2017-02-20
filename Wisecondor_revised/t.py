import sys
import argparse
import glob


parser = argparse.ArgumentParser(description='Create a new reference table from a set of reference samples, outputs table as pickle to a specified output file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
'''
parser.add_argument('refdir', type=str,
                    help='directory containing samples to be used as reference')
parser.add_argument('refout', type=str,
                    help='reference table output, used for sample testing (pickle)')

parser.add_argument('-ignore', type=int, default=0,
                    help='ignore x highest scoring control samples per bin distance calculation, 0 to keep all, use if you are unsure about what samples have aberrations and assume at most x have the same aberration, beware: results become less accurate')

args = parser.parse_args()

print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
    print '\t'.join([arg,str(argsDict[arg])])

print 'Loading reference samples'
controls = dict()
referenceFiles = glob.glob(args.refdir + '/*.correct')  # was .gcc
chromList = [str(chrom) for chrom in range(1,23)]       # Add X and Y chromosomes!!!!!!!
chromList.append('X')
chromList.append('Y')
for refFile in referenceFiles:
    print '\tLoading:\t' + refFile
    curFile = dict()
    for chrom in chromList:
        curFile[chrom] = []
    f = open(refFile, 'r')
    next(f)
    for line in f:
        words = line.split()
        curFile[words[0][3:]].append(int(words[5]))
    controls[refFile] = curFile
    f.close()

topRanks = [('',-1,sys.maxint)] * 10 # plenty of spots to go around...

bd1 = [1,3,15,14,2,0]
bd2 = [7,9,5,4,8,6]
bd3 = [13,15,17,16,14,12]
bd4 = [19,21,23,22,20,19]
bd5 = [25,27,29,28,26,24]
bd6 = [31,33,35,34,32,30]

jChroms = dict()
jChroms['chr1'] = bd1
jChroms['chr2'] = bd2
jChroms['chr3'] = bd3
jChroms['chr4'] = bd4
jChroms['chr5'] = bd5
jChroms['chr6'] = bd6
# Read and process line for each file we have of autosomal chromosomes
for j in range(len(jChroms)-1,-1,-1): 
    print j

#for jChrom in jChroms.keys():
#    updateBestBins(jChroms[jChrom], jChrom)
tempTable = dict()
for jChrom in jChroms.keys():
    tempTable[jChrom] = [1] * 6
    bins = jChroms[jChrom]
    maxPos = 0
    maxPos = topRanks.index(max(topRanks, key=lambda x : x[2]))
    for i in range(len(bins)):      # for each reference bin
        # Zero bins are marked by -1, ignore them
        if bins[i] >= 0:
            if bins[i] < topRanks[maxPos][2]:
                topRanks[maxPos] = (jChrom, i, bins[i])
                maxPos = topRanks.index(max(topRanks, key=lambda top : top[2]))

#print topRanks, '\n'
#topRanks = sorted(topRanks, key=lambda top : top[1])
#topRanks = sorted(topRanks, key=lambda top : top[2])
print tempTable
# Don't take bins close to eachother, take the best one instead
topRanks = sorted(topRanks, key=lambda x : x[2])
    
tempTop = []
for jbin in topRanks:
    if jbin[0] != '' and jbin[1] >=0:
        if tempTable[jbin[0]][jbin[1]] == 1:
            tempTop.append(jbin)
            if jbin[1] > 0:
                tempTable[jbin[0]][jbin[1]-1] = 0
            if jbin[1] < len(tempTable[jbin[0]]) - 1:
                tempTable[jbin[0]][jbin[1]+1] = 0


print tempTop
'''
'''
parser.add_argument('-new', '--n', action='store_false', default=True
    help='build a new complete reference table if True, otherwise, only recalculates optimal cutoff')

args = parser.parse_args()
print '\n# Settings used:'
argsDict = args.__dict__
argsKeys = argsDict.keys()
argsKeys.sort()
for arg in argsKeys:
    print '\t'.join([arg,str(argsDict[arg])])
print args.n
'''
'''
refTable =dict()
lookUp = dict()
for chrom in range(0,4):
    refTable[chrom] = [[(0,1,0),(0,1,1),(0,1,2),(0,1,3)], \
        [(0,1,0),(0,1,1),(0,1,2),(0,1,3)],[(0,1,0),(0,1,1),(0,1,2),(0,1,3)]]
    lookUp[chrom] = []
for chrom in refTable:
    for ibin in range(0, len(refTable[chrom])):
        lookUp[chrom].append([])
        for tbin in refTable[chrom][ibin]:
            if tbin[2] > 1:
                lookUp[chrom][ibin].append(tbin)
print lookUp
'''


parser.add_argument('sample', type=str,
                    help='directory containing samples (.correct) to be used as reference')
parser.add_argument('refout', type=str,
                    help='reference table output, used for sample testing (pickle)')
args = parser.parse_args()

'''
zScores = [1] * 61928
zSmooth = [0] * 61928


with open(args.refout,'wb') as outfile:
    with open(args.refIn,'rU') as infile:
        header = next(infile)
        header = header.strip()
        outfile.write(header + '\tZ\tsw_Z\n')
        i = 0
        for line in infile:
            i += 1
            line = line.strip()
            outfile.write(line + '\t' + str(zScores[i]) + '\t' + str(zSmooth[i]) + '\n')
'''
chromList = [str(chrom) for chrom in range(1,23)]
chromList.append('X')
chromList.append('Y')

readFreq = dict()
for chrom in chromList:
    readFreq[chrom] = []
try:
    outfile = open(args.refout,'wb')
    with open(args.sample,'rb') as sampleFile:
        next(sampleFile)
        last = ''
        i = 0
        for line in sampleFile:
            words = line.split()
            if words[0] == last:
                if i < 100:
                    outfile.write(line)
                    i += 1
            else:
                outfile.write(line)
                last = words[0]
                i = 0
except IOError as err:
    print 'IOError:' + str(err)
    sys.exit()