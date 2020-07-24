#! python2.7

######################
#
# Core Regulatory Circuits
# Young and Bradner Labs
# Version 1.0
# 140724
#
######################

######################           
# Dependencies       
######################

import os
import sys
#sys.path.append('/ark/home/af661/src/utils/')
sys.path.append('/home/rob/Documents/tools/ChIPSeq/young_computation-crcmapper-a547e044fb7f/')
import utils

import string

import numpy
import scipy
import scipy.stats

import subprocess
import os

from string import upper
from random import randrange
from collections import defaultdict

import networkx as nx
from networkx.algorithms.clique import find_cliques_recursive
import pickle


def findValleys(regionsLoci, bamFile, projectName, projectFolder, cutoff = 0.2):
    '''
    takes in the super dict
    returns a dictionary of refseqs with all valley loci that are associated
    '''

    print 'IDENTIFYING VALLEYS IN SUPER ENHANCERS'

    valleyBED = []
    valleyDict = {}

    locusList = regionsLoci.getLoci()
    for region in locusList:
    
        scoreArray = scoreValley(region, bamFile, projectName, projectFolder)
        for index,score in enumerate(scoreArray):
            if score > cutoff:
                valley = utils.Locus(region.chr(), region.start() + index*10,
                                        region.start() + (index+1)*10, '.')
                valleyBED.append([valley.chr(), valley.start(), valley.end()])

    bedfilename = projectFolder + projectName + '_valleys.bed'
    utils.unParseTable(valleyBED, bedfilename, '\t')
    print bedfilename

    return bedfilename


def scoreValley(locus, bamFile, projectName, projectFolder):
    '''
    calculate valley scores for a locus
    based on this refernce:
    http://bioinformatics.oxfordjournals.org/content/26/17/2071.full
    '''

    nbins = locus.len()/10

    #call bamliquidator on the region and store in a temp file
    os.system('bamliquidator ' + bamFile + ' ' + locus.chr() + ' ' + str(locus.start()) + ' '
              + str(locus.end()) + ' . ' + str(nbins) + ' 0 > ' + projectFolder + 'tempBamliquidator_'
              + projectName + '.txt')

    x = utils.parseTable(projectFolder + 'tempBamliquidator_' + projectName + '.txt', '\t')
    density = [int(y[0]) for y in x]
    smoothDensity =  gaussianSmooth(density, 5)

    scoreArray = []
    regionMax = max(smoothDensity)

    #Now take the smooth reads and calaculate a valley score

    for i in range(len(smoothDensity)):
        score = 0
        try:
            leftmax = max(smoothDensity[i-25:i-10])
        except:
            leftmax = 'edge'
        try:
            rightmax = max(smoothDensity[i+10:i+25])
        except:
            rightmax = 'edge'

        if rightmax == 'edge' and leftmax == 'edge':
            shoulderHeightMin = 0
            shoulderHeightMax = 0
        elif leftmax == 'edge':
            shoulderHeightMin = rightmax
            shoulderHeightMax = rightmax
        elif rightmax == 'edge':
            shoulderHeightMin = leftmax
            shoulderHeightMax = leftmax
        else:
            shoulderHeightMin = min(leftmax, rightmax)
            shoulderHeightMax = max(leftmax, rightmax)

        ratio = (shoulderHeightMax-float(smoothDensity[i]))/regionMax
        if ratio > 0.3:
            score = 1
        else:
            score = 0

        scoreArray.append(score)

    return scoreArray


def stitchValleys(valleyList):
    '''
    takes a list of valley loci
    returns a stitched list of valleys to extract seq from
    '''

    valleyCollection = utils.LocusCollection(valleyList,1)
    stitchedValleyCollection = valleyCollection.stitchCollection()
    loci = []
    regions = []
    for valley in stitchedValleyCollection.getLoci():
        if [valley.chr(), valley.start(), valley.end()] not in regions:
            loci.append(valley)
            regions.append([valley.chr(), valley.start(), valley.end()])
    return loci

def gaussianSmooth(readList, degree=5):
    '''
    Smoothing function for raw bamliquidator output
    '''

    window=degree*2-1
    weight=numpy.array([1.0]*window)
    weightGauss=[]

    for i in range(window):

        i = i-degree+1
        frac = i/float(window)
        gauss = 1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)

    weight=numpy.array(weightGauss)*weight
    smoothed=[0.0]*(len(readList)-window)

    for i in range(len(smoothed)):
        smoothed[i]=sum(numpy.array(readList[i:i+window])*weight)/sum(weight)

    smoothed = [0,0,0,0,0] + smoothed + [0,0,0,0] # return an array of the same length

    return smoothed

def generateSubpeakFASTA(regionsLoci, subpeaks, genomeDirectory, projectName, projectFolder, constExtension):
    '''
    from a BED file of constituents
    generate a FASTA for the consituients contained within the canidate supers
    '''

    subpeakBED = [['track name=' + projectName + ' color=204,0,204']]
    subpeakTable = utils.parseTable(subpeaks, '\t')

    subpeakLoci = [utils.Locus(l[0], int(l[1]), int(l[2]), '.') for l in subpeakTable]
    subpeakCollection = utils.LocusCollection(subpeakLoci, 50)


    for region in regionsLoci.getLoci():
        overlaps = subpeakCollection.getOverlap(region)
        extendedOverlaps = [utils.makeSearchLocus(x, constExtension, constExtension) for x in overlaps]

        overlapCollectionTemp = utils.LocusCollection(extendedOverlaps, 50)
        overlapCollection = overlapCollectionTemp.stitchCollection()
        for overlap in overlapCollection.getLoci():
            subpeakBED.append([overlap.chr(), overlap.start(), overlap.end(), '.', region.ID()])


    bedfilename = projectFolder + projectName + '_subpeaks.bed'
    utils.unParseTable(subpeakBED, bedfilename, '\t')

    fasta = []
    for subpeak in subpeakBED[1:]:

        fastaTitle = subpeak[4] + '|'  + str(subpeak[0]) + '|' + str(subpeak[1]) + '|' + str(subpeak[2])
        fastaLine = utils.fetchSeq(genomeDirectory, subpeak[0], int(subpeak[1]+1), 
                                    int(subpeak[2]+1))

        fasta.append('>' + fastaTitle)
        fasta.append(upper(fastaLine))

    outname = projectFolder + projectName + '_SUBPEAKS.fa'

    utils.unParseTable(fasta, outname, '')

def findMotifs(projectFolder, projectName, motifConvertFile, motifDatabaseFile):
    '''
    takes the refseq to subpeak seq dict
    returns the networkx object with all connections
    '''

    # Create a dictionary to call motif names keyed on gene names

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]
    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    bgCmd = 'fasta-get-markov -m 1 < ' + projectFolder + projectName + '_SUBPEAKS.fa > ' + projectFolder + projectName + '_bg.meme'
    subprocess.call(bgCmd, shell=True)

    utils.formatFolder(projectFolder + 'FIMO/', True)

    fimoCmd = 'fimo'
    fimoCmd += ' -verbosity 1'  # thanks for that ;)!
    fimoCmd += ' -text'
    fimoCmd += ' -oc ' + projectFolder + 'FIMO'
    fimoCmd += ' --bgfile ' + projectFolder + projectName + '_bg.meme'
    fimoCmd += ' ' + motifDatabaseFile + ' '
    fimoCmd += projectFolder + projectName + '_SUBPEAKS.fa'
    fimoCmd += ' > '+ projectFolder + 'FIMO/fimo.txt'  ##
    print fimoCmd

    fimoOutput = subprocess.call(fimoCmd, shell=True)  #will wait that fimo is done to go on

    return fimoCmd

###################### 
#
# Main Method
#
######################

def main():

    from optparse import OptionParser

    usage = "usage: %prog [options] -e [ENHANCER_FILE] -b [BAM_FILE] -g [GENOME] -o [OUTPUTFOLDER] -n [NAME]" 
    parser = OptionParser(usage = usage)

    #required flags                                                                                                                                        
    parser.add_option("-e","--enhancer_file", dest="enhancers",nargs = 1, default=None,
                      help = "Provide a ROSE generated enhancer table (_AllEnhancers.table.txt)")
    parser.add_option("-B","--bed_file", dest="regions",nargs = 1, default=None,
                      help = "Provide a bed file (.bed)")                      
    parser.add_option("-b","--bam",dest="bam",nargs =1, default = None,
                      help = "Provide a bam that corresponds to the super enhancer table")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "Provide the build of the genome to be used for the analysis. Currently supports HG19, HG18 and MM9")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "Enter an output folder")
    parser.add_option("-n","--name",dest="name",nargs =1, default = None,
                      help = "Provide a name for the job")
    parser.add_option("-k","--background", dest="background",nargs = 1, default=None,
                      help = "Provide a background BAM file")

    #additional options

    parser.add_option("-s","--subpeaks", dest="subpeaks",nargs=1,default=None,
                      help = "Enter a BED file of regions to search for motifs")
    parser.add_option("-l","--extension-length", dest="extension",nargs = 1, default=100,
                      help = "Enter the length to extend subpeak regions for motif finding")
    parser.add_option("-E","--enhancer_number", dest="Enumber",nargs = 1, default=None,
                      help = "Enter the number of top ranked enhancers to include in the anlaysis. Default is all regions. you can submit 'super'")
    parser.add_option("-N", "--number", dest="number",nargs = 1, default=2,
                      help = "Enter the number of motifs required to assign a binding event")     #I have modified the destination of -N option so that it is different from the destination of -E option
    parser.add_option("--motifs", dest="motifs",nargs = 1, default=False,
                      help = "Enter an alternative PWM file for the analysis")
    parser.add_option("-t","--tfs", dest="tfs",nargs=1,default=None,
                      help = "Enter additional TFs (comma separated) to be used in the bindinf analysis")

    (options,args) = parser.parse_args()

    print(options)

    if (options.enhancers or options.regions) and options.genome and options.output and options.name:

        ###
        # Define all global file names
        ###
        
        if options.motifs:
            motifDatabaseFile = options.motifs
        else:
            motifConvertFile = '/home/rob/Documents/tools/ChIPSeq/CLL_TFnetworks_2018/annotations/MotifDictionary.txt'
            motifDatabaseFile = '/home/rob/Documents/tools/ChIPSeq/CLL_TFnetworks_2018/annotations/VertebratePWMs.txt'
        
        
        if options.bam:
            bamFile = options.bam
            bam = utils.Bam(bamFile)

        if options.background:
            background = options.background

        else: 
            background = None

        
        genome = options.genome
        genome = upper(genome)
        if genome == 'HG19':
            #genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
            genomeDirectory = '/home/rob/shares/epistress/Bioinfo/Annotations/human/Homo_sapiens/NCBI/build37.2/Sequence/Chromosomes/'
            #annotationFile = '/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc'
            annotationFile = '/home/rob/shares/epistress/Bioinfo/Annotations/human/Homo_sapiens/NCBI/build37.2/Annotation/Archives/archive-2015-07-17-14-32-01/Genes/refFlat_refseq.ucsc'
            #TFfile = '/ark/home/af661/src/coreTFnetwork/annotations/TFlist_NMid_hg19.txt'
            TFfile = '/home/rob/shares/epistress/Bioinfo/Annotations/human/Homo_sapiens/NCBI/build37.2/Annotation/Archives/archive-2015-07-17-14-32-01/Genes/TFList_NMid_hg19.txt'

        if genome == 'HG18':
            genomeDirectory = '/grail/genomes/Homo_sapiens/human_gp_mar_06_no_random/fasta/'
            annotationFile = '/ark/home/cl512/src/pipeline/annotation/hg18_refseq.ucsc'
            TFfile = '/ark/home/af661/src/coreTFnetwork/annotations/TFlist_NMid_hg19.txt'

        if genome == 'MM9':
            genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
            annotationFile = '/ark/home/cl512/pipeline/annotation/mm9_refseq.ucsc'
            TFfile = '/ark/home/af661/src/coreTFnetwork/annotations/TFlist_NMid_mm9.txt'


        projectFolder = options.output
        projectName = options.name

        if options.subpeaks:
            subpeakFile = options.subpeaks
        else: subpeakFile = None
        
        
        constExtension = int(options.extension)

        enhancerNumber = options.Enumber
        if options.Enumber != 'super':
            enhancerNumber = options.Enumber
        else:
            enhancerNumber = 'super'

        if options.regions:
            regionsLoci = utils.bedToLocusCollection(options.regions,top = options.Enumber)
        elif options.enhancers:
            regionsLoci = utils.makeSECollection(enhancerFile=options.enhancers,top = options.Enumber)


        if subpeakFile == None:
            subpeakFile = findValleys(regionsLoci, bamFile, projectName, projectFolder, cutoff = 0.2)
            
        generateSubpeakFASTA(regionsLoci, subpeakFile, genomeDirectory, projectName, projectFolder, constExtension)        
        subpeakFile = projectFolder + projectName + '_SUBPEAKS.fa'

        findMotifs(projectFolder, projectName, motifConvertFile, motifDatabaseFile)
        

    else:
        parser.print_help()
        sys.exit()



if __name__ == '__main__':
    main()