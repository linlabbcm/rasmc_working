#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2017 YOUR NAME HERE and Charles Lin lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of YOUR PROJECT HERE




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'rasmc_all'
genome ='rn6'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files


#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Rx
chip_data_file = '%sdata_tables/RASMC_CHIP_DATA_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I, LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)

    #assumes macs has already been run and formatted
#    run_macs(chip_data_file)

#    sys.exit()

    print('\n\n')
    print('#======================================================================')
    print('#======================II. PLOTTING BRD4 TRACKS========================')
    print('#======================================================================')
    print('\n\n')

    # for BRD4
    dataFile = chip_data_file

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_1.gff'

    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    
    
    
    plotName = 'rasmc_all_figure_1_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_2.gff'
    
    plotName = 'rasmc_all_figure_2_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3.gff'

    plotName = 'rasmc_all_figure_3_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')

    figureGFF=[['chr16','Btd','',7671879,7881204,'','.','Btd']]
    figureGFFPath='%sBtd.gff' % (gffFolder)
    utils.unParseTable(figureGFF,figureGFFPath,'\t')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Btd.gff'

    plotName = 'rasmc_all_btd_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/enrich_tfs.gff'

    plotName = 'rasmc_all_enrich_tfs_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')
    

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    plotName = 'rasmc_all_edn1_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')



    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    plotName = 'rasmc_all_edn1_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/regions_fosl1_atf4_spry2_thbs1.gff'

    
    plotName = 'spry2_ext_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/regions_fosl1_atf4_spry2_thbs1.gff'

    plotName = 'spry2_ext_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')





    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'


    plotName = 'Acta2_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'

    plotName = 'Acta2_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'


    plotName = 'cluster_b_goi_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'

    plotName = 'cluster_b_goi_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')



    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'


    plotName = 'figure_3_supp_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'

    plotName = 'figure_3_supp_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')








    print('\n\n')
    print('#======================================================================')
    print('#======================III. ENHANCER PROMOTER==========================')
    print('#======================================================================')
    print('\n\n')

    dataDict = pipeline_dfci.loadDataTable(chip_data_file)
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

    allLoci = []

#    for name in namesList:

#          collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#          allLoci += collection.getLoci()

    #do this for each one in the namesList
    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_BRD4_STITCHED_-0_+0.gff','\t')







    dataFile = chip_data_file
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    gff = '%sRN6_RASMC_BRD4_STITCHED_-0_+0.gff' % (gffFolder)
    activityTable = '%sactiveListTable.txt' % (tableFolder)
    outputFolder = '%senhancerPromoter/BRD4/' % (projectFolder)
    desc_string = '0_STITCH_-_JQ1'

#    makeEnhProBash(dataFile,gff,activityTable,namesList,outputFolder,desc_string)


    annotTable = utils.parseTable(annotFile,'\t')
    dictFile = '/storage/cylin/bin/pipeline/crc/annotation/MotifDictionary.txt'
    motifs = utils.parseTable(dictFile,'\t')
    
    TFlist = []
    
#    for line in annotTable[1:]:
#        for motif in motifs:
#            if line[12].upper()==motif[1]:
#                NMID = line[1]
#                name = line[12].upper()
#                new_line=[NMID,name]
#                TFlist.append(new_line)
    
#    utils.unParseTable(TFlist,'/storage/cylin/grail/projects/rasmc_all/tables/TFlist.txt','\t')

    print('\n\n')
    print('#======================================================================')
    print('#======================IV. TF NETWORK ANALYSIS=========================')
    print('#======================================================================')
    print('\n\n')


    #==========================================================================
    #===================ENH TF GFF=============================================
    #==========================================================================
    print('Creating Brd4 Enhancer TF in degree table from crc rasmc_h3k27ac_0_tss_ENHANCER_TF_TABLE.txt')
    
#    EnhTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_ENHANCER_TF_TABLE.txt' % (projectFolder)

#    EnhTable = utils.parseTable(EnhTableFile,'\t')
#    genesOfInterest=[]

#    enhGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region   
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(EnhTable[1])

#    for line in EnhTable:
#        if line != EnhTable[0]:
#            gene = line[4]
#            ID = line[0]
#            chrom = line[1]
#            TSS = ''
#            start = line[2]
#            stop = line[3]
#            line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,gene]
#            enhGFF.append(line_gff_format)

#    enhGFF_Path = '%srasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff' % (gffFolder,)
#    enhGFFTable = utils.unParseTable(enhGFF, enhGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to Enhancer TF in degree gff')
#    dataFile = chip_data_file
#    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff']
#    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making signal table for Brd4 Enhancer TF in degree')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE_signal_table.txt')

    #==========================================================================
    #===================EDGE TABLE GFF=========================================
    #==========================================================================
    print('Creating Edge Table GFF crc rasmc_h3k27ac_0_tss_EDGE_TABLE.txt')



#    peakTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_EDGE_TABLE.txt' % (projectFolder)

#    clusterTable = utils.parseTable(peakTableFile,'\t')
#    genesOfInterest=[]

#    clusterGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(clusterTable[1])

#    for line in clusterTable:
#        if line != clusterTable[0]:
#            gene = line[0]
#            ID = line[5]
#            chrom = line[2]
#            TSS = ''
#            start = line[3]
#            stop = line[4]
#            line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,'']
#            clusterGFF.append(line_gff_format)

#    clusterGFF_Path = '%srasmc_h3k27ac_0_tss_EDGE_TABLE.gff' % (gffFolder,)
#    clusterGFFTable = utils.unParseTable(clusterGFF, clusterGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to Edge Table gff')
    dataFile = chip_data_file
#    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_EDGE_TABLE.gff']
#    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making brd4 signal table for Edge Table gff')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_EDGE_TABLE.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_EDGE_TABLE_Brd4_signal_table.txt')

    #==========================================================================
    #===================EDGE TABLE GFF=========================================
    #==========================================================================
    print('Creating crc subpeak Table')



    subpeakTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_all_subpeak.bed' % (projectFolder)

    subpeakTable = utils.parseTable(subpeakTableFile,'\t')
   
    subpeakGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(subpeakTable[1])

#    for line in subpeakTable:
#        gene = ''
#        ID = line[4]
#        chrom = line[0]
#        TSS = ''
#        start = line[1]
#        stop = line[2]
#        line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,'']
#        subpeakGFF.append(line_gff_format)
    
#    print(subpeakGFF[1])
#    subpeakGFF_Path = '%srasmc_h3k27ac_0_tss_subpeak_table.gff' % (gffFolder,)
#    subpeakGFFTable = utils.unParseTable(subpeakGFF, subpeakGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to crc subpeak Table gff')
    dataFile = chip_data_file
    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_subpeak_table.gff']
    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making brd4 signal table for crc subpeak Table gff')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_subpeak_table.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_subpeak_Brd4_signal_table.txt')


    print('\n\n')
    print('#======================================================================')
    print('#======================IV. POST ANALYSIS===============================')
    print('#======================================================================')
    print('\n\n')





    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to crc subpeak Table gff')
    dataFile = chip_data_file
    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff']
    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW',
'RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)




    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making signal table for Brd4 Enhancer TF in degree')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/RASMC_BRD4_across_0_stitch_h3k27ac_regions.txt')



       
#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RUNNING MACS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def run_macs(dataFile):
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    namesList.sort()
    print(namesList)
    pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')
    os.chdir(projectFolder) # the silly call macs script has to change into the output dir
    #so this takes us back to the project folder

    #to check for completeness, we will try to find all of the peak files
    peak_calling_done = False
    while not peak_calling_done:
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        for name in namesList:
            peak_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
            print('searching for %s' % (peak_path))
            if utils.checkOutput(peak_path,1,180):
                peak_calling_done =True
                print('found %s' % (peak_path))
                continue
            else:
                print('Error: peak calling timed out')
                sys.exit()
    
    #now format the macs output
    print('formatting macs output')
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink ='',useBackground=True)
    print('Finished running Macs 1.4.2')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKING ENHANCER PROMOTER BASH SCRIPTS~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def makeEnhProBash(dataFile,gff='',activityTable='',namesList=[],outputFolder='',desc_string=''):
    dataDict = pipeline_dfci.loadDataTable(dataFile)

    print(dataDict)
    if len(namesList) == 0:
        namesList = dataDict.keys()
    namesList.sort()
    print(namesList)
    
    for name in namesList:
        bashFileName = "%s%s_%s_enhPro.sh" % (outputFolder,name,desc_string)
        print('Writing %s') % bashFileName
        bashFile = open(bashFileName,'w')

        #shebang
        bashFile.write('#!/usr/bin/bash\n')

        line = 'BAMS=\'%s\'' % (dataDict[name]['bam'])
        bashFile.write(line+'\n')

        line = 'CONTROLS=\'%s\'' % (dataDict[name]['background'])
        bashFile.write(line+'\n\n')

        line = 'OUTPUT=\'./\''
        bashFile.write(line+'\n')

        line = 'INPUT=\'%s\'' % (gff)
        bashFile.write(line+'\n')
        line = 'ACTIVITY=\'%s\'' % (activityTable)
        bashFile.write(line+'\n')
        line = 'NAME=\'%s_%s\'' % (name,desc_string)
        bashFile.write(line+'\n\n\n')

        line = 'python /storage/cylin/bin/pipeline/enhancerPromoter.py -b $BAMS -g %s -i $INPUT -o $OUTPUT -a $ACTIVITY --name $NAME' % (dataDict[name]['genome'])
        bashFile.write(line+'\n')
        
        bashFile.close()




#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
