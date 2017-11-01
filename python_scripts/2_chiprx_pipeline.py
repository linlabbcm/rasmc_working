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
chiprx_data_file = '%sdata_tables/RASMC_H3K27AC_CHIPRX_ALL.txt' % (projectFolder)




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
    pipeline_dfci.summary(chiprx_data_file)

    #assumes macs has already been run and formatted
    #run_macs(chiprx_data_file)

    #sys.exit()
    print('\n\n')
    print('#======================================================================')
    print('#=====================II. MAPPING H3K27AC REGIONS======================')
    print('#======================================================================')
    print('\n\n')

    #assumes macs has already been run and formatted
    #run_macs(chiprx_data_file)
    
    #next makes gffs of the k27ac regions individually, plus of intersect and union
    set_name = 'RASMC_H3K27AC'
    names_list = ['RASMC_H3K27AC_UNSTIM_NEW',
                  'RASMC_H3K27AC_PDGF_2H_NEW',
                  'RASMC_H3K27AC_PDGF_24H_NEW',
                  'RASMC_H3K27AC_UNSTIM_REP1',
                  'RASMC_H3K27AC_PDGF_2H_REP2',
                  'RASMC_H3K27AC_PDGF_24H_REP2',
                  ]

    #gff_list= makeStitchedGFF(chiprx_data_file,set_name,names_list,write_individual=True)
    #print(gff_list)
    #sys.exit()

    gff_list =  ['%sRN6_RASMC_H3K27AC_UNSTIM_NEW_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_UNSTIM_REP1_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_INTERSECT_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_UNION_-0_+0.gff' % (gffFolder),
                ]
    #print(gff_list)
    #signal_table_list=pipeline_dfci.map_regions(chiprx_data_file,gff_list,mappedFolder,signalFolder,names_list=[],medianNorm=False,output='')
    #print(signal_table_list)
    
    #sys.exit()
    
    print('\n\n')
    print('#======================================================================')
    print('#================III. CALCULATING CHIPRX SCALE FACTORS=================')
    print('#======================================================================')
    print('\n\n')

    scale_table_path = '%sRN6_RASMC_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    #scale_table_table = writeScaleFactors(chiprx_data_file,namesList=[],primary_genome='rn6',rx_genome='dm6',output=scale_table_path)

    #sys.exit()

    print('\n\n')
    print('#======================================================================')
    print('#======================IV. MAKING BOX PLOTS============================')
    print('#======================================================================')
    print('\n\n')



    #making boxplots
    boxplot_script_path = '%sr_scripts/4_chiprx_plots.R' % (projectFolder)
    output_folder = utils.formatFolder('%sfigures/4_chiprx_plots/' % (projectFolder),True)
    scale_table_path = '%sRN6_RASMC_NEW_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    #=============================================================================

    #for new unstim h3k27ac regions    
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNSTIM_NEW'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNSTIM_NEW_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)

    #for old unstim h3k27ac regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNSTIM'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNSTIM_REP1_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)

    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)


    #for rasmc intersect regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_INTERSECT'
    signal_table_path = '%sRN6_RASMC_H3K27AC_INTERSECT_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)


    #for rasmc union regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNION'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNION_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)

#   sys.exit()


    print('\n\n')
    print('#======================================================================')
    print('#======================V. MAKING RN6 GFFS=============================')
    print('#======================================================================')
    print('\n\n')


#    pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species='RN6')


    print('\n\n')
    print('#======================================================================')
    print('#======================VI. MAPPING -JQ1 BAMS TO TSS REGIONS============')
    print('#======================================================================')
    print('\n\n')
    
    dataFile = chiprx_data_file
    gffList = ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_TSS_ALL_-1000_+1000.gff']
    setName = 'RASMC_H3K27AC_-JQ1'
    cellTypeList = ['RASMC']
    enrichedFolder = macsEnrichedFolder
#    namesList=['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    pipeline_dfci.mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,macs=True,namesList=['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW'],useBackground=True)


#    setList = [['RASMC_H3K27AC_UNSTIM_REP1'],['RASMC_H3K27AC_UNSTIM_NEW'],['RASMC_H3K27AC_PDGF_2H_REP2'],['RASMC_H3K27AC_PDGF_2H_NEW'], ['RASMC_H3K27AC_PDGF_24H_REP2'],['RASMC_H3K27AC_PDGF_24H_NEW']]
#    mappedEnrichedFile = '%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_RASMC_H3K27AC_-JQ1.txt' % (projectFolder)
#    output = '%stables/RN6_H3K27AC_-JQ1_TSS_ACTIVE_-1000_+1000.txt' % (projectFolder)
#    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)






    print('\n\n')
    print('#======================================================================')
    print('#======================VII. MAKING ACTIVITY TABLE=======================')
    print('#======================================================================')
    print('\n\n')
    
    h3k27ac_mapped_path = '%stables/RN6_H3K27AC_-JQ1_TSS_ACTIVE_-1000_+1000.txt' % (projectFolder)
    expr_path = '%scufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt' % (projectFolder)
    
#    activeTable = makeActiveTable(h3k27ac_mapped_path,expr_path)


    print('\n\n')
    print('#======================================================================')
    print('#======================VIII. MAKING H3K27AC TRACKS=====================')
    print('#======================================================================')
    print('\n\n')





#    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_1.gff'
#    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_2.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


#    plotName = 'rasmc_all_figure_1_h3k27ac_tracks_merge'
    plotName= 'rasmc_all_figure_2_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'
    
    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'
    
#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString =scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_figure_3_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Btd.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_btd_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)




    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_edn1_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'

    dataFile = chiprx_data_file


    plotName= 'rasmc_all_Acta2_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/spry2_20kb_ext.gff'

    dataFile = chiprx_data_file


    plotName= 'spry2_ext_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)




    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'

    dataFile = chiprx_data_file


    plotName= 'cluster_b_goi_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'

    dataFile = chiprx_data_file


    plotName= 'figure_3_supp_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)


    print('\n\n')
    print('#======================================================================')
    print('#======================IX. MAKING 0 STITCHED GFFS======================')
    print('#======================================================================')
    print('\n\n')


#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

    #do this for each one in the namesList
    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff','\t')

#######################Making Timepoint GFFs############################################

    ##Unstim##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_UNSTIM_NO_JQ1_STITCHED_-0_+0.gff','\t')

    ##2H+PDGF##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_2H_NO_JQ1_STITCHED_-0_+0.gff','\t')


    ##24H+PDGF##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_24H_NO_JQ1_STITCHED_-0_+0.gff','\t')



    print('\n\n')
    print('#======================================================================')
    print('#======================X. ENHANCER PROMOTER============================')
    print('#======================================================================')
    print('\n\n')

    dataFile = chiprx_data_file
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]

    gff = '%sRN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff' % (gffFolder)
    activityTable = '%sactiveListTable.txt' % (tableFolder)
    outputFolder = '%senhancerPromoter/H3K27AC/' % (projectFolder)
    desc_string = '0_STITCH_-_JQ1'

#    makeEnhProBash(dataFile,gff,activityTable,namesList,outputFolder,desc_string)


    print('\n\n')
    print('#======================================================================')
    print('#======================XI. CALLING ROSE2===============================')
    print('#======================================================================')
    print('\n\n')

    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]

#    parentFolder = utils.formatFolder('%srose' % (projectFolder),True)
#    bashFileName = '%sRASMC_ENHANCER_H3K27AC_0_tss_rose.sh' % (parentFolder)


#    pipeline_dfci.callRose2(chiprx_data_file,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=0,stitch='',bashFileName =bashFileName,mask='')


    print('\n\n')
    print('#======================================================================')
    print('#======================XI. CLUSTERING ROSE OUTPUT======================')
    print('#======================================================================')
    print('\n\n')




    #load the datadict
    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)

    #analysis name
    analysisName = 'rasmc_h3k27ac_0_tss'

    #nameslist
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 ]

    namesString = string.join(namesList,',')


    print namesList


    #set up the output
    outputFolder = '%sclustering/%s/' % (projectFolder,analysisName)
    outputFolder = utils.formatFolder(outputFolder,True)

    #get the rose folder
    roseFolder ='%srose/' % (projectFolder)

    #set up the bash file
    bashFileName = '%s%s_0_tss_clustering.sh' % (outputFolder,analysisName)
    bashFile = open(bashFileName,'w')

    bashFile.write('#!/usr/bin/bash\n')

    #for now change into pipelinedir just to be safe
    bashFile.write('cd /storage/cylin/bin/pipeline/\n')

    #if you want to do w/ all do this
    #this will run w/ default parameters
    #see the documentation for available flags
    clusterCmd = 'python /storage/cylin/bin/pipeline/clusterEnhancer.py -d %s -i %s -r %s -o %s -n %s' % (chiprx_data_file,namesString,roseFolder,outputFolder,analysisName)

    bashFile.write(clusterCmd + '\n')
    print(clusterCmd)

    bashFile.close()


    
    
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~WRITING SCALE FACTORS~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def writeScaleFactors(dataFile,namesList=[],primary_genome='hg19',rx_genome='dm6',output=''):

    '''
    creates a table of scale factors based on the rx genome read depth
    '''
    
    #first set up the output folder
    #rpm scale factor is what the rpm/bp should be MULTIPLIED by
    #mouse mapped reads give the denominator for what raw r/bp should be divided by
    outputTable = [['NAME','RN6_MAPPED_READS','DM6_MAPPED_READS','RPM_SCALE_FACTOR']]

    
    dataDict=pipeline_dfci.loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = [name for name in dataDict.keys()]
    namesList.sort()
    print('scaling the following datasets')


    for name in namesList:
        
        print('WORKING ON %s' % (name))
        bam_path = dataDict[name]['bam']
        bam = utils.Bam(bam_path)
        bam_mmr = float(bam.getTotalReads())/1000000
        scale_path = string.replace(bam_path,primary_genome,rx_genome)
        scaleBam = utils.Bam(scale_path)
        scale_mmr = float(scaleBam.getTotalReads())/1000000
        #print(bam_path)
        #print(scale_path)
        rpm_scale = bam_mmr/scale_mmr
        scale_line = [bam_mmr,scale_mmr,rpm_scale]
        scale_line = [round(x,4) for x in scale_line]
        outputTable.append([name] + scale_line)

    if len(output) == 0:
        return outputTable
    else:
        utils.unParseTable(outputTable,output,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~SCALING WIGGLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def scaleWiggles(dataFile,scaleTableFile,names_list=[]):

    '''
    first unzips each wiggle
    then scales each line by the rpm
    and rounds to a reasonable number (2 decimal)
    '''

    dataDict=pipeline_dfci.loadDataTable(dataFile)
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('WCE') ==0 and name.count('INPUT') == 0]
    names_list.sort()
    print(names_list)

    print('loading scale factors')

    scale_table = utils.parseTable(scaleTableFile,'\t')
    scale_dict = {}
    for line in scale_table[1:]:
        scale_dict[line[0]] = float(line[2])
    os.chdir(wiggleFolder)
                   
    for name in names_list:

        print('scaling %s' % (name))
        scale_factor = scale_dict[name]
        wig_path_gz = '%swiggles/%s_treat_afterfiting_all.wig.gz' % (projectFolder,name)
        wig_path = '%swiggles/%s_treat_afterfiting_all.wig' % (projectFolder,name)
        

        wig_out = '%swiggles/%s_scaled.wig' % (projectFolder,name)
        wig_out_final ='%swiggles/%s_scaled.wig.gz' % (projectFolder,name)
        if utils.checkOutput(wig_out_final,0,0):
            print('Found scaled wiggle for %s at %s' % (name,wig_out_final))
            continue
        cmd = 'gunzip %s' % (wig_path_gz)
        print(cmd)

        #this should run to completion
        os.system(cmd)

        #now open up the new wig
        wig = open(wig_path,'r')
        wig_scaled = open(wig_out,'w')

        ticker = 0
        for line in wig:

            if ticker % 1000000 == 0:
                print(ticker)
            ticker+=1
            if line[0] == 't' or line[0] == 'v':
                wig_scaled.write(line)
            else:
                line = line.rstrip().split('\t')
                line[1] = str(round(float(line[1])/scale_factor,2))
                line_string = '\t'.join(line) + '\n'
                wig_scaled.write(line_string)


        wig.close()
        wig_scaled.close()
        cmd = 'gzip %s' % (wig_out)
        print(cmd)
        os.system(cmd)

        cmd = 'gzip %s' % (wig_path)
        print(cmd)
        os.system(cmd)
    os.chdir(projectFolder)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING REGION GFFS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeStitchedGFF(dataFile,set_name,names_list,write_individual=False):

    '''
    makes a stitched gff and dumps it into the gff folder
    '''

    dataDict = pipeline_dfci.loadDataTable(dataFile)
    
    loci = []
    collection_dict = {}
    gff_list = []
    for name in names_list:
        print(name)
        macsEnrichedFile = '%s%s_peaks.bed' % (macsEnrichedFolder,name)
        collection = utils.importBoundRegion(macsEnrichedFile,name)
        collection_dict[name]=collection
        loci+= collection.getLoci()
        if write_individual:
            gff_path = '%s%s_%s_-0_+0.gff' % (gffFolder,genome.upper(),name)
            gff = utils.locusCollectionToGFF(collection)
            utils.unParseTable(gff,gff_path,'\t')
            gff_list.append(gff_path)

    all_collection = utils.LocusCollection(loci,50)
    stitched_collection = all_collection.stitchCollection()

    gff = utils.locusCollectionToGFF(stitched_collection)

    out_path = '%s%s_%s_UNION_-0_+0.gff' % (gffFolder,genome.upper(),set_name)
    print(out_path)
    utils.unParseTable(gff,out_path,'\t')
    gff_list.append(out_path)

    #now get the intersect gff
    print('getting intersection gff')
    stitched_loci = stitched_collection.getLoci()
    intersect_loci = []
    ticker = 0
    for locus in stitched_loci:
        if ticker%1000==0:
            print(ticker)
        ticker+=1
        overlap = True
        for name in names_list:
            if len(collection_dict[name].getOverlap(locus,'both')) == 0:
                overlap = False

        if overlap == True:
            intersect_loci.append(locus)

            
    print('identified %s stitched loci' % (len(stitched_loci)))
    print('identified %s intersect loci' % (len(intersect_loci)))

    intersect_collection = utils.LocusCollection(intersect_loci,50)

    intersect_gff = utils.locusCollectionToGFF(intersect_collection)

    intersect_path = '%s%s_%s_INTERSECT_-0_+0.gff' % (gffFolder,genome.upper(),set_name)
    print(intersect_path)
    utils.unParseTable(intersect_gff,intersect_path,'\t')
    gff_list.append(intersect_path)
    return gff_list



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING ACTIVE GENE LIST~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def makeActiveTable(h3k_mapped='',expr_table=''):
    #Takes in path of h3k27ac mapped enriched file for all h3k27ac bams at the TSS +- 1000 bp
    #and the path to the rna expression fpkm means table
    #Any gene found in both is included in the list of active genes

    h3k27ac_mapped_Table = utils.parseTable(h3k_mapped,'\t')
    fpkm_expr_table = utils.parseTable(expr_table,'\t')
    activeList = []

    for x in h3k27ac_mapped_Table:
        for y in fpkm_expr_table:
            if y[0] == x[2]:
                new_line = [x[1],x[2]]
                activeList.append(new_line)

    out_path = '%stables/activeListTable.txt' % (projectFolder)
    activeListTable = utils.unParseTable(activeList,out_path,'\t')
    print('active table located at %stables/activeListTable.txt') % (projectFolder)


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
