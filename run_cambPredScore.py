#!/usr/bin/env python

# source /proj/relibs/relib00/conda/bin/activate
# echo $PATH
# source /proj/relibs/relib00//conda/env/mypy3 #./conda/env/mypy3
# source /proj/relibs/relib00/conda/etc/profile.d/conda.sh
# source activate mypy3
# cd netZooPy
# python
"""
Description:
  Run predScore after benchmark.

Usage:
  -h, --help: help
  -i, --indir: directory where intersections are location
  -o, --outdir
  -c, --cell: if specific cell line wanted
  -t, --TF: if specific TF wanted
  bash data/MotifPipeline/ENCODE/sthlm_motif_pipeline_beta.sh -b50 -c'A549 K562 GM12878 SKNSH HepG2 HeLa'

Example:
  source /proj/relibs/relib00/conda/bin/activate
  source activate mypy3
  python netZooPy/milipeed/benchmark/run_predScore.py -i data/MotifPipeline/sthlm_motif_5_QCbeta -o data/MotifPipeline/sthlm_motif_5_QCbeta/red/test
"""

import sys,os
sys.path.append("/udd/redmo/netZooPy") 
import getopt
import tests
import shutil
# from netZooPy.lioness.lioness import Lioness
from cambPredScore import cambPredScore#, camb_plotPredScore#, camb_allPredScore

# from netZooPy.milipeed.benchmark.predRegion import predRegion
# from netZooPy.milipeed.benchmark.buffer_distr_comp import buffer_distr_comp
# from netZooPy.milipeed.benchmark.plot_predScore import plot_predScore
# from netZooPy.milipeed.benchmark.plot_allPredScore import plot_allPredScore



def main(argv):
    #Create variables
    indir = None
    outdir = None
    cell = None
    TF = None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['help', 'indir=','outdir='])
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        print(__doc__)
        return 2

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            return 0
        elif opt in ('-i', '--indir'):
            indir = arg
        elif opt in ('-o', '--outdir'):
            outdir = arg
        # elif opt in ('-c', '--cell'):
        #     cell = arg
        # elif opt in ('-t','--TF'):
        #     TF = arg

        else:
            print('Unknown option', opt)
            return 1

    #Check if required options are given
    if indir is None or outdir is None:
        print('Missing argument!')
        print(__doc__)
        return 1
    else:
        print('indir: ', indir)
        # print('outdir: ', outdir)
    if TF is not None and cell is not None:
        print('TF:   ', TF)
        print('cell:   ', cell)
    # elif TF is not None and cell is None:
    #     print('TF:   ', TF)
    # elif TF is None and cell is not None:
    #     print('cell:   ', cell)
    else:
        print('all cell and TF combos')

    # Run panda
    print('Starting generic all X all analysis ...')
    
    if os.path.exists(indir+'/tmp'):
        shutil.rmtree(indir+'/tmp')
    cambPredScore(indir)
    # RR='all'
    # DD='all'
    # camb_plotPredScore(outdir,'auroc')
    # camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    print('Finished generic analysis all X all')

#     print('Starting regional analysis ...')
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'N_Shore',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'S_Shore',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'OpenSea',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'N_Shelf',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
# ##STILL NEED TO RUN
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'S_Shelf',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)

#     RR,DD=predScore(indir,outdir,'mean',cell,TF,'Island',None)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     print('Finished regional analysis ...')
    
#     print('Starting depth analysis ...')
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,5)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,10)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,15)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)
#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,20)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)

#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,30)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)

#     RR,DD=predScore(indir,outdir,'mean',cell,TF,None,50)
#     plot_predScore(outdir,outdir,'mean','auroc',RR,DD)
#     plot_allPredScore(outdir,'mean','auroc',RR,DD)
#     print('Finished '+RR+' '+DD)


#     print('Finished depth analysis ...')

    # print('Start predRegion ...')
    # predRegion(indir,outdir,cell,TF)
    # print('Start buffer_distr_comp ...')
    # buffer_distr_comp('mean',indir)
    # print('Start plot_predScore ...')
    
    
    
    
    
    # plot_predScore(outdir,outdir,'mean','aupr')
    # print('Start plot_allPredScore ...')
    
    
    
    
    

    # plot_allPredScore(outdir,'mean','aupr')
    # # print('All done!')
    
    # print('Start predScore ...')
    # predScore(indir,outdir,'median',cell,TF)
    # # print('Start buffer_distr_comp ...')
    # # buffer_distr_comp('median',indir)
    # print('Start plot_predScore ...')
    # plot_predScore(outdir,outdir,'median','auroc')
    # plot_predScore(outdir,outdir,'median','aupr')
    # print('Start plot_allPredScore ...')
    # plot_allPredScore(outdir,'median','auroc')
    # plot_allPredScore(outdir,'median','aupr')
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
