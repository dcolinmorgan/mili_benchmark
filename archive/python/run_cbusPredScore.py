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
    cambPredScore(indir,None)
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    print('Finished generic analysis all X all')
    
    cambPredScore(indir,"TSS1500")
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    
    cambPredScore(indir,"TSS200")
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    
    cambPredScore(indir,"5'UTR")
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    
    cambPredScore(indir,"3'UTR")
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    
    cambPredScore(indir,'1stExon')
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
    
    cambPredScore(indir,'ExonBnd')
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')
  
    cambPredScore(indir,'ExonBnd')
    camb_plotPredScore(outdir,'auroc')
    camb_allPredScore('../../pc/redmo/data/MotifPipeline/**/',outdir,'auroc')

    print('Starting regional analysis ...')
    print('Finished regional analysis ...')
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
