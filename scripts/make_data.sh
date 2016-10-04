#!/bin/bash

PARALLEL=
# select parallelization method, if any
#PARALLEL=moreutils
PARALLEL=gnuparallel
# number of parallel jobs, e.g., number of CPU cores
PARJOBS=6

NUMRUNS=24

make_fIs() {
    NUMBRANCHES=$1
    SIMSCRIPT=scripts/make_fIs.py
    # Measure fI curves
    case $PARALLEL in 
        moreutils) 
            parallel -j${PARJOBS} $SIMSCRIPT $NUMBRANCHES -- $(seq $NUMRUNS)
        ;;
        gnuparallel)
            parallel -j${PARJOBS} $SIMSCRIPT $NUMBRANCHES ::: $(seq $NUMRUNS)
        ;;
        *)
            for RUN in $(seq $NUMRUNS); do $SIMSCRIPT $NUMBRANCHES $RUN; done
        ;;
    esac

    # post-process (collect) fI recordings
    scripts/collect_fIs.py
}

run_sims() {
    SIMSCRIPT=$1
    case $PARALLEL in 
        moreutils) 
            parallel -j${PARJOBS} $SIMSCRIPT -- $(seq $NUMRUNS)
        ;;
        gnuparallel)
            parallel -j${PARJOBS} $SIMSCRIPT ::: $(seq $NUMRUNS)
        ;;
        *)
            for RUN in $(seq $NUMRUNS); do $SIMSCRIPT $RUN; done
        ;;
    esac
}

figure2() {
    # ~ 1 min
    time make_fIs 2
    # ~ 10 min
    time run_sims scripts/figure2.py
}

figureS5() {
    # ~ 7 min
    time make_fIs 8
    # ~ 43 min
    time run_sims scripts/figureS5.py
}

figure2
figureS5
