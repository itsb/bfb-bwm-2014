#!/bin/bash

case $1 in
    # cleanup environment to start fresh
    --cleanup)
        rm -f figs/Figure*.png
        rm -f data/*.h5 data/*.pkl
        rm -f -r mods/{i686,x86_64,powerpc,umac}
        exit
    ;;
    *)
esac

makeoutputdirs() {
mkdir -p data figs
}

makedll() {
cd mods
nrnivmodl
cd ..
}

generate_data() {
export NRN_NMODL_PATH=${PWD}/mods
export HOC_LIBRARY_PATH=${PWD}/hoc
scripts/make_data.sh
}

postprocess_data() {
py/postprocess.py
}

generate_plots() {
py/make_figures.py
}

figure4() {
export NRN_NMODL_PATH=${PWD}/mods
export HOC_LIBRARY_PATH=${PWD}/hoc
scripts/figure4.py
}

makeoutputdirs
makedll

# for figure 4
# run time is ~2 seconds
figure4

# for figures 2 and S5
# run time is ~11 min w/6x parallelization on i7-3960X for Figure 2
# run time is ~50 min w/6x parallelization on i7-3960X for Figure S5
# see PARALLEL and PARJOBS options in scripts/make_data.sh
generate_data
postprocess_data
generate_plots
