#!/bin/bash

# ==============================================================
# === Declarations
# Make sure script exits if error encountered
set -o pipefail
set -o errexit

# Set defaults for options
RMPNG='false'
RMDAT='false'
FRAMERATE=3


print_help() {
    echo "Cleans, makes, and runs simulation. Then, it creates animations."
    echo
    echo "-h     : Prints this message."
    echo "-r     : (default=false) to clean directory of png files"
    echo "-f OPT : (default=3) sets the frames per second of output video" 
    echo "-c     : (default=false) to clean data directory of *.csv files"
}

while getopts o:f:rh flag
do
    case "${flag}" in
        c) RMDAT='true';;
        f) FRAMERATE="${OPTARG}";;
        r) RMPNG='true';;
        h) print_help
            exit 1;;
    esac
done

# ==============================================================
# === Main script
# Clean old data and figures
#echo $'\nCleaning old data and figures..\n'
#rm data/*; rm figs/*.png

# Make and run simulation.
make && echo $'\nStarting program..\n\n' && ./test.out

# Plot outputs and create mp4 files.
echo $'\n\nStarting plotting..'
(python pyscripts/timeplots.py)
(python pyscripts/freqplots2.py)

cd figs; pwd
bash ./makemp4.sh -f $FRAMERATE

if [[ $RMPNG = 'true' ]]
then
    echo $'Cleaning figs..'
    rm figs/output*.png
fi

if [[ $RMDAT = 'true' ]]
then
    echo $'Cleaning data..'
    rm data/*.csv
fi

echo $'\nExiting\n'
