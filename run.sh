#!/bin/bash

# Check if at least two parameters are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: ./run.sh <data_file> <number_of_frames> [-gpu]"
    exit 1
fi

mkdir -p outputs

# Use the parameters in the script
if [ "$3" == "-gpu" ]; then
    ./build/lbm "data/$1.txt" $2 -gpu
else
    ./build/lbm "data/$1.txt" $2
fi

source env/bin/activate
python scripts/plotting2D.py outputs/velocity_out.txt