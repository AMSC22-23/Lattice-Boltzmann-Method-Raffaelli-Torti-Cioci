#!/bin/bash

# Check if two parameters are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./run.sh <data_file> <number>"
    exit 1
fi

mkdir -p outputs

# Use the parameters in the script
./build/lbm "data/$1.txt" $2
source env/bin/activate
python scripts/plotting2D.py outputs/velocity_out.txt