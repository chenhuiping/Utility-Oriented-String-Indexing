#!/bin/bash

cd ./topk

for makefile in Makefile*; do
    echo "Running make with $makefile"
    make -f "$makefile" || { echo "Error with $makefile"; exit 1; }
done

cd ../usi

for makefile in Makefile*; do
    echo "Running make with $makefile"
    make -f "$makefile" || { echo "Error with $makefile"; exit 1; }
done