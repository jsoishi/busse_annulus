#!/usr/bin/env bash

name=$1
if [ -d $name ]; then
    echo "$1: directory exists."
    exit
fi
mkdir $name
cd $name
ln -s ../simulation.py .
cp ../parameters.py .
cp ../run_dss_busse.sh .
ln -s ../plot_profiles.py .
ln -s ../plot_snapshots.py .
ln -s ../plot_scalars.py .
ln -s ../merge.sh .
