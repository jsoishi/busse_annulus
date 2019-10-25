#!/usr/bin/env bash

old=$1
new=$2
if [ -d $new ]; then
    echo "$2: directory exists."
    exit
fi
mkdir $new
cd $new
ln -s ../simulation.py .
cp ../$old/parameters.py .
cp ../$old/run_dss_busse.sh .
ln -s ../plot_profiles.py .
ln -s ../plot_snapshots.py .
ln -s ../plot_scalars.py .
ln -s ../merge.sh .
