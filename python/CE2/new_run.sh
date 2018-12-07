#!/usr/bin/env bash

name=$1
mkdir $name
cd $name
ln -s ../simulation.py .
cp ../parameters.py .
cp ../run_dss_busse.sh .
ln -s ../plot_profiles.py .
ln -s ../plot_snapshots.py .
ln -s ../plot_scalars.py .
