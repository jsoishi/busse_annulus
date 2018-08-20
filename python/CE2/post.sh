
python3 ../scripts/merge_procs.py data_scalars/
python3 ../scripts/merge_sets.py data_scalars.h5 data_scalars/*.h5
python3 ../scripts/plot_scalars.py data_scalars.h5


python3 ../scripts/merge_procs.py data_profiles/
python3 ../scripts/merge_sets.py data_profiles.h5 data_profiles/*.h5
python3 ../scripts/plot_profiles.py data_profiles.h5

python3 ../scripts/merge_procs.py data_snapshots/
python3 ../scripts/plot_snapshots.py data_snapshots/*.h5
