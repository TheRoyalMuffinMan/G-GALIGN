# G-GALIGN
GPU and CPU power global sequence alignment

## Single-Threaded Run Command:
From root directory and after setting `LD_LIBRARY_PATH`, run the following:
``` ./galign_benchmark -q test/pfizer_mrna.fna -r test/sars_spike_protein.fna -o output.txt -g -2 -p -1 -m 1```