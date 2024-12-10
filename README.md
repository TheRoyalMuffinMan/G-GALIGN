# G-GALIGN
GPU and CPU powered global sequence alignment \
Authors: Andrew Hoyle & Zachary Benson

## Build
```bash
make all
```
The above line will compile all the files (assuming pthreads and nvcc is installed) using the Makefile
After compilation, manually linking must be done.
Reference the configuration file (`configurations.sh`) and run the following command to handle linking:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./build
```

## Single-Threaded Run Command:
From root directory, after compiling and setting `LD_LIBRARY_PATH`, run the following:
```bash
./galign_benchmark -q test/pfizer_mrna.fna -r test/sars_spike_protein.fna -o output.txt -g -2 -p -1 -m 1
```

## Multi-Threaded Run Command:
From root directory, after compiling and setting `LD_LIBRARY_PATH`, run the following:
```bash
./multi_galign_benchmark -q test/pfizer_mrna.fna -r test/sars_spike_protein.fna -o output.txt -g -2 -p -1 -m 1
```

## GPU Run Command:
From root directory, after compiling and setting `LD_LIBRARY_PATH`, run the following:
```bash
./gpu_galign_benchmark -q test/pfizer_mrna.fna -r test/sars_spike_protein.fna -o output.txt -g -2 -p -1 -m 1
```