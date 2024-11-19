export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./build
# This needs to be ran after building (ldd the compiled file to see if libgalign.so is linked)
# If this doesn't work, copy it directly to the terminal and run the command (after build or make all)