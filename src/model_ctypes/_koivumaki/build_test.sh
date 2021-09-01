make clean && make
cp koivumaki.so libko.so
gcc -Wall -L$PWD -o out test.c -lko -lm
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
./out
