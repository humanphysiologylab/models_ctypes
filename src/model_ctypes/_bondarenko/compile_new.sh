rm ./bondarenko.o
rm ./bondarenko.so
gcc  -Wall -std=c99 -O3 -ffast-math -march=native -DCFODE_STATIC -fPIC -DNDEBUG -c bondarenko_control.c -o bondarenko.o
gcc -shared -o bondarenko.so bondarenko.o
