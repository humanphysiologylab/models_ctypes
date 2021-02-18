gcc  -Wall -std=c99 -O3 -ffast-math -march=native -DCFODE_STATIC -fPIC -DNDEBUG -c model.c -o model.o
gcc -shared -o model.so model.o
