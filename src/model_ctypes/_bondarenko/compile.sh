rm ./*.o
rm ./*.so
gcc  -Wall -std=c99 -O3 -ffast-math -march=native -DCFODE_STATIC -fPIC -DNDEBUG -c bondarenko_control.c -o bondarenko.o
# gcc  -Wall -std=c99 -O3 -ffast-math -march=native -DCFODE_STATIC -fPIC -DNDEBUG -c model.c -o model.o
gcc -shared -o bondarenko.so bondarenko.o # model.o
