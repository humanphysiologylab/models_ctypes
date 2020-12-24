include ../../liblsoda/Make.inc

SRCS =  \
../../liblsoda/src/cfode.c \
../../liblsoda/src/cfode_static.c \
../../liblsoda/src/common.c \
../../liblsoda/src/corfailure.c \
../../liblsoda/src/correction.c \
../../liblsoda/src/daxpy.c \
../../liblsoda/src/ddot.c \
../../liblsoda/src/dgefa.c \
../../liblsoda/src/dgesl.c \
../../liblsoda/src/dscal.c \
../../liblsoda/src/fnorm.c \
../../liblsoda/src/idamax.c \
../../liblsoda/src/intdy.c \
../../liblsoda/src/lsoda.c \
../../liblsoda/src/methodswitch.c \
../../liblsoda/src/orderswitch.c \
../../liblsoda/src/prja.c \
../../liblsoda/src/scaleh.c \
../../liblsoda/src/solsy.c \
../../liblsoda/src/stoda.c \
../../liblsoda/src/vmnorm.c \
../../liblsoda/src/strdup_printf.c \
./maleckar.c \
./tor.c\
./run.c

OBJS = $(SRCS:%.c=%.o)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -DNDEBUG -c $< -o $@

release debug: maleckar.$(SHLIB_EXT)

maleckar.$(SHLIB_EXT): $(OBJS)
	rm -rf $@
	$(CC) $(LDFLAGS) $(LDFLAGS-add) -shared -o $@ $^

clean:
	rm -f ./*.o ../../liblsoda/src/*.o *.$(SHLIB_EXT)
