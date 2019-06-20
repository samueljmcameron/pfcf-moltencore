ATLAS_HOME = /usr/local/atlas/lib
CFLAGS = -L$(ATLAS_HOME)/lib -I. -lcblas -latlas -lm

spline: spline.c qromb.c trapzd.c polint.c nrutil.c
	gcc -o spline spline.c qromb.c trapzd.c polint.c nrutil.c $(CFLAGS)
