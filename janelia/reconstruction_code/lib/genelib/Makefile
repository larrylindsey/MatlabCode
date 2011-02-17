all: libmyers.a

modules = cdf.o fct_min.o image_filters.o lu_matrix.o xfig.o \
          contour_lib.o fct_root.o image_lib.o utilities.o   \
          draw_lib.o hash.o level_set.o water_shed.o fft_my.o tiff_io.o tiff_image.o

%.o: %.p
	awk -f manager.awk $< > $*.c
	gcc -g -c -fPIC $*.c

libmyers.a: $(modules)
#libtool --mode=link gcc -static -o libmyers.a *.o
	ar -r libmyers.a *.o
clean:
	rm $(modules) *.c

package:
	tar -zcf myers.pack.tar.gz README *.h *.p *.pdf manager.awk Makefile
