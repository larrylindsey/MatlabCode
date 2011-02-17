CFLAGS = -Wall -fpic -I../include -L../lib
FC = gfortran
FFLAGS = -Wall -fpic -I../include -L../lib
LIBNPY = -lnpy
CDATA_FILES =  a.npy  b.npy  c.npy
FDATA_FILES = fa.npy fb.npy fc.npy fd.npy

vpath %.h       ../include
vpath %.mod     ../include
vpath test%.c   ../test
vpath test%.f95 ../test
vpath %.a %.so  ../lib

libnpy.a: npy.o fnpy.o
	ar rc $@ $^
	mv $@ ../lib

libnpy.so.0.2: npy.o fnpy.o
	gfortran -o $@ -shared -Wl,-soname=libnpy.so.0 $^ 
	mv $@ ../lib

check: check_npy.py $(CDATA_FILES) $(FDATA_FILES)
	./$<

$(CDATA_FILES): test_npy
	./$<

$(FDATA_FILES): test_fnpy
	./$<

test_metadata: test_metadata.o npy.o
	gcc $(CFLAGS) -o $@ $^

test_npy: test_npy.o libnpy.a
	gcc $(CFLAGS) -o $@ $^ $(LIBNPY)

test_fnpy: test_fnpy.o libnpy.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBNPY)

npy.o: npy.c npy.h
	gcc $(CFLAGS) -c $^

test_npy.o: npy.h test_npy.c
test_metadata.o: npy.h test_metadata.c libnpy.a

fnpy.o fnpy.mod: fnpy.f95
	$(FC) $(FFLAGS) -c $^
	mv fnpy.mod ../include

test_fnpy.o: test_fnpy.f95 fnpy.mod libnpy.a
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f test_npy test_metadata test_fnpy *.npy mtd *.o 

reallyclean: clean
	rm -f ../lib/libnpy.* ../include/fnpy.mod ../include/*.h.gch
