CC = gcc
INCLUDE_DIR = -I/usr/local/include -I./
CFLAGS = $(INCLUDE_DIR) -g -fPIC -O3
CPPFLAGS = $(INCLUDE_DIR) -g -fPIC -O3
LIB_DIR = -L/usr/local/lib -L./
LIB = -lgut
LDFLAGS = $(LIB_DIR) $(LIB) -pg

all: libgut.a testbed 

libgut_modules = digraph.o

testbed: testbed.o libgut.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $@

libgut.a: $(libgut_modules)
	ar -r libgut.a $(libgut_modules)
	rm $(libgut_modules)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
