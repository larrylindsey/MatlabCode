CC = gcc
INCLUDE_DIR = -I/usr/local/include -I./
CFLAGS = $(INCLUDE_DIR) -g -fPIC -O3 -Wno-deprecated
CPPFLAGS = $(INCLUDE_DIR) -g -fPIC -O3 -Wno-deprecated
LIB_DIR = -L/usr/local/lib -L./
LIB = -liput
LDFLAGS = $(LIB_DIR) $(LIB) -pg

all: libiput.a testbed testbed_array raw2ascii ascii2raw

libiput = file_input_output.o split_param_string.o

testbed: testbed.o libiput.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $@

raw2ascii: raw2ascii.o libiput.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $@

ascii2raw: ascii2raw.o libiput.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $@

testbed_array: testbed_array.o libiput.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $@

libiput.a: $(libiput)
	ar -r libiput.a $(libiput)
	rm $(libiput)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o