all: liblevelsets.a

modules = levelsets_ChanVese.o

%.o: %.cpp
	g++ -g -c -fPIC -O3 -I./ $*.cpp

liblevelsets.a: $(modules)
#libtool --mode=link gcc -static -o libmyers.a *.o
	ar -r liblevelsets.a *.o
clean:
	rm $(modules)

package:
	tar -zcf levelsets.tar.gz *.h *.cpp Makefile

foo: foo.cpp
	g++ -g -c -fPIC -I./ foo.cpp
