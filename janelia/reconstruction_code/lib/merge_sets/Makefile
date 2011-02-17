all: libmerge_sets.a testbed testbed_h

modules= merge_sets.o

%.o: %.cpp
	gcc -O3 -c -fPIC $*.cpp -I./

libmerge_sets.a: $(modules)
	ar -r libmerge_sets.a merge_sets.o

testbed: testbed_merge_sets.cpp $(modules)
	gcc -o testbed_merge_sets testbed_merge_sets.cpp -fPIC -I./ -L./ -lmerge_sets -lstdc++

testbed_h: testbed_merge_sets_h.cpp
	gcc -o testbed_merge_sets_h testbed_merge_sets_h.cpp -fPIC -I./ -L./ -lstdc++

clean:
	rm $(modules) testbed_merge_sets testbed_merge_sets_h

