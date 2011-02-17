// An example program to write an npy-format binary file. 
// Requires npy.h and libnpy.a from Bill McLean's libnpy library.

#include "npy.h"
#include <string>

using std::string;

int main(int argc, char** argv) {
    double a[2][4] = { { 1, 2, 3, 4 }, {5, 6, 7, 8 } };
    int shape[2] = {2, 4};
    int fortran_order = 0;
    string filename("ca.npy");
    if( argc == 2 ) {
        filename = argv[1];
    }

    npy_save_double(const_cast<char*>(filename.c_str()),
                                            fortran_order, 2, shape, &a[0][0]);
    return 0;
}
