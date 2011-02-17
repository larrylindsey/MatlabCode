#!/usr/bin/env python
#
# ilastik2png -- Convert Ilastik HDF5 to sequence of PNG
# Copyright 2010 HHMI.  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following disclaimer
# in the documentation and/or other materials provided with the
# distribution.
#     * Neither the name of HHMI nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author: katzw@janelia.hhmi.org (Bill Katz)
#  Written as part of the FlyEM Project at Janelia Farm Research Center.

import os
import sys
import errno
import getopt
import h5py
import Image
import numpy

HELP_MESSAGE = """
Usage:
    ilastik2png <hdf5 file> <directory-to-create>

options:
    -g, --group:  Follow with group name.  (assumes "prediction" if this option is not supplied)
    -c, --class:  Follow with class number.  (0, 1, etc.)

"""

class Error(Exception):
    """Base-class for exceptions in this module."""

class UsageError(Error):
    def __init__(self, msg):
        self.msg = msg

def main(argv):
    try:
        try:
            opts, args = getopt.gnu_getopt(argv, 'h:', ["help"])
        except getopt.error, msg:
            raise UsageError(msg)

        # option processing
        group_name = 'prediction'
        class_number = 0
        for option, value in opts:
            print "Looking at option:", str(option), str(value)
            if option in ("-h", "--help"):
                raise UsageError(HELP_MESSAGE)
            if option in ("-g", "--group"):
                group_name = value
            if option in ("-c", "--class"):
                class_number = int(value)

        if len(args) < 3 or len(args) > 3:
            raise UsageError(HELP_MESSAGE)
        else:
            hdf5_filename = args[1]
            png_directory = args[2]
            print "\nConverting %s to sequence of PNG files in %s... " % (hdf5_filename, png_directory)
            # Make directory if it doesn't exist
            try:
                os.makedirs(png_directory)
            except OSError as exc:
                if exc.errno == errno.EEXIST:
                    pass
                else:
                    raise
            # Store png files
            f = h5py.File(hdf5_filename, 'r')
            dset = f[group_name]
            for z in xrange(0, dset.shape[3]):
                float32_slice = dset[0, 0:dset.shape[1], 0:dset.shape[2], z, class_number]
                float32_slice *= 255.0 / float32_slice.max()  # Normalize from 0 to 255
                uint8_slice = float32_slice.astype(numpy.uint8)
                img = Image.fromarray(uint8_slice)
                filename = '/'.join([png_directory, "%s-%d.png" % (group_name, z)])
                img.save(filename, "PNG")

    except UsageError, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
