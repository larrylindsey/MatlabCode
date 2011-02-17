#!/bin/env python
"""

short script to examine the annotation.pickle file that Raveler 
produces on export

usage:

annotation-viewer.py picklefile


djo, 3/16/10
snv, 3/17/10

"""

# ------------------------- imports -------------------------

import cPickle as pickle
import os
import sys


# ------------------------- script starts here -------------------------
if len(sys.argv) < 2:
    print "usage: annotation_pickle_to_ascii.py picklefile"
    sys.exit(2)

filename = sys.argv[1]
if not os.path.exists(filename):
    print "file %s seems not to exist" % filename
    sys.exit(1)

f = open(filename, 'rb')
data = pickle.load(f)
f.close()


# at this point, "data" contains one big Python dictionary with
#   a hierarchy of subdictionaries

# data["point"] = dictionary of point annotations; keys = (x, y, z) points
# data["body"] = dictionary of body annotations; keys = body ids

# the dictionaries get messy...and they don't have all the info
#   you need; here are examples of what you can do:

# output corrected body IDs
bodydict = data["body"]
for bodyid in bodydict:
    if "corrected" in bodydict[bodyid]:
        if bodydict[bodyid]["corrected"] == "corrected":
            print 'BODY %d "corrected"' %(bodyid)
    if "comment" in bodydict[bodyid]:
        print 'BODY %d "COMMENT:\n__ID__(%d)\n%s"' %(bodyid, bodyid, bodydict[bodyid]["comment"])

# output bookmarks:
pointdict = data["point"]
for x, y, z in pointdict:
    annotation = pointdict[x, y, z]
    if annotation["kind"] == "bookmark":
        print 'VOXEL %d %d %d "%s"' % (z, y, x, annotation["value"])

# output T-bars
# NOTE: at this time, T-bars are still somewhat incomplete; the T-bar
#   location is exported, but the potential post-synaptic partners are
#   NOT (they are only included in the separate text-file export); there
#   is a placeholder for an "edited" list of post-synaptic partners, but
#   that is not being used yet (shifting priorities, etc., etc.)

# output T-bar locations:
pointdict = data["point"]
for x, y, z in pointdict:
    annotation = pointdict[x, y, z]
    if annotation["kind"] == "T-bar":
        print 'VOXEL %d %d %d "T-bar: %s"' % (z, y, x, annotation["value"])


