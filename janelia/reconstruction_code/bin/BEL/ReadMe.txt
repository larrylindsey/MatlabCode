BEL v1.5.  
This software package is provided for research purpose only.
It a slightly improved version of the one reported in:
Piotr Dollar, Zhuowen Tu, and Serge Belongie, "Supervised Learning of Edges and Object Boundaries", CVPR, June, 2006.
But you may use the above paper for reference.

This package is for Windows ONLY, and it includes:

(1) An executable program BEL.exe
(2) IPL dlls (developed by the Intel lab)
(3) The Berkeley training images (not complete) and testing images.
David Martin, Charless Fowlkes and Jitendra Malik,
"Learning to Detect Natural Image Boundaries Using Local Brightness, Color and Texture Cues", PAMI,
26(5), 530-549, May 2004.
(It may require some additional dlls related to Visual Studio) 
(4) Two classifiers have been trained on the Berkeley training images and
they are (name clf) saved in .\Data\bsd\edge and .\Data\bsd\edge_onCanny

You simply execute "BEL yoursetting.txt" and the yoursetting.txt file specifies where are the training/testing images and how to tune some the parameters for training/testing.
Run "BEL bsd_fast_edge.txt" is a fast version (around 12 seconds on an image of 481x321) with slightly
degraded performance, but still within the same range of BEL.

To train, simply set the training=1 in "setting.txt".
To test, set training=0

All the images are saved in a directory as names "I00000.tif" with the corresponding annotation images (same size) named "I00000_label.tif". If you want to work with color images, simply name the original images as "I00000C.tif". These ".tif" images are saved in an un-compressed format.

In addition, you can also give a mask image, "I00000_premask.tif", if you wish to only classify those pixels interested. But you should make sure that your taining and testing stages are consistent. The parameter "mask" specifies whether or not a mask image is used. If "mask=1", then the algorithm automatically computes an edge map by Canny detector at a low scale (sigma=1.5), and trains/tests those edge pixels. If "mask=2", the program trains/tests on those specified by "_premask.tif" defined by users.

The most important parameters to tune in the training are the maximum depth of the tree, size of image patch, how many cascade levels enforced, and the number of weak classifiers for each node. The final classifier "clf" is saved at the directory specified by "path results=". You can also specify the stopping criterion on a leaf node: the threshold to allow the portion of pixels miss-classified. Increasing the number of tree levels and lowering the thresholds may result in longer training time and over-fitting, depending upon the complexity of your taining data.

Some system functions will be called to delete some intermediate files, which may trigger some warnings. You may simply ignore them.

You can play with it and let us know if you have any questions at
pdollar@cs.ucsd.edu and zhuowen.tu@gmail.com.
