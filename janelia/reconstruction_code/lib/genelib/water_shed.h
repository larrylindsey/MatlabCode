/*****************************************************************************************\
*                                                                                         *
*  Watershed Partitioning                                                                 *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  June 2007                                                                     *
*  Mod   :  Aug 2008 -- Added Waterhshed Graph and Tree concepts                          *
*                                                                                         *
\*****************************************************************************************/

#ifndef _WATERSHED

#define _WATERSHED

#include "image_lib.h"

#ifdef __cplusplus
extern "C"{
#endif
  
#define WSHED    0
#define ONQUEUE -1
#define MASK    -2
#define INIT    -3
#define MINIM   -4
#define OUTBOUND -5

typedef struct
  { int seed;     // a pixel in the catchment basin
    int first;    // index of first edge in "edges" array
    int depth;    // minimum pixel value in the CB corresponding to this vertex;
    int size;     // # of pixels in the basin
  } Basin;

typedef struct
  { int head;     // index of CB on the other side of the damn
    int height;   // least barrier height separating the minimima of the two CBs of the edge
  } Dam;

typedef struct
  { int    iscon4;    // 4-connected (vs. 8-connected) regions?
    int    nlabels;   // The number of distinct catchment basin (CB) labels (1..nlabels)
    Image *labels;    // A GREY image whose values are the CB labels
    int    nbasins;   // The number of CBs in the partition
    Basin *verts;     // An array [0..nbasins-1] of CB records
    Dam   *edges;     // The edges for CB v are v->edges[j].head of weight edges[j].height,
                      //   for all j in [verts[v].first,verts[v+1].first-1]
  } Watershed_2D;

//  Watershed partition (we use a modification of Vincent & Soille, IEEE. Trans. on Pat.
//    Analysis & Mach. Learning, 13, 6 (1991)) an image using either 4- (con4 on) or 8-
//    (con4 off) connectivity.  The image must be either GREY or GREY16.  In the "labels"
//    image that is returned the set of pixels in a given catchment basin (CB) have the
//    same GREY value and values are assigned to CBs in such a way that two CBs with the
//    same label are never adjacent.  Typically no more than a handful of labels are used
//    so encoding them in the GREY image "labels" above is safe (no attempt is made to
//    optimize, we just greedily color the 10's of thousands of CBs that can result, but
//    typically never need more than 8 or 9 colors/labels to color the watershed partition)
//    
//  In addition, we build an adjacency list representation of what I call the dam graph.
//    Each vertex corresponds to a CB and there is an undirected edge between two CB's
//    v and w if they share a boundary in watershed partition.  The weight of the edge
//    corresponds to the barrier height between the CB's, namely it is the elevation from
//    the higher of the two CB minimums to the saddle height between them.  The encoding
//    of the graph is explained by the comments for the fields of the Basin and Dam structs.

Watershed_2D *Build_2D_Watershed(Image *image, int con4);

Watershed_2D *Copy_Watershed_2D(Watershed_2D *part);
void          Pack_Watershed_2D(Watershed_2D *part);
void          Free_Watershed_2D(Watershed_2D *part);
void          Kill_Watershed_2D(Watershed_2D *part);
void          Reset_Watershed_2D();
int           Watershed_2D_Usage();

// Collapse_2D_Watershed progressively merges watersheds in order of barrier height as
//   directed by "handler".  The inputs are a watershed and a user-suppled vector "base"
//   of shed->nbasins elements of size "size" in bytes.  When a pair of watersheds is
//   considered for merging by "handler", pointers to the two "base" elements for those
//   watershed is given as well as the barrier height between the watersheds.  Handler
//   should return a non-zero value if the two elements are to be merged, and zero if they
//   are not.  In the event of a merge, the left argument "base" record will represent the
//   combined watershed and may appropriately be modified by the user.  A pointer to an array M
//   of shed->nbasins integers is returned where M[i] gives the number of the new basin that
//   basin i now belongs to.  The union of M[i] = 0.."newbasins"-1.

int *Collapse_2D_Watershed(Watershed_2D *shed, void *base, int size,
                           int (*handler)(void *,int,void *), int *newcbs);

// Merge_2D_Watershed merges the watersheds in shed according the to mapping in "labels" such
//   as produced by Collpase_2D_Watershed above.  Newbasins must be the number of basins in
//   the new watershed.  A new watershed object is returned in which all the watersheds of
//   the input watershed are merged.

Watershed_2D *Merge_2D_Watershed(Watershed_2D *shed, int *labels, int newbasins);

// Produce a colored image of the watershed shed of the original image source

Image *Color_Watersheds(Watershed_2D *shed, Image *source);

// 3D Watershed (untested !!)

typedef struct
  { int    nlabels;   // The number of distinct catchment basin (CB) labels (1..nlabels)
    Stack *labels;    // A GREY stack whose values are the CB labels
    int    nbasins;   // The number of CBs in the partition
    int   *seeds;     // A pixel in each CB from which one can flood fill or contour
    int    iscon6;    // 6-connected (vs. 26-connected) regions?
  } Watershed_3D;

Watershed_3D *Build_3D_Watershed(Stack *stack, int con6);
Stack *Build_3D_Watershed2(Stack *stack, int con6);
Stack *Build_3D_Watershed3(Stack *stack, int con6);

Watershed_3D *Copy_Watershed_3D(Watershed_3D *part);
void          Pack_Watershed_3D(Watershed_3D *part);
void          Free_Watershed_3D(Watershed_3D *part);
void          Kill_Watershed_3D(Watershed_3D *part);
void          Reset_Watershed_3D();
int           Watershed_3D_Usage();

Stack *Color_Watersheds_3D(Watershed_3D *shed, Stack *source);

Watershed_3D *Build_3D_Watershed_Masked(Stack *frame, int iscon6, Stack * mask, unsigned int mask_label);

Watershed_3D *Build_3D_Watershed_Masked_Seeded(Stack *frame, int image_stack_scale, int iscon6, unsigned short int ** superpixel_maps, int original_width, double ** superpixel_2_seg_maps, unsigned int * pmap, unsigned int mask_label, Image * seed_image, int seed_depth, Stack * label_new);

int *Build_3D_Watershed_PreMasked_Seeded(Stack *frame, int iscon6, Image * seed_image, int seed_depth);

#ifdef __cplusplus
}
#endif // extern "C"

#endif
