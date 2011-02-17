/*****************************************************************************************\
*                                                                                         *
*  Watershed Partitioning                                                                 *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  June 2007                                                                     *
*  Mod   :  Aug 2008 -- Added Waterhshed Graph and Tree concepts                          *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "image_lib.h"
#include "draw_lib.h"
#include "water_shed.h"

#undef  DEBUG
#undef  SHOW_GRAPH
#define DRAW_WATERSHED
#undef  DEBUG_COLLAPSE

/**************************************************************************************\
*                                                                                      *
*  WATERSHED SPACE MANAGEMENT ROUTINES                                                 *
*                                                                                      *
\**************************************************************************************/

#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)

//  Workspace routines for building a component tree

static int *get_chord(int area, char *routine)
{ static int  Max_Area = 0;
  static int *Array = NULL;

  if (area < 0)
    { free(Array);
      Max_Area = 0;
      Array    = NULL;
    }
  else if (area > Max_Area)
    { Max_Area = area;
      Array    = (int *) Guarded_Realloc(Array,sizeof(int)*Max_Area,routine);
      bzero(Array,sizeof(int)*Max_Area);  
    }
  return (Array);
}

static int *get_label(int nbasins, char *routine)
{ static int  Max_Len = 0;
  static int *Array = NULL;

  if (nbasins < 0)
    { free(Array);
      Max_Len = 0;
      Array   = NULL;
    }
  else if (nbasins > Max_Len)
    { Max_Len = nbasins;
      Array   = (int *) Guarded_Realloc(Array,sizeof(int)*Max_Len,routine);
      bzero(Array,sizeof(int)*Max_Len);  
    }
  return (Array);
}

static inline int watershed_2d_vsize(Watershed_2D *part)
{ return (sizeof(Basin)*(part->nbasins+1)); }

static inline int watershed_2d_esize(Watershed_2D *part)
{ return (sizeof(Dam)*(part->verts[part->nbasins+1].first)); }


typedef struct __Watershed_2D
  { struct __Watershed_2D *next;
    int                    vsize;
    int                    esize;
    Watershed_2D           watershed_2d;
  } _Watershed_2D;

static _Watershed_2D *Free_Watershed_2D_List = NULL;
static int    Watershed_2D_Offset, Watershed_2D_Inuse;

static inline void allocate_watershed_2d_verts(Watershed_2D *watershed_2d, int vsize, char *routine)
{ _Watershed_2D *object  = (_Watershed_2D *) (((char *) watershed_2d) - Watershed_2D_Offset);
  if (object->vsize < vsize)
    { if (object->vsize == 0)
        object->watershed_2d.verts = NULL;
      object->watershed_2d.verts  = Guarded_Realloc(object->watershed_2d.verts,vsize,routine);
      object->vsize = vsize;
    }
}

static inline void allocate_watershed_2d_edges(Watershed_2D *watershed_2d, int esize, char *routine)
{ _Watershed_2D *object  = (_Watershed_2D *) (((char *) watershed_2d) - Watershed_2D_Offset);
  if (object->esize < esize)
    { if (object->esize == 0)
        object->watershed_2d.edges = NULL;
      object->watershed_2d.edges  = Guarded_Realloc(object->watershed_2d.edges,esize,routine);
      object->esize = esize;
    }
}

static inline Watershed_2D *new_watershed_2d(int vsize, int esize, char *routine)
{ _Watershed_2D *object;

  if (Free_Watershed_2D_List == NULL)
    { object = (_Watershed_2D *) Guarded_Malloc(sizeof(_Watershed_2D),routine);
      Watershed_2D_Offset = ((char *) &(object->watershed_2d)) - ((char *) object);
      object->vsize = 0;
      object->esize = 0;
    }
  else
    { object = Free_Watershed_2D_List;
      Free_Watershed_2D_List = object->next;
    }
  Watershed_2D_Inuse += 1;
  allocate_watershed_2d_verts(&(object->watershed_2d),vsize,routine);
  allocate_watershed_2d_edges(&(object->watershed_2d),esize,routine);
  object-> watershed_2d.labels = NULL;
  return (&(object->watershed_2d));
}

static inline Watershed_2D *copy_watershed_2d(Watershed_2D *watershed_2d)
{ Watershed_2D *copy = new_watershed_2d(watershed_2d_vsize(watershed_2d),watershed_2d_esize(watershed_2d),"Copy_Watershed_2D");
  Watershed_2D  temp = *copy;
  *copy = *watershed_2d;
  copy->verts = temp.verts;
  if (watershed_2d_vsize(watershed_2d) != 0)
    memcpy(copy->verts,watershed_2d->verts,watershed_2d_vsize(watershed_2d));
  copy->edges = temp.edges;
  if (watershed_2d_esize(watershed_2d) != 0)
    memcpy(copy->edges,watershed_2d->edges,watershed_2d_esize(watershed_2d));
  if (watershed_2d->labels != NULL)
    copy->labels = Copy_Image(temp.labels);
  return (copy);
}

Watershed_2D *Copy_Watershed_2D(Watershed_2D *watershed_2d)
{ return (copy_watershed_2d(watershed_2d)); }

static inline void pack_watershed_2d(Watershed_2D *watershed_2d)
{ _Watershed_2D *object  = (_Watershed_2D *) (((char *) watershed_2d) - Watershed_2D_Offset);
  if (object->vsize > watershed_2d_vsize(watershed_2d))
    { object->vsize = watershed_2d_vsize(watershed_2d);
      if (object->vsize != 0)
        object->watershed_2d.verts = Guarded_Realloc(object->watershed_2d.verts,
                                                     object->vsize,"Pack_Watershed_2D");
      else
        { free(object->watershed_2d.verts);
          object->vsize = 0;
        }
    }
  if (object->esize > watershed_2d_esize(watershed_2d))
    { object->esize = watershed_2d_esize(watershed_2d);
      if (object->esize != 0)
        object->watershed_2d.edges = Guarded_Realloc(object->watershed_2d.edges,
                                                     object->esize,"Pack_Watershed_2D");
      else
        { free(object->watershed_2d.edges);
          object->esize = 0;
        }
    }
  if (watershed_2d->labels != NULL)
    Pack_Image(watershed_2d->labels);
}

void Pack_Watershed_2D(Watershed_2D *watershed_2d)
{ pack_watershed_2d(watershed_2d); }

static inline void free_watershed_2d(Watershed_2D *watershed_2d)
{ _Watershed_2D *object  = (_Watershed_2D *) (((char *) watershed_2d) - Watershed_2D_Offset);
  object->next = Free_Watershed_2D_List;
  Free_Watershed_2D_List = object;
  if (watershed_2d->labels != NULL)
    Free_Image(watershed_2d->labels);
  Watershed_2D_Inuse -= 1;
}

void Free_Watershed_2D(Watershed_2D *watershed_2d)
{ free_watershed_2d(watershed_2d); }

static inline void kill_watershed_2d(Watershed_2D *watershed_2d)
{ _Watershed_2D *object  = (_Watershed_2D *) (((char *) watershed_2d) - Watershed_2D_Offset);
  if (watershed_2d->labels != NULL)
    Kill_Image(watershed_2d->labels);
  if (object->esize != 0)
    free(watershed_2d->edges);
  if (object->vsize != 0)
    free(watershed_2d->verts);
  free(((char *) watershed_2d) - Watershed_2D_Offset);
  Watershed_2D_Inuse -= 1;
}

void Kill_Watershed_2D(Watershed_2D *watershed_2d)
{ kill_watershed_2d(watershed_2d); }

static inline void reset_watershed_2d()
{ _Watershed_2D *object;
  while (Free_Watershed_2D_List != NULL)
    { object = Free_Watershed_2D_List;
      Free_Watershed_2D_List = object->next;
      kill_watershed_2d(&(object->watershed_2d));
      Watershed_2D_Inuse += 1;
    }
}

int Watershed_2D_Usage()
{ return (Watershed_2D_Inuse); }

void Reset_Watershed_2D()
{ reset_watershed_2d();
  get_chord(-1,NULL);
  get_label(-1,NULL);
}

/**************************************************************************************\
*                                                                                      *
*  BUILD A 2D WATERSHED DECOMPOSITION                                                  *
*                                                                                      *
\**************************************************************************************/

static int cwidth   = -1;
static int cheight  = -1;
static int carea;

static int chk_width;
static int chk_height;
static int chk_iscon4  = -1;

static int n_nbrs4;       // # of directions (4 or 8)
static int neighbor4[8];  // offset from current pixel of a step in direction 'i'

// return vector v s.t. v[i] is true iff a step in direction 'i' is OK

static inline int *boundary_pixels_2d(int p)
{ static int bound[8];
  int x, xn, xp;
  int y, yn, yp;

  y = p/cwidth;
  x = p%cwidth;    
  
  yn = (y > 0);
  xn = (x > 0);
  xp = (x < chk_width);
  yp = (y < chk_height);
  if (chk_iscon4)
    { bound[0] = yn;
      bound[1] = xn;
      bound[2] = xp;
      bound[3] = yp;
    }
  else
    { bound[0] = yn && xn;
      bound[1] = yn;
      bound[2] = yn && xp;
      bound[3] = xn;
      bound[4] = xp;
      bound[5] = yp && xn;
      bound[6] = yp;
      bound[7] = yp && xp;
    }

  return (bound);                                                                       
}

//  Establish all globals for neighborhood operations

static void setup_globals_2D(int width, int height, int iscon4)
{ if (cwidth != width || chk_iscon4 != iscon4)
    { cwidth = width;
      carea  = width*height;

      if (iscon4)
        { neighbor4[0] = -cwidth;
          neighbor4[1] = -1;
          neighbor4[2] =  1;
          neighbor4[3] =  cwidth;
          n_nbrs4 = 4;
        }
      else
        { int x, y;
    
          n_nbrs4 = 0;
          for (y = -cwidth; y <= cwidth; y += cwidth)
            for (x = -1; x <= 1; x += 1)
              if (! (x==0 && y == 0))
                neighbor4[n_nbrs4++] = y+x;
        }

      chk_width  = cwidth-1;
      chk_iscon4 = iscon4;
    }

  if (cheight != height)
    { cheight    = height;
      carea      = width*height;
      chk_height = cheight-1;
    }
}

Watershed_2D *Build_2D_Watershed(Image *frame, int iscon4)
{ Watershed_2D *result;

  int    nlabels;   // The number of labels (1..nlabels)
  int    nbasins;   // Number of watersheds
  Basin *verts;     // A list of records, one for each catchment basin
  Dam   *edges;     // A list of edges, encoded as described in the .h file

  int index[0x10001];
  int maxval;

  uint8  *value8;
  uint16 *value16;

  int   *chord;
  int   *out;

  if (frame->kind == COLOR || frame->kind == FLOAT32)
    { fprintf(stderr,"Build_2D_Watershed: Cannot apply to non-greyscale images\n");
      exit (1);
    }

  setup_globals_2D(frame->width,frame->height,iscon4);

  chord = get_chord(2*carea,"Build_2D_Watershed");
  out   = chord + carea;

  if (frame->kind == GREY16)
    { maxval  = 0x10000;
      value8  = NULL;
      value16 = (uint16 *) (frame->array);
    }
  else
    { maxval  = 0x100;
      value8  = frame->array;
      value16 = NULL;
    }

  // Establish lists of pixels with each value

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      for (p = 0; p < carea; p++)
        { v        = value8[p];
          chord[p] = index[v];
          index[v] = p;
        }
    else  // value16 != NULL
      for (p = 0; p < carea; p++)
        { v        = value16[p];
          chord[p] = index[v];
          index[v] = p;
        }
  }

#define PUSH(p)			\
{ if (qbot < 0)			\
    qtop = qbot = (p);		\
  else				\
    qtop = chord[qtop] = (p);	\
  chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

  { int v, p;

    for (p = 0; p < carea; p++)
      out[p] = MINIM;

    if (value8 != NULL)

      for (v = 0; v < maxval; v++)
        { int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

          qbot = -1;
          for (p = index[v]; p >= 0; p = t)
            { int *b, j;
  
              t = chord[p];
              b = boundary_pixels_2d(p);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j] && value8[p+neighbor4[j]] < v)
                  { out[p] = INIT;
                    PUSH(p)
                    break;
                  }
            }

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
          while (qbot >= 0)
            { int p, *b, j, q;
              int cdist;
  
              POP(p);
              b = boundary_pixels_2d(p);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j])
                  { q = p + neighbor4[j];
                    if (out[q] != INIT && value8[q] == v)
                      { out[q] = INIT;
                        PUSH(q);
                      }
                  }
            }
        }

    else // value16 != NULL)

      for (v = 0; v < maxval; v++)
        { int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

          qbot = -1;
          for (p = index[v]; p >= 0; p = t)
            { int *b, j;
  
              t = chord[p];
              b = boundary_pixels_2d(p);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j] && value16[p+neighbor4[j]] < v)
                  { out[p] = INIT;
                    PUSH(p)
                    break;
                  }
            }

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
          while (qbot >= 0)
            { int p, *b, j, q;
              int cdist;
  
              POP(p);
              b = boundary_pixels_2d(p);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j])
                  { q = p + neighbor4[j];
                    if (out[q] != INIT && value16[q] == v)
                      { out[q] = INIT;
                        PUSH(q);
                      }
                  }
            }
        }
  }

  // Build pixel value lists again but now with only INIT pixels in each list

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      { for (p = 0; p < carea; p++)
          if (out[p] == INIT)
            { v        = value8[p];
              chord[p] = index[v];
              index[v] = p;
            }
       }

     else // value16 != NULL
      { for (p = 0; p < carea; p++)
          if (out[p] == INIT)
            { v        = value16[p];
              chord[p] = index[v];
              index[v] = p;
            }
      }
  }

  // Find connected components of MINIM pixels and label with a unique catchment basin #.

  { int v, p;
    int qtop, qbot;

    nbasins = 0;

    qbot = -1;
    for (p = 0; p < carea; p++)
      if (out[p] == MINIM)
        { PUSH(p);
          out[p] = ++nbasins;
#ifdef DEBUG
          printf("Basin bottom: %d\n  (%d,%d)\n",nbasins,p%cwidth,p/cwidth);
#endif
          while (qbot >= 0)
            { int q, *b, j, r;

              POP(q);
              b = boundary_pixels_2d(q);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j])
                  { r = q+neighbor4[j];
                    if (out[r] == MINIM)
                      { PUSH(r);
                        out[r] = nbasins;
#ifdef DEBUG
                        printf("  (%d,%d)\n",r%cwidth,r/cwidth);
#endif
                      }
                  }
            }
        }
  }

  { int v, p;

    for (v = 0; v < maxval; v++)
      { int qbot, qtop, t;

#ifdef DEBUG
        printf("Level %d\n",v);
#endif

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

        qbot = -1;
        for (p = index[v]; p >= 0; p = t)
          { int *b, j;

            t      = chord[p];
            out[p] = MASK;
            b      = boundary_pixels_2d(p);
            for (j = 0; j < n_nbrs4; j++)
              if (b[j] && out[p+neighbor4[j]] >= 0)
                { out[p] = ONQUEUE;
                  PUSH(p);
                  break;
                }
          }

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

        while (qbot >= 0)
          { int p, *b, j, q, o;

            POP(p);
            b = boundary_pixels_2d(p);

            for (j = 0; j < n_nbrs4; j++)
              if (b[j])
                { q = p + neighbor4[j];
                  o = out[q];
                  if (o == MASK)
                    { out[q] = ONQUEUE;
                      PUSH(q);
                    }
                  else if (o != ONQUEUE && o != INIT)
                    { if (out[p] == ONQUEUE)
                        out[p] = o;
                      else if (o != out[p] && o != WSHED)
                        if (out[p] == WSHED)
                          out[p] = o;
                        else
                          out[p] = WSHED;
                    }
                }
#ifdef DEBUG
            printf("  Proc (%d,%d) %d\n",p%cwidth,p/cwidth,out[p]);
#endif
          }
      }
  }

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

  { int p, qbot, qtop;

    qbot = -1;
    for (p = 0; p < carea; p++)
      if (out[p] == WSHED)
        { int *b, j;

          b = boundary_pixels_2d(p);
          for (j = 0; j < n_nbrs4; j++)
            if (b[j] && out[p+neighbor4[j]] > 0)
              { PUSH(p);
                out[p] = ONQUEUE;
                break;
              }
        }

    while (qbot >= 0)
      { int *b, j, q, o, min, m;

        POP(p);
        b = boundary_pixels_2d(p);
#ifdef DEBUG
        printf("  Wshed Pixel: (%d,%d)\n",p%cwidth,p/cwidth);
#endif

        min = 0x10000;
        for (j = 0; j < n_nbrs4; j++)
          if (b[j])
            { q = p + neighbor4[j];
              o = out[q];
              if (value8 == NULL)
                m = value16[q];
              else
                m = value8[q];
              if (o > 0)
                { if (m < min)
                    { out[p] = o;
                      min = m;
                    }
                }
              else if (out[q] != ONQUEUE)
                { PUSH(q);
                  out[q] = ONQUEUE;
                }
            }
      }
  }

  //  Make the data structure now

  result = new_watershed_2d(sizeof(Basin)*(nbasins+1),0,"Build_2D_Watershed");
  result->nbasins = nbasins;
  result->iscon4  = iscon4;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typically the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, verts[b].seed is the new label for CB b,
  //    and nlabels is the number of unique labels used.
  //  Also, determine the minimum pixel of each CB b in verts[b].depth and the number of
  //    CB's adjacent to this CB in verts[b].first

#define APUSH(p)	\
{ chord[p] = atop;	\
  atop = (p);		\
}

#define APOP(p)  { atop = chord[(p)=atop]; }

  { int   qbot, qtop;
    int   atop;
    int   olabel, size;
    int   adjacent[256];
    int   p;

    verts   = result->verts;
    nlabels = 0;

    for (p = 0; p < 256; p++)
      adjacent[p] = 0;

    for (p = 0; p <= nbasins; p++)
      { verts[p].first = 0;
        verts[p].depth = -1;
      }

    qbot = -1;
    atop = -1;
    for (p = 0; p < carea; p++)

      if (out[p] > 0)         //  Flood fill each CB, marking it by flipping sign of out
        { PUSH(p);
          olabel = out[p];
          out[p] = 0;
          size   = 0;
#ifdef DEBUG
          printf("Push (%d,%d): %d\n",p%cwidth,p/cwidth,olabel);
#endif

          while (qbot >= 0)
            { int *b, j, q, r, o;

              POP(q);
              out[q] = -olabel;
              size  += 1;
      
              b = boundary_pixels_2d(q);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j])
                  { r = q + neighbor4[j];
                    o = out[r];
                    if (o == olabel)     
                      { PUSH(r);
                        out[r] = 0;
#ifdef DEBUG
                        printf("  Add (%d,%d)\n",r%cwidth,r/cwidth);
#endif
                      }
                    else
                      { o = -o;
                        if (o > 0 && o != olabel)   //  Pixel in adjacent, already filled CB
                          { if (verts[o].depth < 0)  //    do so and push on A-stack
                              { APUSH(r);
                                verts[o].depth = 1;
                                adjacent[verts[o].seed] = 1;
                              }
#ifdef DEBUG
                            printf("  Border (%d,%d) %d\n",r%cwidth,r/cwidth,o);
#endif
                          }
                      }
                  }
            }

          { int j, q;

            for (j = 1; j < 256; j++)    //  Determine the 1st available label for CB
              if ( ! adjacent[j])        //    This is O(edges) so still linear
                break;

#ifdef DEBUG
            printf("Basin %d: %d\n",olabel,j);
#endif
            verts[olabel].seed  = j;
            verts[olabel].size  = size;
            if (j > nlabels)
              { nlabels = j;
                if (j > 255)
                  { fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
                    fprintf(stderr,"         Returning NULL result\n");
                    fflush(stderr);
                    return (NULL);
                  }
              }

            while (atop >= 0)             //  Count edges between this CB and the already-filled
              { APOP(q);                  //    adjacent CBs
                q = -out[q];
                verts[q].first += 1;
                verts[olabel].first += 1;
                verts[q].depth = -1;
                adjacent[verts[q].seed] = 0;
#ifdef DEBUG
                printf("   Edge %d->%d (%d =? %d)\n",olabel,q,j,verts[q].seed);
#endif
              }
          }
        }
  }

  //  verts[p].first contains the # of edges for CB p.  Make it be the index one beyond
  //    the last edge of the CB and now that you know how many edges there are, allocate

  { int p;

#ifdef DEBUG
    for (p = 1; p <= nbasins; p++)
      { printf("  Basin %d: %d %d %d\n",p,verts[p].seed,verts[p].depth,verts[p].first);
        fflush(stdout);
      }
#endif

    for (p = 1; p <= nbasins; p++)
      verts[p].first += verts[p-1].first;
    verts[0].first = verts[nbasins].first;

    allocate_watershed_2d_edges(result,sizeof(Dam)*verts[nbasins].first,"Build_2D_Watershed");

    edges = result->edges;
  }

  //  Flood fill each CB again, this time inverting the sign of out back to a positive value.
  //    In this pass, the edges are filled in and barrier heights computed.

  { int   qbot, qtop;
    int   atop;
    int   olabel;
    int   p;

    qbot = -1;
    atop = -1;
    for (p = 0; p < carea; p++)

      if (out[p] < 0)         //  Flood fill each CB, marking it by flipping sign of out

        { PUSH(p);
          olabel = -out[p];
          out[p] = 0;

          while (qbot >= 0)
            { int *b, j, q, r, o, c;
              int m1, m2;

              POP(q);
              out[q] = olabel;
              if (value8 == NULL)
                m1 = value16[q];
              else
                m1 = value8[q];
      
              b = boundary_pixels_2d(q);
              for (j = 0; j < n_nbrs4; j++)
                if (b[j])
                  { r = q + neighbor4[j];
                    o = out[r];
                    if (o == -olabel)
                      { PUSH(r);
                        out[r] = 0;
                      }
                    else if (o > 0 && o != olabel)  //  Pixel in adjacent, already filled CB
                      { if (value8 == NULL)
                          m2 = value16[r];
                        else
                          m2 = value8[r];
                        if (verts[o].depth < 0)     //  Push on A-stack if first pixel in CB
                          { APUSH(r);
                            verts[o].depth = 0x10000;
                          }
                        if (m2 < m1)                //  Determine implied barrier height and
                          m2 = m1;
                        if (m2 < verts[o].depth)   //    figure into minimum over all pairs
                          verts[o].depth = m2;     //    of pixels between the CB's
                      }
                  }
            }

          while (atop >= 0)         //  Add the new edges and barrier heights for edges
            { int e, q;             //    between the current CB and already-filled adjacent ones

              APOP(q);
              q = out[q];

              e = verts[olabel].first-1;
              verts[olabel].first = e;
              edges[e].height = verts[q].depth;
              edges[e].head   = q;

              e = verts[q].first-1;
              verts[q].first = e;
              edges[e].height = verts[q].depth;
              edges[e].head   = olabel;

              verts[q].depth = -1;
            }
        }
  }

  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smallest pixel in each region.  Also shift CB numbering from 1... to 0...

  { int    p, b, nedges;
    int    d1, d2;
    uint8 *iarray;

    result->nlabels = nlabels;
    result->labels  = Make_Image(GREY,cwidth,cheight);
    iarray = result->labels->array;

    for (b = 1; b <= nbasins; b++)
      verts[b].depth = 0x10000;

    for (p = 0; p < carea; p++)
      { b = out[p];
        iarray[p] = verts[b].seed;
        if (value8 == NULL)
          { if (value16[p] < verts[b].depth)
              verts[b].depth = value16[p];
          }
        else
          { if (value8[p] < verts[b].depth)
              verts[b].depth = value8[p];
          }
      }

    for (p = 0; p < carea; p++)
      { if (verts[b = out[p]].seed > 0) 
          verts[b].seed = -(p+1); 
        out[p] = chord[p] = 0;
      }

    nedges = verts[0].first;
    for (b = 1; b <= nbasins; b++)
      { verts[b-1].seed  = -(verts[b].seed+1);
        verts[b-1].depth = verts[b].depth;
        verts[b-1].first = verts[b].first;
        verts[b-1].size  = verts[b].size;
      }
    verts[nbasins].first = nedges;

    for (p = 0; p < nedges; p++)
      edges[p].head -= 1;
  }

#ifdef SHOW_GRAPH

  { int i, j;

    for (i = 0; i < nbasins; i++)
      { j = verts[i].seed;
        printf("Basin %d(%d,%d,[%d,%d]) -> ",i,verts[i].depth,verts[i].size,j%cwidth,j/cwidth);
        for (j = verts[i].first; j < verts[i+1].first; j++)
          printf(" %d(%d)",edges[j].head,edges[j].height);
        printf("\n");
      }
  }

#endif

#ifdef DRAW_WATERSHED
  { int    i, j, p;
    uint8 *iarray;
    Image *view;

    iarray = result->labels->array;

    view = Color_Watersheds(result,frame);

    Write_Image("watershed.tif",view);

    Kill_Image(view);
  }
#endif

  return (result);
}


/****************************************************************************************\
*                                                                                        *
*  BUILD A WATERSHED TREE MODELING PROGRESSIVE MERGES OF CB's IN ORDER OF BARRIER HEIGHT *
*                                                                                        *
*****************************************************************************************/

static int ESORT(const void *l, const void *r)
{ Dam *x = (Dam *) l;
  Dam *y = (Dam *) r;
  return (x->height - y->height);
}

static int FSORT(const void *l, const void *r)
{ Dam *x = (Dam *) l;
  Dam *y = (Dam *) r;
  return (x[1].head - y[1].head);
}

static Basin *VERTS;

static int FIND(int v)
{ int w, z;

  w = v;
  while (VERTS[w].first >= 0)
    w = VERTS[w].first;
  while (VERTS[v].first >= 0)
    { z = VERTS[v].first;
      VERTS[v].first = w;
      v = z;
    }
  return (w);
}

static int UNION(int v, int w)       // not balanced, but should still be efficient
{ VERTS[w].first = v;
  return (v);
}

int *Collapse_2D_Watershed(Watershed_2D *shed, void *base, int size,
                           int (*handler)(void *,int,void *), int *newcbs)
{ Basin *verts;
  Dam   *edges;
  int    nbasins;
  int    nedges;
  int   *labels;

  int    v, w;
  int    e, f;

  verts   = shed->verts;
  edges   = shed->edges;
  nbasins = shed->nbasins;
  nedges  = verts[nbasins].first;

  if (nbasins == 1)
    return (NULL);

  labels = get_label(2*(nbasins+1),"Collapse_2D_Watershed");

  // Turn directed adjacency list into weight-ordered set of undirected edges
  //   Two consecutive dams encode an undirected edge, where the second "height"
  //   field is unused and the two "head" fields give the two ends of the edge.

  f = 0;
  for (v = 0; v < nbasins; v++)
    { e = verts[v].first;
      for (verts[v].first = f; e < verts[v+1].first; e++)
        if (edges[e].head < v)
          edges[f++] = edges[e];
    }
  verts[nbasins].first = f;

  f = nedges - 2;
  for (v = nbasins-1; v >= 0; v--)
    { for (e = verts[v+1].first-1; e >= verts[v].first; e--)
        { edges[f] = edges[e];
          edges[f+1].head = v;
          f -= 2;
        }
    }

  qsort(edges,nedges/2,2*sizeof(Dam),ESORT);

  // Do UNION/FIND merge to determine sets

  VERTS = verts;
  for (v = 0; v < nbasins; v++)
    verts[v].first = -1;

  f = nbasins-2;
  for (e = 0; e < nedges; e += 2)
    { v = FIND(edges[e].head);
      w = FIND(edges[e+1].head);
      if (v != w)
        { if (handler(base + v*size, edges[e].height, base + w*size))
            { verts[v].first += verts[w].first;
              v = UNION(v,w);
            }
        }
    }

  w = 1;
  for (v = 0; v < nbasins; v++)
    if (verts[v].first < 0)
      verts[v].first = -(w++);

  *newcbs = w-1;

  for (v = 0; v < nbasins; v++)
    labels[v] = -(verts[FIND(v)].first+1);

  // Turn the undirected edge representation back into the directed adjacency
  //   list representation so that it appears that nothing happened !

  qsort(edges,nedges/2,2*sizeof(Dam),FSORT);   // sort by second vertex (tail ultimately)

  for (v = 0; v < nbasins; v++)
    verts[v].first = 0;

  v = 0;
  for (f = 0; f < nedges; f += 2)
    { w = edges[f].head;
      if (verts[w].first >= 0)
        verts[w].first += 1;
      else
        verts[w].first -= 1;
      w = edges[f+1].head;
      if (verts[w].first >= 0)
        verts[w].first += 1;
      else
        verts[w].first -= 1;
      e = (f>>1);
      edges[e] = edges[f];
      if (w != v)
        { for (v++; v < w; v++)
            verts[v].first = -(verts[v].first+1);
          edges[e].head = -(edges[e].head+1);
        }
    }

  w = verts[0].first;
  if (w < 0)
    w = -(w+1);
  for (v = 1; v < nbasins; v++)
    if (verts[v].first < 0)
      { w += -(verts[v].first+1);
        verts[v].first = -(w+1);
      }
    else
      { w += verts[v].first;
        verts[v].first = w;
      }

  v = nbasins-1;
  for (e = nedges/2-1; e >= 0; e--)
    { f = verts[v].first-1; 
      verts[v].first = f;
      edges[f] = edges[e];
      if (edges[f].head < 0)
        { edges[f].head = -(edges[f].head+1);
          v -= 1;
          while (verts[v].first < 0)
            { verts[v].first = -(verts[v].first+1);
              v -= 1;
            }
        }
    }
  verts[nbasins].first = nedges;

  for (v = nbasins-1; v >= 0; v--)
    { for (e = verts[v].first; e < verts[v+1].first; e++)
        if (edges[e].head < v)
          { w = edges[e].head;
            f = verts[w].first-1;
            verts[w].first  = f;
            edges[f].head   = v;
            edges[f].height = edges[e].height;
          }
    }

#ifdef SHOW_GRAPH

  { int i, j;

    printf("\nRestored edges:\n");
    for (i = 0; i < nbasins; i++)
      { printf("Basin %d(%d) -> ",i,verts[i].depth);
        for (j = verts[i].first; j < verts[i+1].first; j++)
          printf(" %d(%d)",edges[j].head,edges[j].height);
        printf("\n");
      }
  }

#endif

  return (labels);
}


/**************************************************************************************\
*                                                                                      *
*  MERGE WATERSHEDS ALONG THE WATERSHED TREE TO PRODUCE A NEW WATERSHED OBJECT         *
*                                                                                      *
\**************************************************************************************/

Watershed_2D *Merge_2D_Watershed(Watershed_2D *shed, int *labels, int cbasins)
{ int    nbasins;
  Basin *verts;
  Dam   *edges;

  int   *ilabels;
  Basin *cverts;
  Dam   *cedges;
  int    clabels;

  Watershed_2D *result;

  nbasins = shed->nbasins;
  verts   = shed->verts;
  edges   = shed->edges;

  ilabels = labels + (nbasins+1);

  result = new_watershed_2d(sizeof(Basin)*(cbasins+1),0,"Collapse_2D_Watershed");
  result->nbasins = cbasins;
  result->iscon4  = shed->iscon4;
  cverts = result->verts;

  { int c, b;

    for (c = 0; c < cbasins; c++)
      cverts[c].size = 0;

    for (b = 0; b < nbasins; b++)
      cverts[labels[b]].size += 1;

    for (c = 1; c < cbasins; c++)
      cverts[c].size += cverts[c-1].size;

    for (b = 0; b < nbasins; b++)
      { c = cverts[labels[b]].size - 1;
        ilabels[c] = b;
        cverts[labels[b]].size = c;
      }

    labels[nbasins]  = cbasins;
    ilabels[nbasins] = nbasins;
  }
  
#ifdef DEBUG_COLLAPSE
  { int b;

    printf("\nNew basins = %d\n",cbasins);
    for (b = 0; b < nbasins; b++)
      printf("  %3d -> %d (%d)\n",b,labels[b],ilabels[b]);
  }
#endif

  { int b, c, d, e, f;
    int p, q;

    for (c = 0; c < cbasins; c++)
      { cverts[c].seed = -1;
        cverts[c].first = 0;
      }

    p = 0;
    for (c = 0; c < cbasins; c++)
      { for (q = p; labels[ilabels[q]] == c; q++) 
          { b = ilabels[q];
            for (e = verts[b].first; e < verts[b+1].first; e++)
              { d = labels[edges[e].head];
                if (d != c && cverts[d].seed < 0)
                  { cverts[d].seed = 1;
                    cverts[c].first += 1;
                  }
              }
          }

        for (q = p; labels[ilabels[q]] == c; q++) 
          { b = ilabels[q];
            for (e = verts[b].first; e < verts[b+1].first; e++)
              cverts[labels[edges[e].head]].seed = -1;
          }

        p = q;
      }

    for (c = 1; c < cbasins; c++)
      cverts[c].first += cverts[c-1].first;
    cverts[cbasins].first = cverts[cbasins-1].first;

    allocate_watershed_2d_edges(result,sizeof(Dam)*cverts[cbasins].first,"Collapse_2D_Watershed");
    cedges = result->edges;

    p = 0;
    for (c = 0; c < cbasins; c++)
      { for (q = p; labels[ilabels[q]] == c; q++) 
          { b = ilabels[q];
            for (e = verts[b].first; e < verts[b+1].first; e++)
              { d = labels[edges[e].head];
                if (d != c)
                  { if (cverts[d].seed < 0)
                      { f = cverts[c].first-1;
                        cverts[c].first  = cverts[d].seed = f;
                        cedges[f].head   = d;
                        cedges[f].height = edges[e].height;
                      }
                    else
                      { f = cverts[d].seed;
                        if (edges[e].height < cedges[f].height)
                          cedges[f].height = edges[e].height;
                      }
                  }
              }
          }

        for (q = p; labels[ilabels[q]] == c; q++) 
          { b = ilabels[q];
            for (e = verts[b].first; e < verts[b+1].first; e++)
              cverts[labels[edges[e].head]].seed = -1;
          }

        p = q;
      }
  }

  { int c, b;
    int area;

    area = shed->labels->width * shed->labels->height;

    for (c = 0; c < cbasins; c++)
      { cverts[c].seed  = area;
        cverts[c].size  = 0;
      }

    for (b = 0; b < nbasins; b++)
      { c = labels[b];
        if (verts[b].seed < cverts[c].seed)
          cverts[c].seed = verts[b].seed;
        cverts[c].size += verts[b].size;
      }
  }

  { int j, c, e, d;
    int colored[256];
    Paint_Brush _brush, *brush = &_brush;

    for (j = 1; j < 256; j++)
      colored[j] = 0;

    clabels = 0;
    for (c = 0; c < cbasins; c++)
      { for (e = cverts[c].first; e < cverts[c+1].first; e++)
          { d = cedges[e].head;
            if (d < c)
              colored[cverts[d].depth] = 1;
          }

        for (j = 1; j < 256; j++)
          if ( ! colored[j])
            break;

        if (j > clabels)
          { clabels = j;
            if (j > 255)
              { fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
                fprintf(stderr,"         Returning NULL result\n");
                fflush(stderr);
                return (NULL);
              }
          }

        cverts[c].depth = j;
        for (e = cverts[c].first; e < cverts[c+1].first; e++)
          { d = cedges[e].head;
            if (d < c)
              colored[cverts[d].depth] = 0;
          }
      }
    
    result->nlabels = clabels;
    result->labels  = Make_Image(GREY,shed->labels->width,shed->labels->height);
  
    brush->op = SET;
    for (c = 0; c < nbasins; c++)
      { brush->red = cverts[labels[c]].depth / 255. + .0001;
        Draw_Watershed(shed,c,brush,result->labels);
      }
  }

#ifdef DEBUG_COLLAPSE_PRINT
  { int i, j;

    for (j = 0; j < cheight; j++)
      { for (i = 0; i < cwidth; i++)
          printf(" %d",result->labels->array[j*cwidth+i]);
        printf("\n");
      }
  }
#endif

  { int c, b;
    int area;

    area = shed->labels->width * shed->labels->height;

    for (c = 0; c < cbasins; c++)
      { cverts[c].seed  = area;
        cverts[c].size  = 0;
        cverts[c].depth = 0x10000;
      }

    for (b = 0; b < nbasins; b++)
      { c = labels[b];
        if (verts[b].seed < cverts[c].seed)
          cverts[c].seed = verts[b].seed;
        if (verts[b].depth < cverts[c].depth)
          cverts[c].depth = verts[b].depth;
        cverts[c].size += verts[b].size;
      }
  }

#ifdef DEBUG_COLLAPSE
  { int i, j;

    printf("\nCollapsed graph:\n");
    for (i = 0; i < cbasins; i++)
      { j = cverts[i].seed;
        printf("Basin %d[%d,%d] -> ",i,j%cwidth,j/cwidth);
        for (j = cverts[i].first; j < cverts[i+1].first; j++)
          printf(" %d(%d)",cedges[j].head,cedges[j].height);
        printf("\n");
      }
  }
#endif

  return (result);
}


/**************************************************************************************\
*                                                                                      *
*  FLOOD-FILL ROUTINE CALLED FROM draw_lib.  Its here because it needs the chord and   *
*    out arrays that are working structures for the watershed algorithms (which are    *
*    also flood fill based).                                                           *
*                                                                                      *
\**************************************************************************************/

void fill_region(void (*painter)(int), Image *canvas, int iscon4, int (*value)(int), int seed)
{ int  *chord, *out;
 
  setup_globals_2D(canvas->width,canvas->height,iscon4);

  chord  = get_chord(2*carea,"fill_region");
  out    = chord + carea;

  { int p, v, *b, j, q;
    int qbot, qtop;

    qbot = -1;
    PUSH(seed);
    out[seed] = 1;
    while (qbot >= 0)     //  BFS flood on chord queue, paint as you go
      { POP(p);

        painter(p);

        b = boundary_pixels_2d(p);
        for (j = 0; j < n_nbrs4; j++)
          if (b[j])
            { q = p + neighbor4[j];
              if (value(q) && out[q] == 0)
                { PUSH(q);
                  out[q] = 1;
                }
            }
      }

    PUSH(seed);           //  BFS flood again to reset all the marks
    out[seed] = 0;
    while (qbot >= 0)
      { POP(p);
        chord[p] = 0;
        b = boundary_pixels_2d(p);
        for (j = 0; j < n_nbrs4; j++)
          if (b[j])
            { q = p + neighbor4[j];
              if (value(q) && out[q] == 1)
                { PUSH(q);
                  out[q] = 0;
                }
            }
      }
  }
}


/**************************************************************************************\
*                                                                                      *
*  PRODUCE A COLORED IMAGE OF A WATERSHED DECOMPOSITION                                *
*                                                                                      *
\**************************************************************************************/

static int black[3]    = {   0,   0,   0 };
static int red[3]      = { 255,   0,   0 };
static int green[3]    = {   0, 255,   0 };
static int blue[3]     = {   0,   0, 255 };
static int yellow[3]   = { 255, 255,   0 };
static int cyan[3]     = {   0, 255, 255 };
static int magenta[3]  = { 255,   0, 255 };
static int orange[3]   = { 255, 127,   0 };
static int brown[3]    = { 255,  65, 127 };
static int *palette[9] = { red, green, yellow, cyan, orange, brown, blue, magenta };

Image *Color_Watersheds(Watershed_2D *shed, Image *source)
{ int     q, p, i, j, v;
  int    *color, carea;
  Image  *view;
  uint8  *iarray;
  uint8  *value8;
  uint16 *value16;

  view  = Make_Image(COLOR,shed->labels->width,shed->labels->height);
  carea = shed->labels->width * shed->labels->height;

  iarray = shed->labels->array;
  if (source->kind == GREY)
    value8 = source->array;
  else
    value16 = (unsigned short *) (source->array);

  q = 0;
  for (p = 0; p < carea; p++)
    { if (iarray[p] <= 0)
        { printf("Unassigned pix\n");
          continue;
        }
      color = palette[(iarray[p]-1)%8];

      if (value8 == NULL)
        v = (value16[p] >> 8);
      else
        v = value8[p];

      v = (v >> 2) * 3;
      view->array[q]   = (color[0] >> 2) + v;
      view->array[q+1] = (color[1] >> 2) + v;
      view->array[q+2] = (color[2] >> 2) + v;

      q += 3;
    }

  return (view);
}


/**************************************************************************************\
*                                                                                      *
*  BUILD A 3D WATERSHED DECOMPOSITION  (UNTESTED !!!)                                  *
*                                                                                      *
\**************************************************************************************/


static inline int watershed_3d_ssize(Watershed_3D *part)
{ return (sizeof(int)*(part->nbasins+1)); }


typedef struct __Watershed_3D
  { struct __Watershed_3D *next;
    int                    ssize;
    Watershed_3D           watershed_3d;
  } _Watershed_3D;

static _Watershed_3D *Free_Watershed_3D_List = NULL;
static int    Watershed_3D_Offset, Watershed_3D_Inuse;

static inline void allocate_watershed_3d_seeds(Watershed_3D *watershed_3d, int ssize, char *routine)
{ _Watershed_3D *object  = (_Watershed_3D *) (((char *) watershed_3d) - Watershed_3D_Offset);
  if (object->ssize < ssize)
    { if (object->ssize == 0)
        object->watershed_3d.seeds = NULL;
      object->watershed_3d.seeds  = Guarded_Realloc(object->watershed_3d.seeds,ssize,routine);
      object->ssize = ssize;
    }
}

static inline Watershed_3D *new_watershed_3d(int ssize, char *routine)
{ _Watershed_3D *object;

  if (Free_Watershed_3D_List == NULL)
    { object = (_Watershed_3D *) Guarded_Malloc(sizeof(_Watershed_3D),routine);
      Watershed_3D_Offset = ((char *) &(object->watershed_3d)) - ((char *) object);
      object->ssize = 0;
    }
  else
    { object = Free_Watershed_3D_List;
      Free_Watershed_3D_List = object->next;
    }
  Watershed_3D_Inuse += 1;
  allocate_watershed_3d_seeds(&(object->watershed_3d),ssize,routine);
  object-> watershed_3d.labels = NULL;
  return (&(object->watershed_3d));
}

static inline Watershed_3D *copy_watershed_3d(Watershed_3D *watershed_3d)
{ Watershed_3D *copy = new_watershed_3d(watershed_3d_ssize(watershed_3d),"Copy_Watershed_3D");
  Watershed_3D  temp = *copy;
  *copy = *watershed_3d;
  copy->seeds = temp.seeds;
  if (watershed_3d_ssize(watershed_3d) != 0)
    memcpy(copy->seeds,watershed_3d->seeds,watershed_3d_ssize(watershed_3d));
  if (watershed_3d->labels != NULL)
    copy->labels = Copy_Stack(temp.labels);
  return (copy);
}

Watershed_3D *Copy_Watershed_3D(Watershed_3D *watershed_3d)
{ return (copy_watershed_3d(watershed_3d)); }

static inline void pack_watershed_3d(Watershed_3D *watershed_3d)
{ _Watershed_3D *object  = (_Watershed_3D *) (((char *) watershed_3d) - Watershed_3D_Offset);
  if (object->ssize > watershed_3d_ssize(watershed_3d))
    { object->ssize = watershed_3d_ssize(watershed_3d);
      if (object->ssize != 0)
        object->watershed_3d.seeds = Guarded_Realloc(object->watershed_3d.seeds,
                                                     object->ssize,"Pack_Watershed_3D");
      else
        { free(object->watershed_3d.seeds);
          object->ssize = 0;
        }
    }
  if (watershed_3d->labels != NULL)
    Pack_Stack(watershed_3d->labels);
}

void Pack_Watershed_3D(Watershed_3D *watershed_3d)
{ pack_watershed_3d(watershed_3d); }

static inline void free_watershed_3d(Watershed_3D *watershed_3d)
{ _Watershed_3D *object  = (_Watershed_3D *) (((char *) watershed_3d) - Watershed_3D_Offset);
  object->next = Free_Watershed_3D_List;
  Free_Watershed_3D_List = object;
  if (watershed_3d->labels != NULL)
    Free_Stack(watershed_3d->labels);
  Watershed_3D_Inuse -= 1;
}

void Free_Watershed_3D(Watershed_3D *watershed_3d)
{ free_watershed_3d(watershed_3d); }

static inline void kill_watershed_3d(Watershed_3D *watershed_3d)
{ _Watershed_3D *object  = (_Watershed_3D *) (((char *) watershed_3d) - Watershed_3D_Offset);
  if (watershed_3d->labels != NULL)
    Kill_Stack(watershed_3d->labels);
  if (object->ssize != 0)
    free(watershed_3d->seeds);
  free(((char *) watershed_3d) - Watershed_3D_Offset);
  Watershed_3D_Inuse -= 1;
}

void Kill_Watershed_3D(Watershed_3D *watershed_3d)
{ kill_watershed_3d(watershed_3d); }

static inline void reset_watershed_3d()
{ _Watershed_3D *object;
  while (Free_Watershed_3D_List != NULL)
    { object = Free_Watershed_3D_List;
      Free_Watershed_3D_List = object->next;
      kill_watershed_3d(&(object->watershed_3d));
      Watershed_3D_Inuse += 1;
    }
}

int Watershed_3D_Usage()
{ return (Watershed_3D_Inuse); }

void Reset_Watershed_3D()
{ reset_watershed_3d();
  get_chord(-1,NULL);
}

static int cdepth;     //  Other globals, e.g. cwidth, etc. already introduced for 2D code
static int cvolume;
static int chk_depth;
static int chk_iscon6;

static int neighbor6[26];
static int n_nbrs6;

static inline int *boundary_pixels_3d(int p)
{ static int bound[26];
  int x, xn, xp;
  int y, yn, yp;
  int z, zn, zp;

  z = p/carea;
  y = p%carea;
  x = y%cwidth;
  y = y/cwidth;

  yn = (y > 0);
  xn = (x > 0);
  zn = (z > 0);
  xp = (x < chk_width);
  yp = (y < chk_height);
  zp = (z < chk_depth);
  if (chk_iscon6)
    { bound[0] = zn;
      bound[1] = yn;
      bound[2] = xn;
      bound[3] = xp;
      bound[4] = yp;
      bound[5] = zp;
    }
  else
    { int t;

      bound[ 1] = t = zn && yn;
      bound[ 0] = t && xn;
      bound[ 2] = t && xp;
      bound[ 3] = zn && xn;
      bound[ 4] = zn;
      bound[ 5] = zn && xp;
      bound[ 7] = t = zn && yp;
      bound[ 6] = t && xn;
      bound[ 8] = t && xp;

      bound[ 9] = yn && xn;
      bound[10] = yn;
      bound[11] = yn && xp;
      bound[12] = xn;
      bound[13] = xp;
      bound[14] = yp && xn;
      bound[15] = yp;
      bound[16] = yp && xp;

      bound[18] = t = zp && yn;
      bound[17] = t && xn;
      bound[19] = t && xp;
      bound[20] = zp && xn;
      bound[21] = zp;
      bound[22] = zp && xp;
      bound[24] = t = zp && yp;
      bound[23] = t && xn;
      bound[25] = t && xp;
    }

  return (bound);
}

static void setup_globals_3D(int width, int height, int depth, int iscon6)
{ cwidth   = width;
  cheight  = height;
  cdepth   = depth;
  carea    = cwidth*cheight;
  cvolume  = carea*cdepth;

  if (iscon6)
    { neighbor6[0] = -carea;
      neighbor6[1] = -cwidth;
      neighbor6[2] = -1;
      neighbor6[3] =  1;
      neighbor6[4] =  cwidth;
      neighbor6[5] =  carea;
      n_nbrs6 = 6;
    }
  else
    { int x, y, z;

      n_nbrs6 = 0;
      for (z = -carea; z <= carea; z += carea)
        for (y = -cwidth; y <= cwidth; y += cwidth)
          for (x = -1; x <= 1; x += 1)
            if (! (x==0 && y == 0 && z == 0))
              neighbor6[n_nbrs6++] = z+y+x;
    }

  chk_width  = cwidth-1;
  chk_height = cheight-1;
  chk_depth  = cdepth-1;
  chk_iscon6 = iscon6;
}

Watershed_3D *Build_3D_Watershed(Stack *frame, int iscon6)
{ Watershed_3D *result;

  int    nlabels;   // The number of labels (1..nlabels)
  int    nbasins;   // Number of watersheds
  int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

  int index[0x10001];
  int maxval;

  uint8  *value8;
  uint16 *value16;

  int   *chord;
  int   *out;

  if (frame->kind == COLOR)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
      exit (1);
    }
  if (frame->kind == FLOAT32)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
      exit (1);
    }

  setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);

  chord = get_chord(2*cvolume,"Build_3D_Watershed");
  out   = chord + cvolume;

  if (frame->kind == GREY16)
    { maxval  = 0x10000;
      value8  = NULL;
      value16 = (uint16 *) (frame->array);
    }
  else
    { maxval  = 0x100;
      value8  = frame->array;
      value16 = NULL;
    }

  // Establish lists of pixels with each value

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      for (p = 0; p < cvolume; p++)
        { v        = value8[p];
          chord[p] = index[v];
          index[v] = p;
        }
    else  // value16 != NULL
      for (p = 0; p < cvolume; p++)
        { v        = value16[p];
          chord[p] = index[v];
          index[v] = p;
        }
  }

#define PUSH(p)			\
{ if (qbot < 0)			\
    qtop = qbot = (p);		\
  else				\
    qtop = chord[qtop] = (p);	\
  chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

  { int v, p;

    for (p = 0; p < cvolume; p++)
      out[p] = MINIM;

    if (value8 != NULL)

      for (v = 0; v < maxval; v++)
        { int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

          qbot = -1;
          for (p = index[v]; p >= 0; p = t)
            { int *b, j;
  
              t = chord[p];
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j] && value8[p+neighbor6[j]] < v)
                  { out[p] = INIT;
                    PUSH(p)
                    break;
                  }
            }

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
          while (qbot >= 0)
            { int p, *b, j, q;
              int cdist;
  
              POP(p);
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { q = p + neighbor6[j];
                    if (out[q] != INIT && value8[q] == v)
                      { out[q] = INIT;
                        PUSH(q);
                      }
                  }
            }
        }
  }

  // Build pixel value lists again but now with only INIT pixels in each list

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value8[p];
              chord[p] = index[v];
              index[v] = p;
            }
       }

     else // value16 != NULL
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value16[p];
              chord[p] = index[v];
              index[v] = p;
            }
      }
  }

  // Find connected components of MINIM pixels and label with a unique catchment basin #.

  { int v, p;
    int qtop, qbot;

    nbasins = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == MINIM)
        { PUSH(p);
          out[p] = ++nbasins;
          while (qbot >= 0)
            { int q, *b, j, r;

              POP(q);
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q+neighbor6[j];
                    if (out[r] == MINIM)
                      { PUSH(r);
                        out[r] = nbasins;
                      }
                  }
            }
        }
  }

  { int v, p;

    for (v = 0; v < maxval; v++)
      { int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

        qbot = -1;
        for (p = index[v]; p >= 0; p = t)
          { int *b, j;

            t      = chord[p];
            out[p] = MASK;
            b      = boundary_pixels_3d(p);
            for (j = 0; j < n_nbrs6; j++)
              if (b[j] && out[p+neighbor6[j]] >= 0)
                { out[p] = ONQUEUE;
                  PUSH(p);
                  break;
                }
          }

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

        while (qbot >= 0)
          { int p, *b, j, q, o;

            POP(p);
            b = boundary_pixels_3d(p);

            for (j = 0; j < n_nbrs6; j++)
              if (b[j])
                { q = p + neighbor6[j];
                  o = out[q];
                  if (o == MASK)
                    { out[q] = ONQUEUE;
                      PUSH(q);
                    }
                  else if (o != ONQUEUE && o != INIT)
                    { if (out[p] == ONQUEUE)
                        out[p] = o;
                      else if (o != out[p] && o != WSHED)
                        if (out[p] == WSHED)
                          out[p] = o;
                        else
                          out[p] = WSHED;
                    }
                }
          }
      }
  }

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

  { int p, qbot, qtop;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == WSHED)
        { int *b, j;

          b = boundary_pixels_3d(p);
          for (j = 0; j < n_nbrs6; j++)
            if (b[j] && out[p+neighbor6[j]] > 0)
              { PUSH(p);
                out[p] = ONQUEUE;
                break;
              }
        }

    while (qbot >= 0)
      { int *b, j, q, o;

        POP(p);
        b = boundary_pixels_3d(p);

        for (j = 0; j < n_nbrs6; j++)
          if (b[j])
            { q = p + neighbor6[j];
              o = out[q];
              if (o > 0)
                out[p] = o;
              else if (out[q] != ONQUEUE)
                { PUSH(q);
                  out[q] = ONQUEUE;
                }
            }
      }
  }

  //  Make the data structure now

  result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
  result->nbasins = nbasins;
  result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

  { int   qbot, qtop;
    int   olabel;
    int   adjacent[256];
    int   alargest;
    int   p;

    seeds   = result->seeds;
    nlabels = 0;

    for (p = 1; p < 256; p++)
      adjacent[p] = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] > 0)

        { PUSH(p);
          olabel   = out[p];
          out[p]   = 0;
          alargest = 0;

          while (qbot >= 0)
            { int *b, j, q, r, o;

              POP(q);
              out[q] = -olabel;
      
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q + neighbor6[j];
                    o = out[r];
                    if (o == olabel)
                      { PUSH(r);
                        out[r] = 0;
                      }
                    else
                      { o = -o;
                        if (o > 0 && o != olabel)
                          { o = seeds[o];
                            adjacent[o] = 1;
                            if (o > alargest)
                              alargest = o;
                          }
                      }
                  }
            }

          { int j;

            for (j = 1; j <= alargest; j++)
              if (! adjacent[j])
                break;

            seeds[olabel] = j;
            if (j > nlabels)
              { nlabels = j;
                if (j > 255)
                  { fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
                    fprintf(stderr,"         Returning NULL result\n");
                    fflush(stderr);
                    return (NULL);
                  }
              }

            for (j = 1; j <= alargest; j++)
              adjacent[j] = 0;
          }
        }
  }

  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

  { int    p, b;
    uint8 *iarray;

    result->nlabels = nlabels;
    result->labels  = Make_Stack(GREY,cwidth,cheight,cdepth);
    iarray = result->labels->array;

    for (p = 0; p < cvolume; p++)
      { out[p] = b = -out[p];
        iarray[p] = seeds[b];
      }
    for (p = 0; p < cvolume; p++)
      if (seeds[b = out[p]]  > 0) 
        seeds[b] = -(p+1); 
    for (b = 1; b <= nbasins; b++)
      seeds[b] = -seeds[b];
  }

  { int    p;
    int   *color;
    Stack *view;

    view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

    for (p = 0; p < cvolume; p++)
      { if (out[p] <= 0)
          { printf("Unassigned pix\n");
            continue;
          }
        color = palette[(out[p]-1)%8];
        view->array[3*p]   = color[0];
        view->array[3*p+1] = color[1];
        view->array[3*p+2] = color[2];
      }

    Write_Stack("watershed.tif",view);

    Kill_Stack(view);
  }

  return (result);
}

Stack *Build_3D_Watershed2(Stack *frame, int iscon6)
{ 
  int    nlabels;   // The number of labels (1..nlabels)
  int    nbasins;   // Number of watersheds
  int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

  int index[0x10001];
  int maxval;

  uint8  *value8;
  uint16 *value16;

  int   *chord;
  int   *out;

  if (frame->kind == COLOR)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
      exit (1);
    }
  if (frame->kind == FLOAT32)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
      exit (1);
    }

  setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);
printf("%d\n", cvolume);

  chord = get_chord(2*cvolume,"Build_3D_Watershed");
  out   = chord + cvolume;

  if (frame->kind == GREY16)
    { maxval  = 0x10000;
      value8  = NULL;
      value16 = (uint16 *) (frame->array);
    }
  else
    { maxval  = 0x100;
      value8  = frame->array;
      value16 = NULL;
    }

  // Establish lists of pixels with each value

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      for (p = 0; p < cvolume; p++)
        { v        = value8[p];
          chord[p] = index[v];
          index[v] = p;
        }
    else  // value16 != NULL
      for (p = 0; p < cvolume; p++)
        { v        = value16[p];
          chord[p] = index[v];
          index[v] = p;
        }
  }

#define PUSH(p)			\
{ if (qbot < 0)			\
    qtop = qbot = (p);		\
  else				\
    qtop = chord[qtop] = (p);	\
  chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

  { int v, p;

    for (p = 0; p < cvolume; p++)
      out[p] = MINIM;

    if (value8 != NULL)

      for (v = 0; v < maxval; v++)
        { int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

          qbot = -1;
          for (p = index[v]; p >= 0; p = t)
            { int *b, j;
  
              t = chord[p];
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j] && value8[p+neighbor6[j]] < v)
                  { out[p] = INIT;
                    PUSH(p)
                    break;
                  }
            }

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
          while (qbot >= 0)
            { int p, *b, j, q;
              int cdist;
  
              POP(p);
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { q = p + neighbor6[j];
                    if (out[q] != INIT && value8[q] == v)
                      { out[q] = INIT;
                        PUSH(q);
                      }
                  }
            }
        }
  }

  // Build pixel value lists again but now with only INIT pixels in each list

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value8[p];
              chord[p] = index[v];
              index[v] = p;
            }
       }

     else // value16 != NULL
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value16[p];
              chord[p] = index[v];
              index[v] = p;
            }
      }
  }

  // Find connected components of MINIM pixels and label with a unique catchment basin #.

  { int v, p;
    int qtop, qbot;

    nbasins = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == MINIM)
        { PUSH(p);
          out[p] = ++nbasins;
          while (qbot >= 0)
            { int q, *b, j, r;

              POP(q);
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q+neighbor6[j];
                    if (out[r] == MINIM)
                      { PUSH(r);
                        out[r] = nbasins;
                      }
                  }
            }
        }
  }

  { int v, p;

    for (v = 0; v < maxval; v++)
      { int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

        qbot = -1;
        for (p = index[v]; p >= 0; p = t)
          { int *b, j;

            t      = chord[p];
            out[p] = MASK;
            b      = boundary_pixels_3d(p);
            for (j = 0; j < n_nbrs6; j++)
              if (b[j] && out[p+neighbor6[j]] >= 0)
                { out[p] = ONQUEUE;
                  PUSH(p);
                  break;
                }
          }

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

        while (qbot >= 0)
          { int p, *b, j, q, o;

            POP(p);
            b = boundary_pixels_3d(p);

            for (j = 0; j < n_nbrs6; j++)
              if (b[j])
                { q = p + neighbor6[j];
                  o = out[q];
                  if (o == MASK)
                    { out[q] = ONQUEUE;
                      PUSH(q);
                    }
                  else if (o != ONQUEUE && o != INIT)
                    { if (out[p] == ONQUEUE)
                        out[p] = o;
                      else if (o != out[p] && o != WSHED)
                        if (out[p] == WSHED)
                          out[p] = o;
                        else
                          out[p] = WSHED;
                    }
                }
          }
      }
  }

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

  { int p, qbot, qtop;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == WSHED)
        { int *b, j;

          b = boundary_pixels_3d(p);
          for (j = 0; j < n_nbrs6; j++)
            if (b[j] && out[p+neighbor6[j]] > 0)
              { PUSH(p);
                out[p] = ONQUEUE;
                break;
              }
        }

    while (qbot >= 0)
      { int *b, j, q, o;

        POP(p);
        b = boundary_pixels_3d(p);

        for (j = 0; j < n_nbrs6; j++)
          if (b[j])
            { q = p + neighbor6[j];
              o = out[q];
              if (o > 0)
                out[p] = o;
              else if (out[q] != ONQUEUE)
                { PUSH(q);
                  out[q] = ONQUEUE;
                }
            }
      }
  }

  // Copy the labels for output
  Stack * result = Make_Stack(GREY32, frame->width, frame->height, frame->depth);
{
	int * label_new_array = (int *) result->array;
	int p;
	for(p=0; p<cvolume; p++)
		label_new_array[p] = out[p];
}
return result;

/*
  //  Make the data structure now

  result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
  result->nbasins = nbasins;
  result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

  { int   qbot, qtop;
    int   olabel;
    int   adjacent[256];
    int   alargest;
    int   p;

    seeds   = result->seeds;
    nlabels = 0;

    for (p = 1; p < 256; p++)
      adjacent[p] = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] > 0)

        { PUSH(p);
          olabel   = out[p];
          out[p]   = 0;
          alargest = 0;

          while (qbot >= 0)
            { int *b, j, q, r, o;

              POP(q);
              out[q] = -olabel;
      
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q + neighbor6[j];
                    o = out[r];
                    if (o == olabel)
                      { PUSH(r);
                        out[r] = 0;
                      }
                    else
                      { o = -o;
                        if (o > 0 && o != olabel)
                          { o = seeds[o];
                            adjacent[o] = 1;
                            if (o > alargest)
                              alargest = o;
                          }
                      }
                  }
            }

          { int j;

            for (j = 1; j <= alargest; j++)
              if (! adjacent[j])
                break;

            seeds[olabel] = j;
            if (j > nlabels)
              { nlabels = j;
                if (j > 255)
                  { fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
                    fprintf(stderr,"         Returning NULL result\n");
                    fflush(stderr);
                    return (NULL);
                  }
              }

            for (j = 1; j <= alargest; j++)
              adjacent[j] = 0;
          }
        }
  }

  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

  { int    p, b;
    uint8 *iarray;

    result->nlabels = nlabels;
    result->labels  = Make_Stack(GREY,cwidth,cheight,cdepth);
    iarray = result->labels->array;

    for (p = 0; p < cvolume; p++)
      { out[p] = b = -out[p];
        iarray[p] = seeds[b];
      }
    for (p = 0; p < cvolume; p++)
      if (seeds[b = out[p]]  > 0) 
        seeds[b] = -(p+1); 
    for (b = 1; b <= nbasins; b++)
      seeds[b] = -seeds[b];
  }
*/
  { int    p;
    int   *color;
    Stack *view;

    view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

    for (p = 0; p < cvolume; p++)
      { if (out[p] <= 0)
          { printf("Unassigned pix\n");
            continue;
          }
        color = palette[(out[p]-1)%8];
        view->array[3*p]   = color[0];
        view->array[3*p+1] = color[1];
        view->array[3*p+2] = color[2];
      }

    Write_Stack("watershed.tif",view);

    Kill_Stack(view);
  }

  return (result);
}

Stack *Build_3D_Watershed3(Stack *frame, int iscon6)
{ 
  int    nlabels;   // The number of labels (1..nlabels)
  int    nbasins;   // Number of watersheds
  int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

  int index[0x10001];
  int maxval;

  uint8  *value8;
  uint16 *value16;

  int   *chord;
  int   *out;

  if (frame->kind == COLOR)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
      exit (1);
    }
  if (frame->kind == FLOAT32)
    { fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
      exit (1);
    }

  setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);

  chord = get_chord(2*cvolume,"Build_3D_Watershed");
  out   = chord + cvolume;

  if (frame->kind == GREY16)
    { maxval  = 0x10000;
      value8  = NULL;
      value16 = (uint16 *) (frame->array);
    }
  else
    { maxval  = 0x100;
      value8  = frame->array;
      value16 = NULL;
    }

  // Establish lists of pixels with each value

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      for (p = 0; p < cvolume; p++)
        { v        = value8[p];
          chord[p] = index[v];
          index[v] = p;
        }
    else  // value16 != NULL
      for (p = 0; p < cvolume; p++)
        { v        = value16[p];
          chord[p] = index[v];
          index[v] = p;
        }
  }

#define PUSH(p)			\
{ if (qbot < 0)			\
    qtop = qbot = (p);		\
  else				\
    qtop = chord[qtop] = (p);	\
  chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

  { int v, p;

    for (p = 0; p < cvolume; p++)
      out[p] = MINIM;

    if (value8 != NULL)

      for (v = 0; v < maxval; v++)
        { int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

          qbot = -1;
          for (p = index[v]; p >= 0; p = t)
            { int *b, j;
  
              t = chord[p];
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j] && value8[p+neighbor6[j]] < v)
                  { out[p] = INIT;
                    PUSH(p)
                    break;
                  }
            }

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
          while (qbot >= 0)
            { int p, *b, j, q;
              int cdist;
  
              POP(p);
              b = boundary_pixels_3d(p);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { q = p + neighbor6[j];
                    if (out[q] != INIT && value8[q] == v)
                      { out[q] = INIT;
                        PUSH(q);
                      }
                  }
            }
        }
  }

  // Build pixel value lists again but now with only INIT pixels in each list

  { int p, v;

    for (v = 0; v <= maxval; v++)
      index[v] = -1;

    if (value8 != NULL)
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value8[p];
              chord[p] = index[v];
              index[v] = p;
            }
       }

     else // value16 != NULL
      { for (p = 0; p < cvolume; p++)
          if (out[p] == INIT)
            { v        = value16[p];
              chord[p] = index[v];
              index[v] = p;
            }
      }
  }

  // Find connected components of MINIM pixels and label with a unique catchment basin #.

  { int v, p;
    int qtop, qbot;

    nbasins = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == MINIM)
        { PUSH(p);
          out[p] = ++nbasins;
          while (qbot >= 0)
            { int q, *b, j, r;

              POP(q);
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q+neighbor6[j];
                    if (out[r] == MINIM)
                      { PUSH(r);
                        out[r] = nbasins;
                      }
                  }
            }
        }
  }

  { int v, p;

    for (v = 0; v < maxval; v++)
      { int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

        qbot = -1;
        for (p = index[v]; p >= 0; p = t)
          { int *b, j;

            t      = chord[p];
            out[p] = MASK;
            b      = boundary_pixels_3d(p);
            for (j = 0; j < n_nbrs6; j++)
              if (b[j] && out[p+neighbor6[j]] >= 0)
                { out[p] = ONQUEUE;
                  PUSH(p);
                  break;
                }
          }

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

        while (qbot >= 0)
          { int p, *b, j, q, o;

            POP(p);
            b = boundary_pixels_3d(p);

            for (j = 0; j < n_nbrs6; j++)
              if (b[j])
                { q = p + neighbor6[j];
                  o = out[q];
                  if (o == MASK)
                    { out[q] = ONQUEUE;
                      PUSH(q);
                    }
		  else if(o > 0)
                    {
                      if(out[p] == ONQUEUE)
                        out[p] = o;
		      else
                        { if(out[p]!=o)
                          out[p] = WSHED;
                        }
                    }
/*
                  else if (o != ONQUEUE && o != INIT)
                    { if (out[p] == ONQUEUE)
                        out[p] = o;
		      else if (o != out[p] && o != WSHED)
                        if (out[p] == WSHED)
                          out[p] = o;
                        else
                          out[p] = WSHED;
                    }
*/
                }
          }
      }
  }
  /*
  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

  { int p, qbot, qtop;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] == WSHED)
        { int *b, j;

          b = boundary_pixels_3d(p);
          for (j = 0; j < n_nbrs6; j++)
            if (b[j] && out[p+neighbor6[j]] > 0)
              { PUSH(p);
                out[p] = ONQUEUE;
                break;
              }
        }

    while (qbot >= 0)
      { int *b, j, q, o;

        POP(p);
        b = boundary_pixels_3d(p);

        for (j = 0; j < n_nbrs6; j++)
          if (b[j])
            { q = p + neighbor6[j];
              o = out[q];
              if (o > 0)
                out[p] = o;
              else if (out[q] != ONQUEUE)
                { PUSH(q);
                  out[q] = ONQUEUE;
                }
            }
      }
  }
  */
  // Copy the labels for output
  Stack * result = Make_Stack(GREY32, frame->width, frame->height, frame->depth);
{
	int * label_new_array = (int *) result->array;
	int p;
	for(p=0; p<cvolume; p++)
		label_new_array[p] = out[p];
}
return result;

/*
  //  Make the data structure now

  result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
  result->nbasins = nbasins;
  result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

  { int   qbot, qtop;
    int   olabel;
    int   adjacent[256];
    int   alargest;
    int   p;

    seeds   = result->seeds;
    nlabels = 0;

    for (p = 1; p < 256; p++)
      adjacent[p] = 0;

    qbot = -1;
    for (p = 0; p < cvolume; p++)
      if (out[p] > 0)

        { PUSH(p);
          olabel   = out[p];
          out[p]   = 0;
          alargest = 0;

          while (qbot >= 0)
            { int *b, j, q, r, o;

              POP(q);
              out[q] = -olabel;
      
              b = boundary_pixels_3d(q);
              for (j = 0; j < n_nbrs6; j++)
                if (b[j])
                  { r = q + neighbor6[j];
                    o = out[r];
                    if (o == olabel)
                      { PUSH(r);
                        out[r] = 0;
                      }
                    else
                      { o = -o;
                        if (o > 0 && o != olabel)
                          { o = seeds[o];
                            adjacent[o] = 1;
                            if (o > alargest)
                              alargest = o;
                          }
                      }
                  }
            }

          { int j;

            for (j = 1; j <= alargest; j++)
              if (! adjacent[j])
                break;

            seeds[olabel] = j;
            if (j > nlabels)
              { nlabels = j;
                if (j > 255)
                  { fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
                    fprintf(stderr,"         Returning NULL result\n");
                    fflush(stderr);
                    return (NULL);
                  }
              }

            for (j = 1; j <= alargest; j++)
              adjacent[j] = 0;
          }
        }
  }

  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

  { int    p, b;
    uint8 *iarray;

    result->nlabels = nlabels;
    result->labels  = Make_Stack(GREY,cwidth,cheight,cdepth);
    iarray = result->labels->array;

    for (p = 0; p < cvolume; p++)
      { out[p] = b = -out[p];
        iarray[p] = seeds[b];
      }
    for (p = 0; p < cvolume; p++)
      if (seeds[b = out[p]]  > 0) 
        seeds[b] = -(p+1); 
    for (b = 1; b <= nbasins; b++)
      seeds[b] = -seeds[b];
  }
*/
  { int    p;
    int   *color;
    Stack *view;

    view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

    for (p = 0; p < cvolume; p++)
      { if (out[p] <= 0)
          { printf("Unassigned pix\n");
            continue;
          }
        color = palette[(out[p]-1)%8];
        view->array[3*p]   = color[0];
        view->array[3*p+1] = color[1];
        view->array[3*p+2] = color[2];
      }

    Write_Stack("watershed.tif",view);

    Kill_Stack(view);
  }

  return (result);
}

/**************************************************************************************\
*                                                                                      *
*  PRODUCE A COLORED STACK OF A 3D WATERSHED DECOMPOSITION                             *
*                                                                                      *
\**************************************************************************************/
Stack *Color_Watersheds_3D(Watershed_3D *shed, Stack *source)
{ int     q, p, i, j, v;
  int    *color, carea;
	Stack *view;
	uint8  *iarray;
	uint8  *value8;
	uint16 *value16;

	view  = Make_Stack(COLOR,shed->labels->width,shed->labels->height,shed->labels->depth);
	carea = shed->labels->width * shed->labels->height * shed->labels->depth;

	iarray = shed->labels->array;
	if (source->kind == GREY)
			value8 = source->array;
	else
			value16 = (unsigned short *) (source->array);

	q = 0;
	for (p = 0; p < carea; p++)
	{ if (iarray[p] <= 0)
	{ printf("Unassigned pix\n");
	continue;
	}
	color = palette[(iarray[p]-1)%8];

	if (value8 == NULL)
			v = (value16[p] >> 8);
	else
			v = value8[p];

	v = (v >> 2) * 3;
	view->array[q]   = (color[0] >> 2) + v;
	view->array[q+1] = (color[1] >> 2) + v;
	view->array[q+2] = (color[2] >> 2) + v;

	q += 3;
	}

	return (view);
}

/**************************************************************************************\
*                                                                                      *
*  BUILD A 3D WATERSHED DECOMPOSITION SUCH THAT IT IS CONFINED TO A USER DEFINED 			 *
*	 SUB-VOLUME MASK (POSSIBLY IRREGULAR)																								 *		
*                                                                                      *
\**************************************************************************************/
		Watershed_3D *Build_3D_Watershed_Masked(Stack *frame, int iscon6, Stack * mask, unsigned int mask_label)
{ Watershed_3D *result;

	int    nlabels;   // The number of labels (1..nlabels)
	int    nbasins;   // Number of watersheds
	int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

	int index[0x10001];
	int maxval;

	uint8  *value8;
	uint16 *value16;

	int   *chord;
	int   *out;

	if (frame->kind == COLOR)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
	exit (1);
	}
	if (frame->kind == FLOAT32)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
	exit (1);
	}
	if (mask->kind != GREY && mask->kind != GREY32)
	{ fprintf(stderr,"Mask should be UINT8 or UINT32\n");
	exit (1);
	}

	setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);

	chord = get_chord(2*cvolume,"Build_3D_Watershed");
	out   = chord + cvolume;

	if (frame->kind == GREY16)
	{ maxval  = 0x10000;
	value8  = NULL;
	value16 = (uint16 *) (frame->array);
	}
	else
	{ maxval  = 0x100;
	value8  = frame->array;
	value16 = NULL;
	}
	
  // Establish lists of pixels with each value

	{ int p, v;

	for (v = 0; v <= maxval; v++)
		index[v] = -1;

	if (value8 != NULL)
				for (p = 0; p < cvolume; p++)
	{ if(mask->array[p]!=mask_label)
			continue;
		v        = value8[p];
	chord[p] = index[v];
	index[v] = p;
	}
	else  // value16 != NULL
				for (p = 0; p < cvolume; p++)
	{ if(mask->array[p]!=mask_label)
			continue;
		v        = value16[p];
	chord[p] = index[v];
	index[v] = p;
	}
	}

#define PUSH(p)			\
{ if (qbot < 0)			\
			qtop = qbot = (p);		\
					else				\
					qtop = chord[qtop] = (p);	\
					chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

{ int v, p;

    for (p = 0; p < cvolume; p++)
			if(mask->array[p]!=mask_label)
				out[p] = OUTBOUND;
			else
				out[p] = MINIM;

		if (value8 != NULL)

					for (v = 0; v < maxval; v++)
		{ int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

							qbot = -1;
					for (p = index[v]; p >= 0; p = t)
					{ int *b, j;
					if(mask->array[p]!=mask_label)
						continue;
					t = chord[p];
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
							if (b[j] && value8[p+neighbor6[j]] < v)
					{ out[p] = INIT;
					PUSH(p)
							break;
					}
					}

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
							while (qbot >= 0)
					{ int p, *b, j, q;
					int cdist;
  
					POP(p);
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
							if (b[j])
					{ q = p + neighbor6[j];
					if (out[q] != INIT && value8[q] == v)
					{ out[q] = INIT;
					PUSH(q);
					}
					}
					}
		}
}

  // Build pixel value lists again but now with only INIT pixels in each list

{ int p, v;

    for (v = 0; v <= maxval; v++)
			index[v] = -1;

		if (value8 != NULL)
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value8[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}

		else // value16 != NULL
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value16[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}
}

  // Find connected components of MINIM pixels and label with a unique catchment basin #.

{ int v, p;
    int qtop, qbot;

		nbasins = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == MINIM)
		{ PUSH(p);
		out[p] = ++nbasins;
		while (qbot >= 0)
		{ int q, *b, j, r;

		POP(q);
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ r = q+neighbor6[j];
		if (out[r] == MINIM)
		{ PUSH(r);
		out[r] = nbasins;
		}
		}
		}
		}
}

{ int v, p;

    for (v = 0; v < maxval; v++)
		{ int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

						qbot = -1;
				for (p = index[v]; p >= 0; p = t)
				{ int *b, j;

				t      = chord[p];
				out[p] = MASK;
				b      = boundary_pixels_3d(p);
				for (j = 0; j < n_nbrs6; j++)
						if (b[j] && out[p+neighbor6[j]] >= 0)
				{ out[p] = ONQUEUE;
				PUSH(p);
				break;
				}
				}

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

						while (qbot >= 0)
				{ int p, *b, j, q, o;

				POP(p);
				b = boundary_pixels_3d(p);

				for (j = 0; j < n_nbrs6; j++)
						if (b[j])
				{ q = p + neighbor6[j];
				o = out[q];
				if (o == MASK)
				{ out[q] = ONQUEUE;
				PUSH(q);
				}
				else if (o != ONQUEUE && o != INIT)
				{ if (out[p] == ONQUEUE)
					out[p] = o;
					else if (o != out[p] && o != WSHED)
							if (out[p] == WSHED)
							out[p] = o;
					else
							out[p] = WSHED;
				}
				}
				}
		}
}

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

{ int p, qbot, qtop;

    qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == WSHED)
		{ int *b, j;

		b = boundary_pixels_3d(p);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j] && out[p+neighbor6[j]] > 0)
		{ PUSH(p);
		out[p] = ONQUEUE;
		break;
		}
		}

		while (qbot >= 0)
		{ int *b, j, q, o;

		POP(p);
		b = boundary_pixels_3d(p);

		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ q = p + neighbor6[j];
		o = out[q];
		if (o > 0)
				out[p] = o;
		else if (out[q] != ONQUEUE)
		{ PUSH(q);
		out[q] = ONQUEUE;
		}
		}
		}
}

  //  Make the data structure now

			result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
	result->nbasins = nbasins;
	result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

{ int   qbot, qtop;
    int   olabel;
		int   adjacent[256];
		int   alargest;
		int   p;

		seeds   = result->seeds;
		nlabels = 0;

		for (p = 1; p < 256; p++)
				adjacent[p] = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] > 0)

		{ PUSH(p);
		olabel   = out[p];
		out[p]   = 0;
		alargest = 0;

		while (qbot >= 0)
		{ int *b, j, q, r, o;

		POP(q);
		out[q] = -olabel;
      
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ r = q + neighbor6[j];
		o = out[r];
		if (o == olabel)
		{ PUSH(r);
		out[r] = 0;
		}
		else
		{ o = -o;
		if (o > 0 && o != olabel)
		{ o = seeds[o];
		adjacent[o] = 1;
		if (o > alargest)
				alargest = o;
		}
		}
		}
		}

		{ int j;

		for (j = 1; j <= alargest; j++)
			if (! adjacent[j])
				break;

		seeds[olabel] = j;
		if (j > nlabels)
		{ nlabels = j;
		if (j > 255)
		{ fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
		fprintf(stderr,"         Returning NULL result\n");
		fflush(stderr);
		return (NULL);
		}
		}

		for (j = 1; j <= alargest; j++)
			adjacent[j] = 0;
		}
		}
}

  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

{ int    p, b;
    uint8 *iarray;

		result->nlabels = nlabels;
		result->labels  = Make_Stack(GREY,cwidth,cheight,cdepth);
		iarray = result->labels->array;

		for (p = 0; p < cvolume; p++)
		{ out[p] = b = -out[p];
		iarray[p] = seeds[b];
		}
		for (p = 0; p < cvolume; p++)
			if (seeds[b = out[p]]  > 0) 
				seeds[b] = -(p+1); 
		for (b = 1; b <= nbasins; b++)
						seeds[b] = -seeds[b];
}

{ int    p;
    int   *color;
		Stack *view;

		view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

		for (p = 0; p < cvolume; p++)
		{ if (out[p] <= 0)
		{ printf("Unassigned pix\n");
		continue;
		}
		view->array[3*p]   = out[p]==0?0:(((out[p])%255) + 1);
		view->array[3*p+1] = out[p]==0?0:(((out[p])%255) + 1);
		view->array[3*p+2] = out[p]==0?0:(((out[p])%255) + 1);
		}

		Write_Stack("watershed.tif",view);

		Kill_Stack(view);
}

  return (result);
}

/**************************************************************************************\
*                                                                                      *
*  BUILD A SEEDED 3D WATERSHED DECOMPOSITION SO THAT IT IS CONFINED TO A USER DEFINED  *
*	 SUB-VOLUME MASK (POSSIBLY IRREGULAR)																								 *		
*                                                                                      *
\**************************************************************************************/
		Watershed_3D *Build_3D_Watershed_Masked_Seeded(Stack *frame, int image_stack_scale, int iscon6, unsigned short int ** superpixel_maps, int original_width, double ** superpixel_2_seg_maps, unsigned int * pmap, unsigned int mask_label, Image * seed_image, int seed_depth, Stack * label_new)
{ Watershed_3D *result;

	int    nlabels;   // The number of labels (1..nlabels)
			int    nbasins;   // Number of watersheds
			int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

			int index[0x10001];
	int maxval;

	//(P)%plane_area
	
	int plane_area = frame->width * frame->height;
	int frame_width = frame->width;
	int frame_height = frame->height;
#define PLANE_MAPPING(Q) (image_stack_scale*(original_width*((Q)/frame_width) + ((Q)%frame_width)))
#define BODY_ID(P) (pmap[(int)(superpixel_2_seg_maps[(P)/plane_area][superpixel_maps[(P)/plane_area][PLANE_MAPPING((P)%plane_area)]])])	
#define IS_BODY_VOXEL(P) (BODY_ID(P)==mask_label || BODY_ID(P)==0)
#define IS_NOT_BODY_VOXEL(P) (BODY_ID(P)!=mask_label && BODY_ID(P)!=0)
			
	uint8  *value8;
	uint16 *value16;

	int   *chord;
	int   *out;

	if (frame->kind == COLOR)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
	exit (1);
	}
	if (frame->kind == FLOAT32)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
	exit (1);
	}

	setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);
	/*
	{
	FILE * fo = fopen("debug.txt", "w");
	int p;
	for(p=0; p<10; p++)
		fprintf(fo, "%u ", pmap[(int)(superpixel_2_seg_maps[p][superpixel_maps[p][200]])]);
	fclose(fo);
	}
	*/
	
	chord = get_chord(2*cvolume,"Build_3D_Watershed");
	out   = chord + cvolume;

	if (frame->kind == GREY16)
	{ maxval  = 0x10000;
	value8  = NULL;
	value16 = (uint16 *) (frame->array);
	}
	else
	{ maxval  = 0x100;
	value8  = frame->array;
	value16 = NULL;
	}
	
  // Establish lists of pixels with each value

	{ int v, p;

	for (p = 0; p < cvolume; p++)
		if(IS_NOT_BODY_VOXEL(p))
			out[p] = OUTBOUND;
	else
					out[p] = MINIM;
	}
	
	{ int p, v;

	for (v = 0; v <= maxval; v++){
		index[v] = -1;
	}

	if (value8 != NULL)
				for (p = 0; p < cvolume; p++)
	{ //if(IS_NOT_BODY_VOXEL(p))
		//  continue;
		if(out[p]==OUTBOUND)
			continue;
		v        = value8[p];
		chord[p] = index[v];
		index[v] = p;
	}
	else  // value16 != NULL
				for (p = 0; p < cvolume; p++)
	{// if(IS_NOT_BODY_VOXEL(p))
		//	continue;
		if(out[p]==OUTBOUND)
			continue;
		v        = value16[p];
		chord[p] = index[v];
		index[v] = p;
	}
	}
	
#define PUSH(p)			\
{ if (qbot < 0)			\
			qtop = qbot = (p);		\
					else				\
					qtop = chord[qtop] = (p);	\
					chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

{ int v, p;

		
if (value8 != NULL)

						for (v = 0; v < maxval; v++)
		{ int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

							qbot = -1;
					for (p = index[v]; p >= 0; p = t)
					{ int *b, j;
					t = chord[p];
					//if(IS_NOT_BODY_VOXEL(p))
					//		continue;
					if(out[p]==OUTBOUND)
						continue;
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
								if (b[j] && value8[p+neighbor6[j]] < v)
					{ out[p] = INIT;
					PUSH(p)
							break;
					}
					}

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
							while (qbot >= 0)
					{ int p, *b, j, q;
					int cdist;
  
					POP(p);
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
							if (b[j])
					{ q = p + neighbor6[j];
					if (out[q] != INIT && out[q] != OUTBOUND && value8[q] == v)
					{ out[q] = INIT;
					PUSH(q);
					}
					}
					}
		}
		
	{	// Force connected components in the seed image to have connected components of MINIM labels
		uint8 * seed_array = seed_image->array;
		int start_seed_mask = frame->width*frame->height*seed_depth;
		int end_seed_mask = frame->width*frame->height*(seed_depth+1);
		int p1;
		for(p1=0,p=start_seed_mask; p<end_seed_mask; p++,p1++)
				if(seed_array[p1]>0 && out[p]!=OUTBOUND)
				out[p] = MINIM;
	}		
		
}
	

  // Build pixel value lists again but now with only INIT pixels in each list

{ int p, v;

    for (v = 0; v <= maxval; v++)
			index[v] = -1;

		if (value8 != NULL)
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value8[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}

		else // value16 != NULL
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value16[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}
}

  // Find connected components of MINIM pixels and label with a unique catchment basin #.
{ int v, p;
    int qtop, qbot;

		nbasins = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == MINIM)
		{ PUSH(p);
			out[p] = ++nbasins;
		while (qbot >= 0)
		{ int q, *b, j, r;

		POP(q);
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ r = q+neighbor6[j];
		if (out[r] == MINIM)
		{ PUSH(r);
			out[r] = nbasins;
		}
		}
		}
		}
}

	int * seed_label_ref = (int *) calloc(nbasins+1, sizeof(int));
	uint8 * is_seeded_basin = (uint8 *) calloc(nbasins+1, sizeof(uint8));
{
	int i;
	for(i=0; i<=nbasins; i++){
		seed_label_ref[i] = i;
		is_seeded_basin[i]=0;
	}
	
	uint8 * seed_array = seed_image->array;
	int start_seed_mask = frame->width*frame->height*seed_depth;
	int end_seed_mask = frame->width*frame->height*(seed_depth+1);
	int p,p1;
	for(p1=0,p=start_seed_mask; p<end_seed_mask; p++,p1++)
		if(seed_array[p1]>0){
			if(out[p]>0 && out[p]<=nbasins)
				is_seeded_basin[out[p]] = 1;
		}

}

{ int v, p;

    for (v = 0; v < maxval; v++)
		{ int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

						qbot = -1;
				for (p = index[v]; p >= 0; p = t)
				{ int *b, j;

				t      = chord[p];
				if(out[p]==OUTBOUND)
					continue;
				out[p] = MASK;
				b      = boundary_pixels_3d(p);
				for (j = 0; j < n_nbrs6; j++)
						if (b[j] && out[p+neighbor6[j]] >= 0)
				{ out[p] = ONQUEUE;
				PUSH(p);
				break;
				}
				}

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

				while (qbot >= 0){ 
					int p, *b, j, q, o;

					POP(p);
					b = boundary_pixels_3d(p);
	
					for (j = 0; j < n_nbrs6; j++)
						if (b[j]){ 
							q = p + neighbor6[j];
							o = out[q];
							if(o==OUTBOUND)
									continue;
							if (o == MASK){ 
								out[q] = ONQUEUE;
								PUSH(q);
							}
							else 
								if (o != ONQUEUE && o != INIT){ 
									if (out[p] == ONQUEUE)
										out[p] = o;
									else
										if (o != out[p] && o != WSHED)
											if (out[p] == WSHED)
												out[p] = o;
											else{
												if(o<=0 || out[p]<=0)
													out[p] = WSHED;
												else{
													int adam_o = o;
													while(seed_label_ref[adam_o]!=adam_o)
															adam_o = seed_label_ref[adam_o];
													seed_label_ref[o] = adam_o;
													int adam_out_p = out[p];
													while(seed_label_ref[adam_out_p]!=adam_out_p)
															adam_out_p = seed_label_ref[adam_out_p];
													seed_label_ref[out[p]] = adam_out_p;
													
													if(is_seeded_basin[adam_o]==1 && is_seeded_basin[adam_out_p]==1) // both o and out[p] are seeded
														out[p] = WSHED;
													if(is_seeded_basin[adam_o]!=1 && is_seeded_basin[adam_out_p]!=1){ //neither o nor out[p] is seeded
														if(adam_out_p<=adam_o){
															seed_label_ref[out[p]] = adam_o;
															seed_label_ref[adam_out_p] = adam_o;
														}
														else{
															seed_label_ref[o] = adam_out_p;
															seed_label_ref[adam_o] = adam_out_p;
														}
													}
													if(is_seeded_basin[adam_o]==1 && is_seeded_basin[adam_out_p]!=1){ //o is seeded & out[p] is not
														// out[p]'s basin is going to be merged with o's
														seed_label_ref[out[p]] = adam_o;
														seed_label_ref[adam_out_p] = adam_o;
														is_seeded_basin[out[p]] = 1;
														is_seeded_basin[adam_out_p] = 1;
														out[p] = o;
													}
													if(is_seeded_basin[adam_o]!=1 && is_seeded_basin[adam_out_p]==1){ //o is not seeded & out[p] is
														// o's basin is going to be merged with out[p]'s
														seed_label_ref[o] = adam_out_p;
														seed_label_ref[adam_o] = adam_out_p;
														is_seeded_basin[o] = 1;
														is_seeded_basin[adam_o] = 1;
														out[p] = o;
													}
												}
											}
								}
						}
				}
		}
}

{
	int p;
	for(p=0; p<cvolume; p++){
		if(out[p]>0 && out[p]<=nbasins){
			int a = out[p];
			while(seed_label_ref[a]!=a)
				a = seed_label_ref[a];
			seed_label_ref[out[p]]=a;			
			out[p] = a;
		}
	}
}

	free(seed_label_ref);
	free(is_seeded_basin);
	 		

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

{ int p, qbot, qtop;

    qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == WSHED)
		{ int *b, j;

		b = boundary_pixels_3d(p);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j] && out[p+neighbor6[j]] > 0)
		{ PUSH(p);
		out[p] = ONQUEUE;
		break;
		}
		}

		while (qbot >= 0)
		{ int *b, j, q, o;

		POP(p);
		b = boundary_pixels_3d(p);

		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ q = p + neighbor6[j];
		o = out[q];
		if (o > 0)
				out[p] = o;
		else if (out[q] != ONQUEUE)
		{ PUSH(q);
		out[q] = ONQUEUE;
		}
		}
		}
}

{
	int * label_new_array = (int *) label_new->array;
	int p;
	for(p=0; p<cvolume; p++)
		label_new_array[p] = BODY_ID(p)!=mask_label?0:out[p];
}

/*	
  //  Make the data structure now

			result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
	result->nbasins = nbasins;
	result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

{ int   qbot, qtop;
    int   olabel;
		int   adjacent[256];
		int   alargest;
		int   p;

		seeds   = result->seeds;
		nlabels = 0;

		for (p = 1; p < 256; p++)
				adjacent[p] = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] > 0)

		{ PUSH(p);
		olabel   = out[p];
		out[p]   = 0;
		alargest = 0;

		while (qbot >= 0)
		{ int *b, j, q, r, o;

		POP(q);
		out[q] = -olabel;
      
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ r = q + neighbor6[j];
		o = out[r];
		if (o == olabel)
		{ PUSH(r);
		out[r] = 0;
		}
		else
		{ o = -o;
		if (o > 0 && o != olabel)
		{ o = seeds[o];
		adjacent[o] = 1;
		if (o > alargest)
				alargest = o;
		}
		}
		}
		}

		{ int j;

		for (j = 1; j <= alargest; j++)
			if (! adjacent[j])
				break;

		seeds[olabel] = j;
		if (j > nlabels)
		{ nlabels = j;
		if (j > 255)
		{ fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
		fprintf(stderr,"         Returning NULL result\n");
		fflush(stderr);
		return (NULL);
		}
		}

		for (j = 1; j <= alargest; j++)
			adjacent[j] = 0;
		}
		}
}


  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

{ int    p, b;
    uint8 *iarray;

		result->nlabels = nlabels;
		result->labels  = label_new; //Make_Stack(GREY,cwidth,cheight,cdepth);
		iarray = result->labels->array;

		for (p = 0; p < cvolume; p++)
		{ out[p] = b = -out[p];
		iarray[p] = seeds[b];
		}
		for (p = 0; p < cvolume; p++)
			if (seeds[b = out[p]]  > 0) 
				seeds[b] = -(p+1); 
		for (b = 1; b <= nbasins; b++)
						seeds[b] = -seeds[b];
}
*/

/*{ int    p;
    int   *color;
		Stack *view;

		view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

		for (p = 0; p < cvolume; p++)
		{ if (out[p] <= 0)
		{ printf("Unassigned pix\n");
		continue;
		}
		color = palette[(out[p]-1)%8];
		view->array[3*p]   = color[0];
		view->array[3*p+1] = color[1];
		view->array[3*p+2] = color[2];
		}

		Write_Stack("watershed.tif",view);

		Kill_Stack(view);
}
*/

{ int     q, p, i, j, v;
  int    *color, carea;
	Stack *view;
	uint8  *value8;
	uint16 *value16;

	view  = Make_Stack(COLOR,frame->width,frame->height,frame->depth);
	carea = frame->width * frame->height * frame->depth;

	if (frame->kind == GREY)
			value8 = frame->array;
	else
			value16 = (unsigned short *) (frame->array);

	q = 0;
	for (p = 0; p < carea; p++, q+=3)
	{ 
		if (out[p] <= 0)
		{ printf("Unassigned pix\n");
		view->array[q]   = 0;
		view->array[q+1] = 0;
		view->array[q+2] = 0;
		continue;
		}
	color = palette[(out[p]-1)%8];

	if (value8 == NULL)
			v = 0xFF - (value16[p] >> 8);
	else
			v = 0xFF - value8[p];

	v = (v >> 2) * 3;
	view->array[q]   = BODY_ID(p)!=mask_label?0:((color[0] >> 2) + v);
	view->array[q+1] = BODY_ID(p)!=mask_label?0:((color[1] >> 2) + v);
	view->array[q+2] = BODY_ID(p)!=mask_label?0:((color[2] >> 2) + v);
	}

	Write_Stack("watershed.tif",view);

	Kill_Stack(view);
}

  return (result);
}



/**************************************************************************************\
*                                                                                      *
		*  BUILD A SEEDED 3D WATERSHED DECOMPOSITION SO THAT IT IS CONFINED TO A USER DEFINED  *
		*	 SUB-VOLUME MASK (POSSIBLY IRREGULAR)																								 *		
		*                                                                                      *
		\**************************************************************************************/
int *Build_3D_Watershed_PreMasked_Seeded(Stack *frame, int iscon6, Image * seed_image, int seed_depth)
{ Watershed_3D *result;

	int    nlabels;   // The number of labels (1..nlabels)
			int    nbasins;   // Number of watersheds
			int   *seeds;     // A pixel in each watershed from which one can flood fill or contour

			int index[0x10001];
	int maxval;

	//(P)%plane_area
	
			int plane_area = frame->width * frame->height;
	int frame_width = frame->width;
	int frame_height = frame->height;
	/*
//#define PLANE_MAPPING(Q) (image_stack_scale*(original_width*((Q)/frame_width) + ((Q)%frame_width)))
//#define BODY_ID(P) (pmap[(int)(superpixel_2_seg_maps[(P)/plane_area][superpixel_maps[(P)/plane_area][PLANE_MAPPING((P)%plane_area)]])])	
//#define IS_BODY_VOXEL(P) (BODY_ID(P)==mask_label || BODY_ID(P)==0)
//#define IS_NOT_BODY_VOXEL(P) (BODY_ID(P)!=mask_label && BODY_ID(P)!=0)
			*/
	
			uint8  *value8;
	uint16 *value16;

	int   *chord;
	int   *out;

	if (frame->kind == COLOR)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to COLOR images\n");
	exit (1);
	}
	if (frame->kind == FLOAT32)
	{ fprintf(stderr,"Build_3D_Watershed: Cannot apply to FLOAT32 images\n");
	exit (1);
	}

	setup_globals_3D(frame->width,frame->height,frame->depth,iscon6);
	/*
	{
	FILE * fo = fopen("debug.txt", "w");
	int p;
	for(p=0; p<10; p++)
	fprintf(fo, "%u ", pmap[(int)(superpixel_2_seg_maps[p][superpixel_maps[p][200]])]);
	fclose(fo);
}
	*/
	
	chord = get_chord(2*cvolume,"Build_3D_Watershed");
	out   = chord + cvolume;

	if(seed_image==NULL)
			return out;
	
	if (frame->kind == GREY16)
	{ maxval  = 0x10000;
	value8  = NULL;
	value16 = (uint16 *) (frame->array);
	}
	else
	{ maxval  = 0x100;
	value8  = frame->array;
	value16 = NULL;
	}
	
  // Establish lists of pixels with each value

/*	{ int v, p;
		unsigned char * mask_array = mask->array;
	for (p = 0; p < cvolume; p++)
		if(mask_array[p]==)
			out[p] = OUTBOUND;
	else
					out[p] = MINIM;
	}
*/	
	{ int p, v;

	for (v = 0; v <= maxval; v++){
		index[v] = -1;
	}

	if (value8 != NULL)
		for (p = 0; p < cvolume; p++)
	{ //if(IS_NOT_BODY_VOXEL(p))
		//  continue;
		if(out[p]==OUTBOUND)
				continue;
		v        = value8[p];
		chord[p] = index[v];
		index[v] = p;
	}
	else  // value16 != NULL
				for (p = 0; p < cvolume; p++)
	{// if(IS_NOT_BODY_VOXEL(p))
		//	continue;
		if(out[p]==OUTBOUND)
				continue;
		v        = value16[p];
		chord[p] = index[v];
		index[v] = p;
	}
	}
	
#define PUSH(p)			\
{ if (qbot < 0)			\
			qtop = qbot = (p);		\
					else				\
					qtop = chord[qtop] = (p);	\
					chord[qtop] = -1;		\
}

#define POP(p)  { qbot = chord[p = qbot]; }

  //  Label all pixels that are not part of a catchment basin minimum with INIT
  //    and those that are with MINIM.

{ int v, p;

		
if (value8 != NULL)

	for (v = 0; v < maxval; v++)
{ int qtop, qbot, t;
  
          //  Push all c-pixels that are adjacent to a pixel of lesser value and mark INIT

							qbot = -1;
					for (p = index[v]; p >= 0; p = t)
					{ int *b, j;
					t = chord[p];
					//if(IS_NOT_BODY_VOXEL(p))
					//		continue;
					if(out[p]==OUTBOUND)
							continue;
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
							if (b[j] && value8[p+neighbor6[j]] < v)
					{ out[p] = INIT;
					PUSH(p)
							break;
					}
					}

          //  Find all c-pixels connected to another labeled INIT and mark as INIT
  
							while (qbot >= 0)
					{ int p, *b, j, q;
					int cdist;
  
					POP(p);
					b = boundary_pixels_3d(p);
					for (j = 0; j < n_nbrs6; j++)
							if (b[j])
					{ q = p + neighbor6[j];
					if (out[q] != INIT && out[q] != OUTBOUND && value8[q] == v)
					{ out[q] = INIT;
					PUSH(q);
					}
					}
					}
}
		
{	// Force connected components in the seed image to have connected components of MINIM labels
		uint8 * seed_array = seed_image->array;
		int start_seed_mask = frame->width*frame->height*seed_depth;
		int end_seed_mask = frame->width*frame->height*(seed_depth+1);
		int p1;
		for(p1=0,p=start_seed_mask; p<end_seed_mask; p++,p1++)
				if(seed_array[p1]>0 && out[p]!=OUTBOUND)
				out[p] = MINIM;
}		
		
}
	

  // Build pixel value lists again but now with only INIT pixels in each list

{ int p, v;

    for (v = 0; v <= maxval; v++)
			index[v] = -1;

		if (value8 != NULL)
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value8[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}

		else // value16 != NULL
		{ for (p = 0; p < cvolume; p++)
			if (out[p] == INIT)
			{ v        = value16[p];
			chord[p] = index[v];
			index[v] = p;
			}
		}
}

  // Find connected components of MINIM pixels and label with a unique catchment basin #.
{ int v, p;
    int qtop, qbot;

		nbasins = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == MINIM)
		{ PUSH(p);
		out[p] = ++nbasins;
		while (qbot >= 0)
		{ int q, *b, j, r;

		POP(q);
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ r = q+neighbor6[j];
		if (out[r] == MINIM)
		{ PUSH(r);
		out[r] = nbasins;
		}
		}
		}
		}
}

	int * seed_label_ref = (int *) calloc(nbasins+1, sizeof(int));
	uint8 * is_seeded_basin = (uint8 *) calloc(nbasins+1, sizeof(uint8));
{
	int i;
	for(i=0; i<=nbasins; i++){
		seed_label_ref[i] = i;
		is_seeded_basin[i]=0;
	}
	
	uint8 * seed_array = seed_image->array;
	int start_seed_mask = frame->width*frame->height*seed_depth;
	int end_seed_mask = frame->width*frame->height*(seed_depth+1);
	int p,p1;
	for(p1=0,p=start_seed_mask; p<end_seed_mask; p++,p1++)
			if(seed_array[p1]>0){
		if(out[p]>0 && out[p]<=nbasins)
			is_seeded_basin[out[p]] = 1;
			}

}

{ int v, p;

    for (v = 0; v < maxval; v++)
		{ int qbot, qtop, t;

        // Push all INIT c-pixels that are adjacent to a previously labeled pixel and
        //   label ONQUEUE, label all others MASK

						qbot = -1;
				for (p = index[v]; p >= 0; p = t)
				{ int *b, j;

				t      = chord[p];
				if(out[p]==OUTBOUND)
						continue;
				out[p] = MASK;
				b      = boundary_pixels_3d(p);
				for (j = 0; j < n_nbrs6; j++)
						if (b[j] && out[p+neighbor6[j]] >= 0)
				{ out[p] = ONQUEUE;
				PUSH(p);
				break;
				}
				}

        //  In BFS, reach all c-pixels (have MASK set) in order of distance, and label
        //    with catchment basin or WSHED

						while (qbot >= 0){ 
					int p, *b, j, q, o;

					POP(p);
					b = boundary_pixels_3d(p);
	
					for (j = 0; j < n_nbrs6; j++)
							if (b[j]){ 
						q = p + neighbor6[j];
						o = out[q];
						if(o==OUTBOUND)
								continue;
						if (o == MASK){ 
							out[q] = ONQUEUE;
							PUSH(q);
						}
						else 
							if (o != ONQUEUE && o != INIT){ 
							if (out[p] == ONQUEUE)
								out[p] = o;
							else
										if (o != out[p] && o != WSHED)
										if (out[p] == WSHED)
										out[p] = o;
							else{
								if(o<=0 || out[p]<=0)
									out[p] = WSHED;
								else{
									int adam_o = o;
									while(seed_label_ref[adam_o]!=adam_o)
											adam_o = seed_label_ref[adam_o];
									seed_label_ref[o] = adam_o;
									int adam_out_p = out[p];
									while(seed_label_ref[adam_out_p]!=adam_out_p)
											adam_out_p = seed_label_ref[adam_out_p];
									seed_label_ref[out[p]] = adam_out_p;
													
									if(is_seeded_basin[adam_o]==1 && is_seeded_basin[adam_out_p]==1) // both o and out[p] are seeded
											out[p] = WSHED;
									if(is_seeded_basin[adam_o]!=1 && is_seeded_basin[adam_out_p]!=1){ //neither o nor out[p] is seeded
										if(adam_out_p<=adam_o){
											seed_label_ref[out[p]] = adam_o;
											seed_label_ref[adam_out_p] = adam_o;
										}
										else{
											seed_label_ref[o] = adam_out_p;
											seed_label_ref[adam_o] = adam_out_p;
										}
									}
									if(is_seeded_basin[adam_o]==1 && is_seeded_basin[adam_out_p]!=1){ //o is seeded & out[p] is not
														// out[p]'s basin is going to be merged with o's
																seed_label_ref[out[p]] = adam_o;
														seed_label_ref[adam_out_p] = adam_o;
														is_seeded_basin[out[p]] = 1;
														is_seeded_basin[adam_out_p] = 1;
														out[p] = o;
									}
									if(is_seeded_basin[adam_o]!=1 && is_seeded_basin[adam_out_p]==1){ //o is not seeded & out[p] is
														// o's basin is going to be merged with out[p]'s
																seed_label_ref[o] = adam_out_p;
														seed_label_ref[adam_o] = adam_out_p;
														is_seeded_basin[o] = 1;
														is_seeded_basin[adam_o] = 1;
														out[p] = o;
									}
								}
							}
							}
							}
						}
		}
}

{
	int p;
	for(p=0; p<cvolume; p++){
		if(out[p]>0 && out[p]<=nbasins){
			int a = out[p];
			while(seed_label_ref[a]!=a)
					a = seed_label_ref[a];
			seed_label_ref[out[p]]=a;			
			out[p] = a;
		}
	}
}

	free(seed_label_ref);
	free(is_seeded_basin);
	 		

  //  Find all WSHED pixels in order of distance from a non-WSHED pixel and randomly
  //    assign to one of the nearest catchment basins

{ int p, qbot, qtop;

    qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] == WSHED)
		{ int *b, j;

		b = boundary_pixels_3d(p);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j] && out[p+neighbor6[j]] > 0)
		{ PUSH(p);
		out[p] = ONQUEUE;
		break;
		}
		}

		while (qbot >= 0)
		{ int *b, j, q, o;

		POP(p);
		b = boundary_pixels_3d(p);

		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
		{ q = p + neighbor6[j];
		o = out[q];
		if (o > 0)
				out[p] = o;
		else if (out[q] != ONQUEUE && out[q]!=OUTBOUND)
		{ PUSH(q);
		out[q] = ONQUEUE;
		}
		}
		}
}

/*{
	int * label_new_array = (int *) label_new->array;
	int p;
	for(p=0; p<cvolume; p++)
			label_new_array[p] = out[p];
//	for(p=0; p<cvolume; p++)
			//label_new_array[p] = BODY_ID(p)!=mask_label?0:out[p];
}
*/
/*	
  //  Make the data structure now

			result = new_watershed_3d(sizeof(int)*(nbasins+1),"Build_2D_Watershed");
	result->nbasins = nbasins;
	result->iscon6  = iscon6;

  //  Greedily assign labels to catchment basins (CBs) so that no basin is adjacent to
  //    another with the same number.  Typicall the number of labels required is very,
  //    very small, 6 or 7 in my experience.  At the end of this segment, out[p] is
  //    the negative of its label reassignment, seeds[b] is the new label for CB b,
  //    and nlabels is the number of unique labels used.

{ int   qbot, qtop;
    int   olabel;
		int   adjacent[256];
		int   alargest;
		int   p;

		seeds   = result->seeds;
		nlabels = 0;

		for (p = 1; p < 256; p++)
				adjacent[p] = 0;

		qbot = -1;
		for (p = 0; p < cvolume; p++)
				if (out[p] > 0)

{ PUSH(p);
		olabel   = out[p];
		out[p]   = 0;
		alargest = 0;

		while (qbot >= 0)
{ int *b, j, q, r, o;

		POP(q);
		out[q] = -olabel;
      
		b = boundary_pixels_3d(q);
		for (j = 0; j < n_nbrs6; j++)
				if (b[j])
{ r = q + neighbor6[j];
		o = out[r];
		if (o == olabel)
{ PUSH(r);
		out[r] = 0;
}
		else
{ o = -o;
		if (o > 0 && o != olabel)
{ o = seeds[o];
		adjacent[o] = 1;
		if (o > alargest)
				alargest = o;
}
}
}
}

{ int j;

		for (j = 1; j <= alargest; j++)
			if (! adjacent[j])
				break;

		seeds[olabel] = j;
		if (j > nlabels)
{ nlabels = j;
		if (j > 255)
{ fprintf(stderr,"Warning: More than 255 colors needed to label watershed\n");
		fprintf(stderr,"         Returning NULL result\n");
		fflush(stderr);
		return (NULL);
}
}

		for (j = 1; j <= alargest; j++)
			adjacent[j] = 0;
}
}
}


  //  Make an image with the assigned labels, and make seeds contain the lexicographically
  //    smalles pixel in each region.

{ int    p, b;
    uint8 *iarray;

		result->nlabels = nlabels;
		result->labels  = label_new; //Make_Stack(GREY,cwidth,cheight,cdepth);
		iarray = result->labels->array;

		for (p = 0; p < cvolume; p++)
{ out[p] = b = -out[p];
		iarray[p] = seeds[b];
}
		for (p = 0; p < cvolume; p++)
			if (seeds[b = out[p]]  > 0) 
				seeds[b] = -(p+1); 
		for (b = 1; b <= nbasins; b++)
						seeds[b] = -seeds[b];
}
*/

/*{ int    p;
    int   *color;
		Stack *view;

		view = Make_Stack(COLOR,frame->width,frame->height,frame->depth);

		for (p = 0; p < cvolume; p++)
{ if (out[p] <= 0)
{ printf("Unassigned pix\n");
		continue;
}
		color = palette[(out[p]-1)%8];
		view->array[3*p]   = color[0];
		view->array[3*p+1] = color[1];
		view->array[3*p+2] = color[2];
}

		Write_Stack("watershed.tif",view);

		Kill_Stack(view);
}*/

{ int     q, p, i, j, v;
  int    *color, carea;
	Stack *view;
	uint8  *value8;
	uint16 *value16;

	view  = Make_Stack(COLOR,frame->width,frame->height,frame->depth);
	carea = frame->width * frame->height * frame->depth;

	if (frame->kind == GREY)
			value8 = frame->array;
	else
			value16 = (unsigned short *) (frame->array);

	q = 0;
	for (p = 0; p < carea; p++, q+=3)
	{ 
/*		if (out[p] <= 0)
		{ printf("Unassigned pix\n");
		view->array[q]   = 0;
		view->array[q+1] = 0;
		view->array[q+2] = 0;
		continue;
		}
*/
		color = out[p]==OUTBOUND?0:palette[(out[p]-1)%8];

		if (value8 == NULL)
				v = 0xFF - (value16[p] >> 8);
		else
				v = 0xFF - value8[p];

		v = (v >> 2) * 3;
		view->array[q]   = out[p]==OUTBOUND?0:((color[0] >> 2) + v);
		view->array[q+1] = out[p]==OUTBOUND?0:((color[1] >> 2) + v);
		view->array[q+2] = out[p]==OUTBOUND?0:((color[2] >> 2) + v);
	}

	Write_Stack("watershed.tif",view);

	Kill_Stack(view);
}

  return out;
}
