/*****************************************************************************************\
*                                                                                         *
*  Xfig drawing package                                                                   *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  August 2006                                                                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "xfig.h"

typedef struct
  { char     *name;
    FILE     *hfile;
    FILE     *tfile;
    Transform xform;
    int       colortop;

    int       pen_color;
    int       thickness;
    int       line_style;
    double    dash_width;
 
    int       join_style;
    int       cap_style;
 
    int       fill_color;
    int       fill_style;
 
    int       depth;

    int       point_size;
    int       font;
    int       flag;
    int       alignment;
 
    int       arrow_beg;
    int       arrow_end;
    int       arrow_style;
    int       arrow_fill;
    double    arrow_thick;
    double    arrow_width;
    double    arrow_height;
  } Myfig_Canvas;

static Myfig_Canvas *current_canvas;
static double        xscale;
static double        xoffset;
static double        yscale;
static double        yoffset;

static double        PI;


typedef struct __Point
  { struct __Point *next;
    Point           point;
  } _Point;

static _Point *Free_Point_List = NULL;
static int    Point_Offset, Point_Inuse;

static inline Point *new_point(char *routine)
{ _Point *object;

  if (Free_Point_List == NULL)
    { object = (_Point *) Guarded_Malloc(sizeof(_Point),routine);
      Point_Offset = ((char *) &(object->point)) - ((char *) object);
    }
  else
    { object = Free_Point_List;
      Free_Point_List = object->next;
    }
  Point_Inuse += 1;
  return (&(object->point));
}

static inline Point *copy_point(Point *point)
{ Point *copy = new_point("Copy_Point");
  *copy = *point;
  return (copy);
}

Point *Copy_Point(Point *point)
{ return (copy_point(point)); }

static inline void free_point(Point *point)
{ _Point *object  = (_Point *) (((char *) point) - Point_Offset);
  object->next = Free_Point_List;
  Free_Point_List = object;
  Point_Inuse -= 1;
}

void Free_Point(Point *point)
{ free_point(point); }

static inline void kill_point(Point *point)
{
  free(((char *) point) - Point_Offset);
  Point_Inuse -= 1;
}

void Kill_Point(Point *point)
{ kill_point(point); }

static inline void reset_point()
{ _Point *object;
  while (Free_Point_List != NULL)
    { object = Free_Point_List;
      Free_Point_List = object->next;
      kill_point(&(object->point));
      Point_Inuse += 1;
    }
}

void Reset_Point()
{ reset_point(); }

int Point_Usage()
{ return (Point_Inuse); }

Point *Make_Point(double x, double y)
{ Point *point = new_point("Make_Point");
  point->x = x;
  point->y = y;
  return (point);
}


typedef struct __Box
  { struct __Box *next;
    Box           box;
  } _Box;

static _Box *Free_Box_List = NULL;
static int    Box_Offset, Box_Inuse;

static inline Box *new_box(char *routine)
{ _Box *object;

  if (Free_Box_List == NULL)
    { object = (_Box *) Guarded_Malloc(sizeof(_Box),routine);
      Box_Offset = ((char *) &(object->box)) - ((char *) object);
    }
  else
    { object = Free_Box_List;
      Free_Box_List = object->next;
    }
  Box_Inuse += 1;
  return (&(object->box));
}

static inline Box *copy_box(Box *box)
{ Box *copy = new_box("Copy_Box");
  *copy = *box;
  return (copy);
}

Box *Copy_Box(Box *box)
{ return (copy_box(box)); }

static inline void free_box(Box *box)
{ _Box *object  = (_Box *) (((char *) box) - Box_Offset);
  object->next = Free_Box_List;
  Free_Box_List = object;
  Box_Inuse -= 1;
}

void Free_Box(Box *box)
{ free_box(box); }

static inline void kill_box(Box *box)
{
  free(((char *) box) - Box_Offset);
  Box_Inuse -= 1;
}

void Kill_Box(Box *box)
{ kill_box(box); }

static inline void reset_box()
{ _Box *object;
  while (Free_Box_List != NULL)
    { object = Free_Box_List;
      Free_Box_List = object->next;
      kill_box(&(object->box));
      Box_Inuse += 1;
    }
}

void Reset_Box()
{ reset_box(); }

int Box_Usage()
{ return (Box_Inuse); }

Box *Make_Box(double w, double h)
{ Box *box = new_box("Make_Box");
  box->w = w;
  box->h = h;
  return (box);
}


typedef struct __Frame
  { struct __Frame *next;
    Frame           frame;
  } _Frame;

static _Frame *Free_Frame_List = NULL;
static int    Frame_Offset, Frame_Inuse;

static inline Frame *new_frame(char *routine)
{ _Frame *object;

  if (Free_Frame_List == NULL)
    { object = (_Frame *) Guarded_Malloc(sizeof(_Frame),routine);
      Frame_Offset = ((char *) &(object->frame)) - ((char *) object);
    }
  else
    { object = Free_Frame_List;
      Free_Frame_List = object->next;
    }
  Frame_Inuse += 1;
  return (&(object->frame));
}

static inline Frame *copy_frame(Frame *frame)
{ Frame *copy = new_frame("Copy_Frame");
  *copy = *frame;
  return (copy);
}

Frame *Copy_Frame(Frame *frame)
{ return (copy_frame(frame)); }

static inline void free_frame(Frame *frame)
{ _Frame *object  = (_Frame *) (((char *) frame) - Frame_Offset);
  object->next = Free_Frame_List;
  Free_Frame_List = object;
  Frame_Inuse -= 1;
}

void Free_Frame(Frame *frame)
{ free_frame(frame); }

static inline void kill_frame(Frame *frame)
{
  free(((char *) frame) - Frame_Offset);
  Frame_Inuse -= 1;
}

void Kill_Frame(Frame *frame)
{ kill_frame(frame); }

static inline void reset_frame()
{ _Frame *object;
  while (Free_Frame_List != NULL)
    { object = Free_Frame_List;
      Free_Frame_List = object->next;
      kill_frame(&(object->frame));
      Frame_Inuse += 1;
    }
}

void Reset_Frame()
{ reset_frame(); }

int Frame_Usage()
{ return (Frame_Inuse); }

Frame *Make_Frame(Point *corner, Box *window)
{ Frame *frame = new_frame("Make_Frame");
  frame->corner = *corner;
  frame->window = *window;
  Free_Point(corner);
  Free_Box(window);
  return (frame);
}


typedef struct __Transform
  { struct __Transform *next;
    Transform           transform;
  } _Transform;

static _Transform *Free_Transform_List = NULL;
static int    Transform_Offset, Transform_Inuse;

static inline Transform *new_transform(char *routine)
{ _Transform *object;

  if (Free_Transform_List == NULL)
    { object = (_Transform *) Guarded_Malloc(sizeof(_Transform),routine);
      Transform_Offset = ((char *) &(object->transform)) - ((char *) object);
    }
  else
    { object = Free_Transform_List;
      Free_Transform_List = object->next;
    }
  Transform_Inuse += 1;
  return (&(object->transform));
}

static inline Transform *copy_transform(Transform *transform)
{ Transform *copy = new_transform("Copy_Transform");
  *copy = *transform;
  return (copy);
}

Transform *Copy_Transform(Transform *transform)
{ return (copy_transform(transform)); }

static inline void free_transform(Transform *transform)
{ _Transform *object  = (_Transform *) (((char *) transform) - Transform_Offset);
  object->next = Free_Transform_List;
  Free_Transform_List = object;
  Transform_Inuse -= 1;
}

void Free_Transform(Transform *transform)
{ free_transform(transform); }

static inline void kill_transform(Transform *transform)
{
  free(((char *) transform) - Transform_Offset);
  Transform_Inuse -= 1;
}

void Kill_Transform(Transform *transform)
{ kill_transform(transform); }

static inline void reset_transform()
{ _Transform *object;
  while (Free_Transform_List != NULL)
    { object = Free_Transform_List;
      Free_Transform_List = object->next;
      kill_transform(&(object->transform));
      Transform_Inuse += 1;
    }
}

void Reset_Transform()
{ reset_transform(); }

int Transform_Usage()
{ return (Transform_Inuse); }

Transform *Make_Transform(Frame *from, Frame *to)
{ Transform *xform = new_transform("Make_Transform");
  double xscale, yscale;

  xform->xscale  = xscale = (to->window.w) / (from->window.w);
  xform->yscale  = yscale = (to->window.h) / (from->window.h);
  xform->xoffset = to->corner.x - from->corner.x * xscale;
  xform->yoffset = to->corner.y - from->corner.y * yscale;
  Free_Frame(from);
  Free_Frame(to);
  return (xform);
}

void Finish_Xfig_Drawing(Xfig_Canvas *canvas)
{ Myfig_Canvas *this_canvas = (Myfig_Canvas *) canvas;
  char         *command;
  char         *name;

  fclose(this_canvas->hfile);
  fclose(this_canvas->tfile);
  name    = this_canvas->name;
  command = (char *) Guarded_Malloc(strlen(name)*3 + 100,"Finish_Xfig_Drawing");
  sprintf(command,"cat .%s.head .%s.tail >%s.fig",name,name,name);
  system(command);
  sprintf(command,"rm -f .%s.head .%s.tail",name,name);
  system(command);
  if (current_canvas == this_canvas)
    current_canvas = NULL;
  free(name);
  free(this_canvas);
}

Xfig_Canvas *Begin_Xfig_Drawing(char *name, Transform *xform)
{ Myfig_Canvas *canvas;
  FILE         *hfile;

  canvas = (Myfig_Canvas *) Guarded_Malloc(sizeof(Myfig_Canvas),"Begin_Xfig_Drawing");

  { char *buffer = (char *) Guarded_Malloc(strlen(name)+7,"Begin_Xfig_Drawing"); 

    sprintf(buffer,".%s.head",name);
    canvas->hfile = hfile = Guarded_Fopen(buffer,"w","Begin_Xfig_Drawing"); 
    sprintf(buffer,".%s.tail",name);
    canvas->tfile = Guarded_Fopen(buffer,"w","Begin_Xfig_Drawing"); 
    free(buffer);
  }

  PI = acos(-1.);

  xform->xscale  *= 1200;   //  Don't forget the target space is in inches!  1200pts/inch
  xform->xoffset *= 1200;
  xform->yscale  *= 1200;
  xform->yoffset *= 1200;

  canvas->name     = Guarded_Strdup(name,"Begin_Xfig_Drawing");
  canvas->xform    = *xform;
  canvas->colortop = 32;

  current_canvas  = canvas;
  xscale          = xform->xscale;
  xoffset         = xform->xoffset;
  yscale          = xform->yscale;
  yoffset         = xform->yoffset;

  fprintf(hfile,"#FIG 3.2\n");
  fprintf(hfile,"Landscape\n");
  fprintf(hfile,"Center\n");                                                          
  fprintf(hfile,"Inches\n");                                                          
  fprintf(hfile,"Letter\n");                                                          
  fprintf(hfile,"100.\n");
  fprintf(hfile,"Single\n");
  fprintf(hfile,"-2\n");
  fprintf(hfile,"1200 2\n");

  canvas->pen_color  = BLACK;
  canvas->thickness  = 1;
  canvas->line_style = SOLID;
  canvas->dash_width = 3.;

  canvas->join_style = POINTED_JOIN;
  canvas->cap_style  = SQUARED_CAP;

  canvas->fill_style = -1;
  canvas->fill_color = WHITE;

  canvas->font       = TIMES;
  canvas->flag       = 4;
  canvas->point_size = 12;
  canvas->alignment  = LEFT_ALIGN;

  canvas->depth = 50;

  canvas->arrow_beg    = 0;
  canvas->arrow_end    = 0;
  canvas->arrow_style  = LINE_STYLE;
  canvas->arrow_fill   = 0;
  canvas->arrow_thick  = 1.;
  canvas->arrow_width  = 120.;
  canvas->arrow_height = 60.;

  return ((Xfig_Canvas *) canvas);
}

void Set_Current_Canvas(Xfig_Canvas *canvas)
{ current_canvas = (Myfig_Canvas *) canvas;
  xscale          = current_canvas->xform.xscale;
  xoffset         = current_canvas->xform.xoffset;
  yscale          = current_canvas->xform.yscale;
  yoffset         = current_canvas->xform.yoffset;
}

int Add_Color_Xfig_Canvas(double red, double green, double blue)
{ int ired, igreen, iblue, hex;

  ired   = red * 255;
  igreen = green * 255;
  iblue  = blue * 255;
  hex    = (ired << 16) | (igreen << 8) | iblue;
  fprintf(current_canvas->hfile,"0 %d #.6x\n",current_canvas->colortop,hex);
  return (current_canvas->colortop++);
}

void Set_Pen_Color(int color)
{ current_canvas->pen_color = color; }

void Set_Thickness(int thick)
{ current_canvas->thickness = thick; }

void Set_Line_Style(int line_style, double dash_width)
{ current_canvas->line_style = line_style;
  if (dash_width > 0.)
    current_canvas->dash_width = dash_width;
}

void Set_Join_Style(int join_style)
{ current_canvas->join_style = join_style; }

void Set_Cap_Style(int cap_style)
{ current_canvas->cap_style = cap_style; }

void Set_Fill_Style(int fill_style)
{ current_canvas->fill_style = fill_style; }

void Set_Fill_Color(int color)
{ current_canvas->fill_color = color; }

void Set_Depth(int depth)
{ current_canvas->depth = depth; }

void Set_Arrow_Mode(int beg, int end)
{ current_canvas->arrow_beg = beg;
  current_canvas->arrow_end = end;
}

void Set_Arrow_Style(int a_style, int fill)
{ current_canvas->arrow_style = a_style;
  current_canvas->arrow_fill  = fill;
}

void Set_Arrow_Sizes(double thick, double width, double height)
{ current_canvas->arrow_thick  = thick;
  current_canvas->arrow_width  = 15.*width;
  current_canvas->arrow_height = 15.*width;
}

static Point *Point_Array = NULL;
static int    Point_Top;
static int    Point_Max = 0;

void First_Line_Point(Point *p)
{ if (Point_Max == 0)
    { Point_Max = 100;
      Point_Array = (Point *) Guarded_Malloc(sizeof(Point)*Point_Max,"First_Line_Point");
    }
  Point_Top = 0;
  Point_Array[Point_Top++] = *p;
}

void Next_Line_Point(Point *p)
{ if (Point_Top >= Point_Max)
    { Point_Max = 1.2*Point_Top + 50;
      Point_Array = (Point *)
                      Guarded_Realloc(Point_Array,sizeof(Point)*Point_Max,"Next_Line_Point");
    }
  Point_Array[Point_Top++] = *p;
}

static inline int mapx(double x)
{ return ((int) (x*xscale + xoffset)); } 

static inline int mapy(double y)
{ return ((int) (y*yscale + yoffset)); } 

static inline int scalex(double x)
{ return ((int) (x*xscale)); } 

static inline int scaley(double y)
{ return ((int) (y*yscale)); } 

static void emit_line(int type, int roundness)
{ Point *p;

  fprintf(current_canvas->tfile,"2 %d %d %d %d %d %d -1 %d %f %d %d %d %d %d %d\n",
          type,
          current_canvas->line_style,
          current_canvas->thickness,
          current_canvas->pen_color,
          current_canvas->fill_color,
          current_canvas->depth,
          current_canvas->fill_style,
          current_canvas->dash_width,
          current_canvas->join_style,
          current_canvas->cap_style,
          roundness,
          current_canvas->arrow_end,
          current_canvas->arrow_beg,
          Point_Top);
  if (current_canvas->arrow_end)
    fprintf(current_canvas->tfile,"\t%d %d %f %f %f\n",
            current_canvas->arrow_style,
            current_canvas->arrow_fill,
            current_canvas->arrow_thick,
            current_canvas->arrow_width,
            current_canvas->arrow_height);
  if (current_canvas->arrow_beg)
    fprintf(current_canvas->tfile,"\t%d %d %f %f %f\n",
            current_canvas->arrow_style,
            current_canvas->arrow_fill,
            current_canvas->arrow_thick,
            current_canvas->arrow_width,
            current_canvas->arrow_height);
  fprintf(current_canvas->tfile,"\t");
  for (p = Point_Array; p < Point_Array+Point_Top; p++)
    fprintf(current_canvas->tfile," %d %d",mapx(p->x),mapy(p->y));
  fprintf(current_canvas->tfile,"\n");
}

void Last_Line_Point(Point *p, int closed)
{ Next_Line_Point(p);
  if (closed)
    { Next_Line_Point(Point_Array);
      emit_line(3,-1);
    }
  else
    emit_line(1,-1);
}

void Etch_Rectangle(Point *ul_corner, Point *lr_corner, int roundness)
{ Point corner;

  First_Line_Point(ul_corner);
  corner.x = ul_corner->x;
  corner.y = lr_corner->y;
  Next_Line_Point(&corner);
  Next_Line_Point(lr_corner);
  corner.x = lr_corner->x;
  corner.y = ul_corner->y;
  Next_Line_Point(&corner);
  Next_Line_Point(ul_corner);
  if (roundness > 0)
    emit_line(4,roundness);
  else
    emit_line(2,-1);
}

void Etch_Ellipse(Point *center, double xradius, double yradius, double angle)
{ fprintf(current_canvas->tfile,"1 1 %d %d %d %d %d -1 %d 0.0 1 %f %d %d %d %d %d %d %d %d\n",
          current_canvas->line_style,
          current_canvas->thickness,
          current_canvas->pen_color,
          current_canvas->fill_color,
          current_canvas->depth,
          current_canvas->fill_style,
          (PI / 360.) * angle,
          mapx(center->x),
          mapy(center->y),
          scalex(xradius),
          scaley(yradius),
          mapx(center->x),
          mapy(center->y),
          mapx(center->x + xradius),
          mapy(center->y + yradius));
}

void Etch_Circle(Point *center, double radius)
{ Etch_Ellipse(center,radius,radius,0.); }

void Set_Point_Size(int size)
{ current_canvas->point_size = size; }

void Set_Font(int font, int style)
{ if (font <= PALATINO)
    { current_canvas->font = font + style;
      current_canvas->flag = 4;
    }
  else if (font >= LATEX_ROMAN)
    { current_canvas->font = (font - LATEX_ROMAN);
      if (style == ITALIC)
        current_canvas->font += 2;
      else if (style == BOLD)
        current_canvas->font += 1;
      current_canvas->flag = 0;
    }
  else
    { current_canvas->font = font;
      current_canvas->flag = 4;
    }
}

void Set_Alignment(int align)
{ current_canvas->alignment = align; }

void Etch_Text(Point *anchor, double angle, char *text)
{ fprintf(current_canvas->tfile,"4 %d %d %d -1 %d %d %f %d 0 0 %d %d %s\\001\n",
          current_canvas->alignment,
          current_canvas->pen_color,
          current_canvas->depth,
          current_canvas->font,
          current_canvas->point_size,
          (PI / 360.) * angle,
          current_canvas->flag,
          mapx(anchor->x),
          mapy(anchor->y),
          text);
}
