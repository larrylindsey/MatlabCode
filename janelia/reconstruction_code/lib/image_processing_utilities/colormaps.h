// Color maps for display label maps

#ifndef __COLORMAPS__
#define __COLORMAPS__
namespace iput
{
  namespace display
  {
    class ColorMap{
    public:
      ColorMap(int n_c, unsigned char * m){
        n_color = n_c;
        map = (unsigned char *) new unsigned char[n_c*3];
        for(int i=0; i<n_c*3; i++)
          map[i] = m[i];
      }
      ~ColorMap(){
        delete [] map;
      }
      int n_color;
      unsigned char * map;
    };

    unsigned char temp[26*3] =
    {
      179,  116,  116,
      179,  147,  116,
      179,  179,  116,
      147,  179,  116,
      116,  179,  116,
      116,  179,  147,
      116,  179,  179,
      116,  147,  179,
      116,  116,  179,
      147,  116,  179,
      179,  116,  179,
      179,  116,  147,
      179,  116,  116,
      179,   80,   80,
      179,  129,   80,
      179,  179,   80,
      129,  179,   80,
      80,  179,   80,
      80,  179,  129,
      80,  179,  179,
      80,  129,  179,
      80,   80,  179,
      129,   80,  179,
      179,   80,  179,
      179,   80,  129,
      179,   80,   80
    };
    ColorMap colormap_26(26, temp);

  }
}
#endif // __COLORMAPS__
