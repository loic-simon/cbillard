// We use the casio graphics library,2018
// documentation: https://www.cairographics.org/manual/cairo-Paths.html
#include <unistd.h>
#include <assert.h>
#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>

typedef struct {                            // Type Particle
  double x, y, vx, vy, colR, colG, colB;      // Position, vitesse, couleur
  char plein, enJeu;                          // Boules pleines ou non, en jeu ou tombées
} Particle;

typedef struct {                            // Type Paroi
  double x1, y1, alpha, l;                    // Position, angle, longeur
} Paroi;

typedef struct {                            // Type Trou
  double xtrou, ytrou, radius;                // Postition, rayon
} Trou;


// ********* Code ci-dessous pratiquement non modifié (modification de quelques arguments publics de Graphics) *********

class Graphics {
 private:
  int Np;
  double lmin,lmax,Lmin,Lmax,Dimx,Dimy,L,l,d,diam;
  cairo_surface_t *sfc;
  cairo_t *cr;
  Display *dsp;
  Drawable da;
  int count;
 public:
  Graphics(int N, int Pix, double dmin, double dmax, double Dmin, double Dmax, double Lcc, double lcc, double dcc, double diameter){
    Np=N;
    count=1000;
    assert( dmin < dmax );
    assert(Np >0) ;
    lmin=dmin;
    lmax=dmax;
    Lmin=Dmin;
    Lmax=Dmax;
    Dimx=Pix*5/9;
    Dimy=Pix+6;
    L=Lcc;
    l=lcc;
    d=dcc;
    diam=diameter;

    if (( dsp = XOpenDisplay(NULL)) == NULL)      exit(1);//window management X11, and cairo graphics
    int screen = DefaultScreen(dsp);
    da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, Dimx, Dimy, 0, 0, 0);
    XMapWindow(dsp, da);
    sfc = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), Dimx , Dimy);
    cairo_xlib_surface_set_size(sfc, Dimx, Dimy);
    cr = cairo_create (sfc);
  } //empty window now on screen

  //some error messages on trying to duplicate the window
  Graphics & operator=(Graphics &g){// stop the program when copying the window
    printf("Don't use = with graphics objects  %p \n", &g ); exit(1);
  }
  Graphics(const Graphics &g ){// stop the program when passing windows as argument
    printf("Don't pass graphics objects: use a pointer  %p \n" , &g); exit(2);
  }
  void writePNG(){//save to png file, using "count" to produce a numbered file;
    char c[100];
    snprintf(c,99,"md%d.png", count++ );
    cairo_surface_write_to_png(sfc, c );
  }

  void draw(Paroi*,Trou*,Particle*,double,double,double,char);//draw the particles
  ~Graphics(){cairo_destroy (cr);cairo_surface_destroy (sfc); } //clean up function
};
