//  https:www.cairographics.org/manual/cairo-Paths.html fix
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "Graphics.h"

void
Graphics::draw(Paroi* par, Trou* T, Particle *p, double FPS, double alpha, double v0, char joueur) {

  struct timespec tm={0,1000000000/FPS};    // Sleep in nanoseconds between frames          (non modifié)
  XEvent event;                             // check if window closed and finish            (non modifié)
  if( XCheckWindowEvent(dsp,da, DestroyNotify , &event)){ XCloseDisplay(dsp); exit(1);}  // (non modifié)

  double gamma = Dimx/(lmax-lmin+4*d);      // Facteur d'échelle : gamma*X(m) = X(px) (largeur du jeu plus bords : 2*2d)
  double talon = gamma*2*d;                 // Taille du cadre du billard : doit être ajouté à tous les positionnements pour que les orignes en m et en px correspondent (équivalent de alpha dans le code initial)
  // char numParoi[10];                     // Utilisé pour debug : print numéro parois

  double xmin = talon + lmin*gamma;         // Coordonnées du billard (px)
  double xmax = talon + lmax*gamma;
  double ymin = talon + Lmin*gamma;
  double ymax = talon + Lmax*gamma;

  cairo_push_group(cr);                     // start drawing
  cairo_set_source_rgb(cr, 0, 0.19, 0.1);   // Vert foncé du tapis de la table
  cairo_paint (cr);                         // Fond vert foncé


  // Boules d'angle
  for(int i=16; i<28; i++){
    cairo_set_source_rgb(cr, 0,0.3,0);      // Vert clair des bandes
    cairo_new_sub_path(cr);
    cairo_arc(cr,(talon + gamma* (p[i].x -lmin)),Dimy-(talon + gamma*(p[i].y - lmin)),gamma*diam/2, 0, 2 * M_PI);
    cairo_fill(cr);
  }


  // Parois
  for(int i=0; i<18; i++){
    cairo_set_source_rgb(cr,0,0.3,0);
    cairo_set_line_width (cr, 1);
    cairo_move_to (cr,talon + gamma*par[i].x1,Dimy - talon - gamma*par[i].y1);
    // printf("%i, %lf\n", i, gamma*par[i].l*cos(par[i].alpha));
    cairo_rel_line_to (cr,gamma*par[i].l*cos(par[i].alpha), - gamma*par[i].l*sin(par[i].alpha));
    cairo_rel_line_to (cr,- gamma*d*sin(par[i].alpha), -gamma*d*cos(par[i].alpha));
    cairo_rel_line_to (cr,- gamma*par[i].l*cos(par[i].alpha), gamma*par[i].l*sin(par[i].alpha));
    cairo_rel_line_to (cr,gamma*d*sin(par[i].alpha), gamma*d*cos(par[i].alpha));
    cairo_fill(cr);

    // DEBUG : affichage du numéro des parois et de la ligne de contact avec le centre des particules

    //cairo_move_to (cr,talon + gamma*par[i].x1,Dimy - talon - gamma*par[i].y1);
    //// printf("%i, %lf\n", i, gamma*par[i].l*cos(par[i].alpha));
    //cairo_rel_line_to (cr,gamma*par[i].l*cos(par[i].alpha), -gamma*par[i].l*sin(par[i].alpha));

    //cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);//white
    //sprintf(numParoi, "%d",i);
    //cairo_set_font_size(cr,15);
    //cairo_show_text (cr, numParoi);

    //cairo_stroke (cr);

    //cairo_set_source_rgb(cr, 1,0,0);//white
    //cairo_move_to (cr,talon + gamma*(par[i].x1 + 0.5*diam*sin(par[i].alpha)),Dimy - talon - gamma*(par[i].y1 //- 0.5*diam*cos(par[i].alpha)));
    //// printf("%i, %lf\n", i, gamma*par[i].l*cos(par[i].alpha));
    //cairo_rel_line_to (cr,gamma*par[i].l*cos(par[i].alpha), -gamma*par[i].l*sin(par[i].alpha));

    //cairo_stroke (cr);

  }

  // Cadre du billard
  cairo_set_source_rgb (cr, 0.42, 0.2, 0.16);                               // Marron style bois
  cairo_rectangle (cr,0,0,xmax+2*talon,0.112*gamma);
  cairo_rectangle (cr,0,0,0.112*gamma,ymax+2*talon);
  cairo_rectangle (cr,xmax-0.012*gamma,0,0.112*gamma,ymax+2*talon);
  cairo_rectangle (cr,0,ymax-0.012*gamma+3,xmax+2*talon,0.112*gamma+6);

  // Cadre sous les trous
  cairo_rectangle (cr,0.1*gamma,0.1*gamma,20,20);
  cairo_rectangle (cr,0.1*gamma+xmax-40,0.1*gamma,20,20);
  cairo_rectangle (cr,0.1*gamma,0.1*gamma+ymax-40,20,20);
  cairo_rectangle (cr,0.1*gamma+xmax-40,0.1*gamma+ymax-40,20,20);
  cairo_fill(cr);
  cairo_stroke (cr);


  // Trous
  for(int i = 0; i < 6; i++){
    cairo_set_source_rgb(cr,0,0,0);         // Noir
    cairo_new_sub_path(cr);
    cairo_arc(cr,talon + T[i].xtrou * gamma, Dimy - talon - T[i].ytrou * gamma, gamma * T[i].radius, 0, 2 * M_PI);
    cairo_fill(cr);
  }


  // Boules réelles (de jeu)
  for(int i=0; i<16; i++){
    cairo_set_source_rgb(cr, p[i].colR, p[i].colG, p[i].colB);
    cairo_new_sub_path(cr);
    cairo_arc(cr,(talon + gamma* (p[i].x -lmin)),Dimy-(talon + gamma*(p[i].y - lmin)),gamma*diam/2, 0, 2 * M_PI);
    cairo_fill(cr);

    if (p[i].plein == 0) {    // Boules rayées : bandes blanches en haut et en bas de la boule
      cairo_set_source_rgb(cr,1,1,1);           // Blanc

      cairo_new_sub_path(cr);
      cairo_arc(cr,(talon + gamma* (p[i].x -lmin)),Dimy-(talon + gamma*(p[i].y - lmin)),gamma*diam/2, 0.6, M_PI-0.6);
      cairo_close_path(cr);
      cairo_fill(cr);

      cairo_new_sub_path(cr);
      cairo_arc(cr,(talon + gamma* (p[i].x -lmin)),Dimy-(talon + gamma*(p[i].y - lmin)),gamma*diam/2, M_PI+0.6, -0.6);
      cairo_close_path(cr);
      cairo_fill(cr);
    }
  }


  // Flèche pour le tir de la blanche
  if(v0 != 0) {         // Pour ne pas afficher de flèche, on envoie v0 = 0
    cairo_set_source_rgb(cr,0.4,0.4,0.4);       // Gris clair
    cairo_set_line_width(cr, 4);
    cairo_move_to(cr,talon + gamma* (p[0].x -lmin),Dimy-(talon + gamma*(p[0].y - lmin)));
    cairo_rel_line_to(cr,(10*v0)*cos(alpha), (10*v0)*(-sin(alpha)));
    cairo_stroke(cr);
  }

  // Affichage du joueur dont c'est le tour
  if (joueur == 1) {
    cairo_set_source_rgb(cr,1,1,1);
    cairo_move_to (cr, talon + gamma*(-3*d/2), Dimy - talon - gamma*(2*d + 20*d));
    cairo_set_font_size(cr,20);
    cairo_show_text (cr, "J1");

    cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
    cairo_new_sub_path(cr);
    cairo_arc(cr,talon + gamma*(-d), Dimy - talon - gamma*(2*d + 18*d),gamma*diam/2, 0, 2 * M_PI);
    cairo_fill(cr);

  } else if (joueur == -1) {
    cairo_set_source_rgb(cr,1,1,1);
    cairo_move_to (cr, talon + gamma*(lmax + d/2), Dimy - talon - gamma*(2*d + 20*d));
    cairo_set_font_size(cr,20);
    cairo_show_text (cr, "J2");

    cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
    cairo_new_sub_path(cr);
    cairo_arc(cr,talon + gamma*(lmax + d), Dimy - talon - gamma*(2*d + 18*d),gamma*diam/2, 0, 2 * M_PI);
    cairo_fill(cr);

    cairo_set_source_rgb(cr,1,1,1);           // Blanc
    cairo_new_sub_path(cr);
    cairo_arc(cr,talon + gamma*(lmax + d), Dimy - talon - gamma*(2*d + 18*d), gamma*diam/2, 0.6, M_PI-0.6);
    cairo_close_path(cr);
    cairo_fill(cr);

    cairo_new_sub_path(cr);
    cairo_arc(cr,talon + gamma*(lmax + d),Dimy - talon - gamma*(2*d + 18*d),gamma*diam/2, M_PI+0.6, -0.6);
    cairo_close_path(cr);
    cairo_fill(cr);
  }

  // ********** Fin du code modifié **********

  cairo_pop_group_to_source(cr);              // finished drawing operations for this set of positions
  cairo_paint(cr);                            // send to screen
  cairo_surface_flush(sfc);                   // send to x11
  nanosleep( &tm , NULL);                     // this sets the animations speed
  XFlush(dsp);                                // sync X11 to cairo

  //https://cairographics.org/Xlib/
}
