// Bibliothèques standard
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>       // Pour initialiser le germe aléatoire

// Pour la lecture des touches du clavier
#include <termios.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/uio.h>

#include "Graphics.h"     // Définition des types Particle, Paroi, Trou et de la classe Graphics

#define EPS 1e-10         // Petite erreur, permet d'éviter les erreurs d'arrondis : utilisés pour les comparaisons d'un double à 0


enum col_type {           // Types de collision possible
    particle,
    paroi
};

typedef struct {           // Type Event : collisions possibles boule-boule et boule-paroi
    enum col_type type;         // type de la collision
    int i, j;                   // particule 1, particule 2 / paroi
    double time;                // temps avant avènement de l'évènement (en considérant la vitesse constante, c.f. rapport)
} Event;


double rad(double deg){                           //Conversion degrés -> radians
    return M_PI*deg/180;
}

void initparois(Paroi* par, double L, double l, double d) {
    // Dimensions du jeu, prises de https://www.ffbillard.com/pages/le-materiel-25.html :
    double dHorExt = 0.15,  dHorInt = 0.11;                           // Largeurs extérieure et intérieure des poches latérales
    double dCoinExt = 0.14, dCoinInt = 0.105;                         // Largeurs extérieure et intérieure des poches d'angle
    double alphaCoin = rad(30);                                       // Angle par rapport à la paroi des poches d'angle
    double lCoin = (dCoinExt - dCoinInt)/(2*sin(M_PI/4 - alphaCoin)); // Longueurs des bandes encadrant les poches d'angle
    double alphaHor = rad(30);                                        // Angle par rapport à la paroi des poches latérales
    double lHor = (dHorExt - dHorInt)/(2*sin(alphaHor));              // Longueurs des bandes encadrant les poches latérales

    // HAUT

    par[0].x1 = d + dCoinExt/sqrt(2);
    par[0].y1 = L + d;
    par[0].alpha = - alphaCoin;
    par[0].l = - lCoin;

    par[1].x1 = d + dCoinExt/sqrt(2);
    par[1].y1 = L + d;
    par[1].alpha = 0;
    par[1].l = l - sqrt(2)*dCoinExt;

    par[2].x1 = d + l - dCoinExt/sqrt(2);
    par[2].y1 = L + d;
    par[2].alpha = alphaCoin;
    par[2].l = lCoin;

    // BAS

    par[3].x1 = d + dCoinExt/sqrt(2);
    par[3].y1 = d;
    par[3].alpha = - M_PI + alphaCoin;
    par[3].l = lCoin;

    par[4].x1 = d + dCoinExt/sqrt(2);
    par[4].y1 = d;
    par[4].alpha = M_PI;
    par[4].l = - (l - sqrt(2)*dCoinExt);

    par[5].x1 = d + l - dCoinExt/sqrt(2);
    par[5].y1 = d;
    par[5].alpha = M_PI - alphaCoin;
    par[5].l = - lCoin;

    // GAUGHE en bas

    par[6].x1 = d;
    par[6].y1 = d + dCoinExt/sqrt(2);
    par[6].alpha = - M_PI/2 - alphaCoin - M_PI;
    par[6].l = - lCoin;

    par[7].x1 = d;
    par[7].y1 = d + dCoinExt/sqrt(2);
    par[7].alpha = M_PI/2;
    par[7].l = L/2 - dHorExt/2 - dCoinExt/sqrt(2);

    par[8].x1 = d;
    par[8].y1 = d + L/2 - dHorExt/2;
    par[8].alpha = M_PI - alphaHor;
    par[8].l = lHor;

    // GAUGHE en haut

    par[9].x1 = d;
    par[9].y1 = d + L/2 + dHorExt/2;
    par[9].alpha = alphaHor;
    par[9].l = - lHor;

    par[10].x1 = d;
    par[10].y1 = d + L/2 + dHorExt/2;
    par[10].alpha =  M_PI/2;
    par[10].l = L/2 - dHorExt/2 - dCoinExt/sqrt(2);

    par[11].x1 = d;
    par[11].y1 = d + L - dCoinExt/sqrt(2);
    par[11].alpha = alphaCoin + M_PI/2;
    par[11].l = lCoin;

    // DROITE en bas

    par[12].x1 = l + d;
    par[12].y1 = d + dCoinExt/sqrt(2);
    par[12].alpha = - M_PI/2 + alphaCoin;
    par[12].l = lCoin;

    par[13].x1 = l + d;
    par[13].y1 = d + dCoinExt/sqrt(2);
    par[13].alpha = - M_PI/2;
    par[13].l = - (L/2 - dHorExt/2 - dCoinExt/sqrt(2));

    par[14].x1 = l + d;
    par[14].y1 = d + L/2 - dHorExt/2;
    par[14].alpha = alphaHor - M_PI;
    par[14].l = - lHor;

    // DROITE en haut

    par[15].x1 = l + d;
    par[15].y1 = d + L/2 + dHorExt/2;
    par[15].alpha = - alphaHor;
    par[15].l = lHor;

    par[16].x1 = l + d;
    par[16].y1 = d + L/2 + dHorExt/2;
    par[16].alpha = - M_PI/2;
    par[16].l = - (L/2 - dHorExt/2 - dCoinExt/sqrt(2));

    par[17].x1 = l + d;
    par[17].y1 = d + L - dCoinExt/sqrt(2);
    par[17].alpha = - M_PI/2 - alphaCoin;
    par[17].l = - lCoin;

}


void initParoiBalls(Particle* p, Paroi* par, double r){      // Initialisation des boules d'angle
    for (int i = 0; i < 12; i++) {              // Particules immobiles
        p[16 + i].vx = 0;
        p[16 + i].vy = 0;
    }

    p[16].x = par[0].x1;                        // Positions : au niveau du croisement des parois formant un angle, puis translatées à l'intérieur des bandes (empiriquement) pour qu'elles dépassent à peine
    p[16].y = par[0].y1 + r;
    p[17].x = par[2].x1;
    p[17].y = par[2].y1 + r;

    p[18].x = par[3].x1;
    p[18].y = par[3].y1 - r;
    p[19].x = par[5].x1;
    p[19].y = par[5].y1 - r;

    p[20].x = par[6].x1 - r;
    p[20].y = par[6].y1;
    p[21].x = par[8].x1 - r;
    p[21].y = par[8].y1 - 0.4*r;
    p[22].x = par[9].x1 - r;
    p[22].y = par[9].y1 + 0.4*r;
    p[23].x = par[11].x1 - r;
    p[23].y = par[11].y1;

    p[24].x = par[12].x1 + r;
    p[24].y = par[12].y1;
    p[25].x = par[14].x1 + r;
    p[25].y = par[14].y1 - 0.4*r;
    p[26].x = par[15].x1 + r;
    p[26].y = par[15].y1 + 0.4*r;
    p[27].x = par[17].x1 + r;
    p[27].y = par[17].y1;
}


void initparticles(Particle *p, double L, double l, double d, double diameter){   // Initialisation des particules
    //Définition des positions initiales des boules et de leurs couleurs
    //  Boule N°       0  1    2    3   4     5    6    7    8   9    10   11  12   13  14    15
    double xsrp[16] = {0, 0,   -1,  -2, 3,    4,   -1,  -2,  0,  1,   2,   -4, 1,   -3,  2,   0};
    double yss[16]  = {0, -8,  -5,  -2, 1,    4,   1,   4,   -2, -5,  -2,  4,  1,   1,   4,   4};
    double colR[16] = {1, 1,   0.1, 1,  0.6,  1,   0,   0.5, 0,  1,   0.1, 1,  0.6, 1,   0,   0.5};
    double colG[16] = {1, 0.7, 0.3, 0,  0.2,  0.5, 0.6, 0.1, 0,  0.7, 0.3, 0,  0.2, 0.5, 0.6, 0.1};
    double colB[16] = {1, 0.1, 0.7, 0,  0.6,  0.2, 0.3, 0.1, 0,  0.1, 0.7, 0,  0.6, 0.2, 0.3, 0.1};
    char plein[16]  = {1, 1,   1,   1,  1,    1,   1,   1,   1,  0,   0,   0,  0,   0,   0,   0};

    double rplus = 1.1*diameter/2;              // Moitié de la distance entre deux boules
    double s = rplus/sqrt(3);                   // Distance verticale entre les centres des boules de deux lignes consécutives
    double xc = d + l/2, yc = d + 3*L/4;        // Coordonnées du centre du triangle

    for (int i = 0; i < 16; i++) {              //Màj du tableau des boules
        if (i == 0) {
            p[i].x = d + l/2;
            p[i].y = L/4;
        } else {
            p[i].x = xc + xsrp[i]*rplus;
            p[i].y = yc + yss[i]*s;
        }
        p[i].vx = 0;
        p[i].vy = 0;
        p[i].colR = colR[i];
        p[i].colG = colG[i];
        p[i].colB = colB[i];
        p[i].plein = plein[i];
        p[i].enJeu = 1;
    }
}


void initTrou(Trou* T, Paroi* par) {              // Initialisation des trous
  double dHorInt = 0.11;     // Largeur intérieure des poches latérales
  double dCoinInt = 0.105;   // Largeur intérieure des poches d'angle

  int bord1[6] = {0, 8, 3, 5,  14, 2};    // Numéros des parois encadrant les trous 0 à 5
  int bord2[6] = {11,9, 6, 12, 15, 17};

  for (int i = 0;i<6;i++) {
    T[i].xtrou = (par[bord1[i]].x1 + par[bord1[i]].l*cos(par[bord1[i]].alpha) + par[bord2[i]].x1 + par[bord2[i]].l*cos(par[bord2[i]].alpha))/2;   //  Positions des trous : au milieu du segment reliant les extrémités des parois encadrant le trou
    T[i].ytrou = (par[bord1[i]].y1 + par[bord1[i]].l*sin(par[bord1[i]].alpha) + par[bord2[i]].y1 + par[bord2[i]].l*sin(par[bord2[i]].alpha))/2;
  }

  T[0].radius = dCoinInt/2;       // Rayon du trou
  T[1].radius = dHorInt/2;
  T[2].radius = dCoinInt/2;
  T[3].radius = dCoinInt/2;
  T[4].radius = dHorInt/2;
  T[5].radius = dCoinInt/2;
}


char colPart(Particle* p, double diameter, int i, int j, Event* e, long n) {    // Gère la collision entre deux particules
  //Calcul des coefficients de l'équation quadratique
  double dx = p[i].x - p[j].x;
  double dy = p[i].y - p[j].y;
  double dvx = p[i].vx - p[j].vx;
  double dvy = p[i].vy - p[j].vy;

  double a = pow(dvx,2) + pow(dvy,2);                   // Coefficients de l'équation
  double b = 2*dx*dvx + 2*dy*dvy;
  double c = pow(dx,2) + pow(dy,2) - pow(diameter,2);

  // Calcul des solutions de l'équation
  double delta = pow(b,2) - 4*a*c;
  if (delta >= 0) {
    double tplus = (-b + sqrt(delta))/(2*a);
    double tmoins = (-b - sqrt(delta))/(2*a);

    if (tmoins > EPS) {                   // Si tmoins est positif (à la précision près), il est nécessairement inférieur à tplus, donc c'est le seul qu'on ajoute à e
      e[n].type = particle;
      e[n].i = i;
      e[n].j = j;
      e[n].time = tmoins;
      return(1);                          // On ajoute 1 au nombre d'évènements
    } else if (tplus > EPS) {             // Si tmoins n'est pas positif, on teste tplus, qui est plus grand
      e[n].type = particle;
      e[n].i = i;
      e[n].j = j;
      e[n].time = tplus;
      return(1);
    }
  }
  return(0);                               // Si tplus et tmoins sont négatifs, il n'y a pas de collisions, on n'ajoute pas d'évènements
}


void actPosVit(Particle* p, int Np, double t) {   // Màj des positions et des vitesses pour un temps donné
    double V0,g = 9.81,mu = 0.3;
    for (int i = 0; i < Np; i++) {
        V0 = sqrt(pow(p[i].vx,2) + pow(p[i].vy,2));
        if (V0 > 1e-2) {
          p[i].x += p[i].vx*t - mu*g*pow(t,2)*p[i].vx/2;    //Loi de frottements
          p[i].y += p[i].vy*t - mu*g*pow(t,2)*p[i].vy/2;
          p[i].vx -= mu*g*t*p[i].vx/2;
          p[i].vy -= mu*g*t*p[i].vy/2;
        } else {
          p[i].vx = 0;
          p[i].vy = 0;
        }
    }
}


void colCalcV(Particle* p, int i, int j) {        //Calcul des vitesses post-collision interparticulaire
    double dx = p[i].x - p[j].x;
    double dy = p[i].y - p[j].y;
    double dvx = p[i].vx - p[j].vx;
    double dvy = p[i].vy - p[j].vy;

    double ps = (dx*dvx + dy*dvy)/(pow(dx, 2) + pow(dy, 2));

    if (j < 16) {       // Vraie boule
        p[i].vx -= ps*dx;
        p[i].vy -= ps*dy;
        p[j].vx += ps*dx;
        p[j].vy += ps*dy;
    } else {            // Particule d'angle
        p[i].vx -= 2*ps*dx;
        p[i].vy -= 2*ps*dy;
    }
}


void colParCalcV(Particle* p, Paroi* par, int i, int k) {    //Calcul des vitesses post-collision avec paroi
    double alpha = par[k].alpha;
    double vx = p[i].vx;
    double vy = p[i].vy;
    double v0 = sqrt(pow(vx,2) + pow(vy,2));

    double beta = atan(vy/vx);

    if (vx < 0) { // atan -> valeur dans ]-pi/2, pi/2[  ==> pb qd vx < 0
        beta += M_PI;
    }

    p[i].vx = v0*cos(2*alpha - beta);
    p[i].vy = v0*sin(2*alpha - beta);
}


int getch(void) {                                 // Récupération des touches du clavier, sans avoir à appuyer sur "Enter" entre chaque
    // Nous avions prévu de faire la récupération des touches avec la bibliothèque SDL, mais nous n'avons pas pu l'installer
    // => Solution de substitution clé en main : https://forums.macrumors.com/threads/getch-command-on-mac.2011972/
    char chbuf[1];
    struct termios oldstate, newstate;
    tcgetattr(0, &oldstate);
    newstate = oldstate;
    newstate.c_lflag &= ~ICANON;
    newstate.c_lflag &= ~ECHO;
    tcsetattr(0, TCSANOW, &newstate);
    read(0, &chbuf, 1);
    tcsetattr(0, TCSANOW, &oldstate);
    return *chbuf;
}


void getarrow(Particle* p, Trou* T, Paroi* par, Graphics* pgw, char joueur) {     // Récupération du tir (angle et puissance) + tracé de la flèche de tir
  printf("> Orientez votre tir avec Q/A et D/E. Gérez la puissance avec Z et S puis validez avec Enter ou Space\n");
  double alpha = M_PI/2 + (drand48()-0.5) * M_PI/6;
  double v0 = 5.;
  (*pgw).draw(par,T,p, 1000, alpha, v0, joueur);
  int ch = 0;
  while (ch != 10 && ch != 32) {      // "Enter" ou "Space" pour tirer
    ch = getch();
    if (ch == 'q') {                  // Rotation lente sens positif
      alpha += M_PI/180;
    }
    if (ch == 'a') {                  // Rotation rapide sens positif
      alpha += 10*M_PI/180;
    }
    if (ch == 'd') {                  // Rotation lente sens négatif
      alpha -= M_PI/180;
    }
    if (ch == 'e') {                  // Rotation rapide sens négatif
      alpha -= 10*M_PI/180;
    }
    if (ch == 'z' && v0 < 12.5) {     // Augmente la puissance (en limitant la vitesse max)
      v0 += 0.25;
    }
    if (ch == 's' && v0 > 0.25) {     // Diminue la puissance (en gardant une vitesse non nulle)
      v0 -= 0.25;
    }
    (*pgw).draw(par,T,p, 1000, alpha, v0, joueur);    // Tracé de la frame avec la flèche
  }
  p[0].vx = v0*cos(alpha);
  p[0].vy = v0*sin(alpha);
}


void placeBlanche(Particle* p, Graphics* pgw, double lmin, double lmax, double L, double d, double r, Paroi* par, Trou* T, char joueur) {   // Replace la boule blanche si elle n'est plus en jeu
  printf("> Placez la boule blanche avec q/a et d/e, puis validez avec Enter ou Space\n");

  p[0].x = (lmax - lmin)/2;           // Par défaut : au milieu
  p[0].y = L/4;
  (*pgw).draw(par,T,p, 1000, 0, 0, joueur);    // Affichage

  int ch = 0;
  while (ch != 10 && ch != 32) {      // "Enter" ou "Space" pour valider le placement
    ch = getch();
    if (ch == 'q' && p[0].x > lmin + r + d + 0.01) {        // Décalage lent gauche
      p[0].x -= 0.01;
    }
    if (ch == 'a' && p[0].x > lmin + r + d + 0.1) {         // Décalage rapide gauche
      p[0].x -= 0.1;
    }
    if (ch == 'd' && p[0].x < lmax - r - d - 0.01) {        // Décalage lent droite
      p[0].x += 0.01;
    }
    if (ch == 'e' && p[0].x < lmax - r - d - 0.1) {        // Décalage rapide droite
      p[0].x += 0.1;
    }
    (*pgw).draw(par,T,p, 1000, 0, 0, joueur);    // Affichage
  }
  p[0].enJeu = 1;

}


char arret(Particle* p) {                         // Test sur la vitesse des 16 boules pour lancer le prochain tir
  for (int i = 0; i<16; i++) {
    if (p[i].vx != 0 || p[i].vy != 0) {
      return(0);
    }
  }
  return(1);
}


double dist(double x1, double y1, double x2, double y2) {     // Distance entre deux points (x1, y1) et (x2, y2)
  return(sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)));
}


char tombe(Particle* p, Trou* T, int i, double lmax, double d) {    // Gère si une boule tombe dans un trou ou non
  char kray = 0,kpl = 0;          // Nombre de boules de chaque type (rayées ou pleines) tombées dans les poches (utile pour l'affichage des boules sur le côté)

  for (int k = 0; k<16; k++) {    // Calcule kray et kpl
    if (p[k].plein == 1 && p[k].enJeu == 0){
      kpl++;
    } else if (p[k].plein == 0 && p[k].enJeu == 0){
      kray++;
    }
  }

  for (int j = 0;j<6;j++) {       // Pour chaque trou
    if (dist(p[i].x,p[i].y,T[j].xtrou,T[j].ytrou) < T[j].radius) {    // Teste si la distance entre la particule et le centre du trou est inférieure au rayon du trou
      p[i].enJeu = 0;             // Si oui, on la met hors jeu et on annule sa vitesse
      p[i].vx = 0;
      p[i].vy = 0;
      if (p[i].plein) {
        if (i == 0) {             // Blanche : enlevée du terrain
          p[i].x = -lmax;
          p[i].y = -lmax;
          return(3);
        } else if (i == 8) {      // Noire : enlevée du terrain
          p[i].x = -lmax;
          p[i].y = -lmax;
          return(4);
        } else {                  // Autre pleines : affichées sur le bord gauche
          p[i].x = -d;
          p[i].y = 2*d + kpl*2*d;
          return(1);
        }
      } else {                    // Rayées : affichées sur le bord droit
        p[i].x = lmax + d;
        p[i].y = 2*d + kray*2*d;
        return(2);
      }
    }
  }
  return(0);
}


char coupsuivant(double FPS, int Np, double diameter, double L, double d, double lmin, double lmax, Graphics * pgw, Particle* p, Paroi* par,Trou* T, char joueur) {
    Event* e = (Event*) malloc( Np*(Np+24) * sizeof(Event) );

    long ne;         // Nombre d'évènements écrits dans e
    long n;

    if (joueur == 1) {
      printf("\n\nTour du joueur 1 (boules pleines)\n");
    } else {
      printf("\n\nTour du joueur 2 (boules rayées)\n");
    }
    printf("-------------------------------\n\n");


    if(p[0].enJeu == 0){      // Si la blanche est hors-jeu
      placeBlanche(p, pgw, lmin, lmax, L, d, diameter/2, par, T, joueur);
    }

    getarrow(p,T,par,pgw,joueur);     // Demande au joueur de tirer la boule

    (*pgw).draw(par,T,p, FPS, 0, 0, joueur);

    double taff = 1/FPS;            // Temps avant prochain affichage
    long naff = 0;                  // Nombre d'affichages

    long ncol;
    double tcol;

    double vx,vy,x0,y0,x1,y1,alpha,lpar,r = diameter/2,qprov,x2,y2,xc,yc,lcosa,xmx2,t;    // Variables temporaires pour les calculs
    int k0,k1;

    char estTombe = 0, nouvJoueur = 0, premTouche = 0;                      // Sert à appliquer les règles du jeu : déterminer qui joue au prochain tour

    while (arret(p) == 0) {
        // Actualisation des collisions à venir

        n = 0;
        for (int i = 0; i < Np; i++) {                                      // On parcourt toutes les boules
            vx = p[i].vx ; vy = p[i].vy ; x0 = p[i].x ; y0 = p[i].y;        // Raccourcis

            if (vx != 0 || vy != 0) {                                       // On ne travaille que sur les partucules mobiles

                for (int k = 0; k < 18; k++) {                              // k parcourt les 18 parois
                    x1 = par[k].x1; y1 = par[k].y1; alpha = par[k].alpha;   lpar = par[k].l;  // Raccourcis

                    x2 = x1 + r*sin(alpha);                                 // (x2, y2) coordonnées limites à r de la barrière (pour tests avec le centre des boules)
                    y2 = y1 - r*cos(alpha);

                    if (k == 7 || k == 10 ){    // MURS VERTICAUX alpha = + pi/2

                        yc = y0 + (x2 - x0)/(vx/vy);                        // Calcul de l'ordonnée de la collision
                        if ((lpar > 0 && yc >= y1 && yc <= y1 + lpar) || (lpar < 0 && yc <= y1 && yc >= y1 + lpar)) {     // Si l'ordonnée de la collision est entre les ordonnées limites du mur
                            t = (yc - y0)/vy;                               // On calcule tcol, et, s'il est positif, on met à jour e
                            if (t >= 0) {
                                e[n].type = paroi;
                                e[n].i = i;
                                e[n].j = k;
                                e[n].time = t;
                                assert(e[n].time >= 0);
                                n++;
                            }
                        }

                    } else if (k == 13 || k == 16) { // MURS VERTICAUX alpha = - pi/2 : idem

                        yc = y0 + (x2 - x0)/(vx/vy);
                        if ((lpar < 0 && yc >= y1 && yc <= y1 - lpar) || (lpar > 0 && yc <= y1 && yc >= y1 - lpar)) {
                            t = (yc - y0)/vy;
                            if (t >= 0) {
                                e[n].type = paroi;
                                e[n].i = i;
                                e[n].j = k;
                                e[n].time = t;
                                assert(e[n].time >= 0);
                                n++;
                            }
                        }

                    } else {      // MURS HORIZONTAUX ET OBLIQUES

                        xc = ( y0 - y2 + x2*tan(alpha) - x0*vy/vx )/( tan(alpha) - vy/vx );       // Calcul de l'abscisse de la collision
                        lcosa = lpar*cos(alpha);      // Raccourcis
                        xmx2 = xc-x2;

                        if ((lcosa > 0 && xmx2 >= 0 && xmx2 <= lcosa) || (lcosa < 0 && xmx2 <= 0 && xmx2 >= lcosa)) {     // Si l'abscisse de la collision est entre les abscisses limites de la paroi
                            t = (xc - x0)/vx;                           // On calcule tcol, et, s'il est positif, on met à jour e
                            if (t >= 0) {
                                e[n].type = paroi;
                                e[n].i = i;
                                e[n].j = k;
                                e[n].time = t;
                                assert(e[n].time >= 0);
                                n++;
                            }
                        }

                    }
                }

                for (int j = 0; j < Np; j++) {
                    if (i != j) {
                      if (colPart(p, diameter, i, j, e, n)) {         // ajoute un évènement dans e et renvoie 1 si collision entre i et j, renvoie 0 sinon)
                        n ++;
                        if(i == 0 && premTouche == 0){                // détermine la première boule touchée par la blanche, si elle n'en a touché aucune déjà
                          if (p[j].plein && j != 8) {
                            premTouche = 1;       // Boule pleine
                          } else if (j == 8) {
                            premTouche = 4;       // Noire
                          } else {
                            premTouche = 2;       // Boule rayée
                          }
                        }
                      }
                    }
                }
            }
        }

        ne = n;                             // Nombre d'évènements écrits dans e

        // Détermination de la prochaine collision
        tcol = e[0].time;
        ncol = 0;
        for (n = 0; n < ne; n++) {          // Si on trouve un temps plus cours que e[0], on le prend
            if (e[n].time < tcol) {
                ncol = n;
                tcol = e[n].time;
            }
        }


        // Détermination si affichage ou collision en premier
        if (ne == 0 || taff < tcol) {
            // Affichage en 1er
            actPosVit(p, Np, taff);                 // Actualisation des positions et vitesses
            for (int i = 0;i<16;i++) {
              if (p[i].plein) {                     // Détermination si une boule tombe dans un trou
                estTombe += tombe(p,T,i,lmax,d);
              } else {
                estTombe += tombe(p,T,i,lmax,d);
              }
            }
            (*pgw).draw(par,T, p, FPS, 0, 0, joueur);    // On affiche
            for (n = 0; n < ne; n++) {
                e[n].time -= taff;                  // Toutes collisions avancées de taff
            }
            taff = 1/FPS;                           // Prochain affichage dans 1/FPS
        } else {
            // Collision en 1er
            actPosVit(p, Np, tcol);                 // Actualisation des positions
            for (int i = 0;i<16;i++) {
              if (p[i].plein) {                     // Détermination si une boule tombe dans un trou
                estTombe += tombe(p,T,i,lmax,d);
              } else {
                estTombe += tombe(p,T,i,lmax,d);
              }
            }
            taff -= tcol;                           // Affichage avancé de tcol

            switch (e[ncol].type) {               // Mise à jour des vitesses
                case particle:
                    colCalcV(p, e[ncol].i, e[ncol].j);
                    break;
                case paroi:
                    colParCalcV(p, par, e[ncol].i, e[ncol].j);
                    break;
                default:
                    printf("Erreur ! Type de collision invalide\n");
                    break;
            }
        }

        // Affichage d'une phrase dans le terminal et mise à jour éventuelle de nouvJoueur (prochain joueur à jouer) lorsqu'une boule tombe
        switch (estTombe) {
          case 1:     // Boule pleine
            if (joueur == 1) {
              printf("Bien joué ! Une boule dans le trou\n");
              if (nouvJoueur == 0) {
                nouvJoueur = joueur;
              }
            } else {
              printf("Aïe !! Une boule adverse est tombée...\n");
            }
            break;
          case 2:     // Boule rayée
            if (joueur == 1) {
              printf("Aïe !! Une boule adverse est tombée...\n");
            } else {
              printf("Bien joué ! Une boule dans le trou\n");
              if (nouvJoueur == 0) {
                nouvJoueur = joueur;
              }
            }
            break;
          case 3:     // Blanche
            printf("Faute !!!! La blanche est tombée... Votre adversaire devra la replacer !\n");
            if (nouvJoueur != -joueur) {
              nouvJoueur = -joueur;
            }
            break;
          case 4:     // Noire
            printf("..............................\n");
            break;
          default:
            break;
        }
        estTombe = 0;

    }   // Fin du while : le coup est fini

    // Détermination du prochin joueur à jouer, en fonction de la première boule à avoir été touchée (et des règles précédentes)
    printf("Première boule touchée : ");
    switch (premTouche) {
      case 0:   // Aucune boule touchée
        printf("aucune... ==> FAUTE\n");
        if (nouvJoueur != -joueur) {
          nouvJoueur = -joueur;
        }
        break;
      case 1:     // Boule pleine
        if (joueur == 1) {
          printf("une boule pleine ==> OK\n");
        } else {
          printf("une boule adverse... ==> FAUTE\n");
          if (nouvJoueur != -joueur) {
            nouvJoueur = -joueur;
          }
        }
        break;
      case 2:     // Boule rayée
        if (joueur == 1) {
          printf("une boule adverse... ==> FAUTE\n");
          if (nouvJoueur != -joueur) {
            nouvJoueur = -joueur;
          }
        } else {
          printf("une boule rayée ==> OK\n");
        }
        break;
      case 4:     // Noire
        printf("la noire... ==> FAUTE\n");
        if (nouvJoueur != -joueur) {
          nouvJoueur = -joueur;
        }
        break;
      default:
        break;
    }

    if (nouvJoueur == 0) {          // Si aucune des règles ci-dessus ne s'est appliquée, l'adversaire joue
      nouvJoueur = -joueur;
    }

    free(e);
    return(nouvJoueur);
}


int main(){
    // GRANDEURS GLOBALES

    double FPS = 200.;
    int Np = 16 + 12;                       // Jeu + boules modélisant les angles
    double diameter = 0.058;                 // Boule de billard américain : d = 57,2 mm
    int Pix = 900;                          // Longeur du billard total, en pixels
    double L = 2.54, l = 1.27, d = 0.05;    // Longeur, largeur de la zone de jeu, d largeur de la bande (2d largeur du bord)
    double Lmin = 0, Lmax = L + 2*d;        // Définition des dimansions de l'aire de jeu
    double lmin = 0, lmax = l + 2*d;
    srand48(time(NULL));                    // Utilisé pour position initiale de la flèche aléatoire

    // INITIALISATION

    Graphics gw(Np, Pix, lmin, lmax, Lmin, Lmax, L, l, d, diameter);  // Fenêtre graphique

    Particle p[28];
    initparticles(p, L, l, d, diameter);                              // Boules de jeu

    Paroi par[18];
    initparois(par, L, l, d);                                         // Parois

    Trou T[6];
    initTrou(T,par);

    initParoiBalls(p, par, diameter/2);                               // Boules d'angle

    // DÉBUT DU JEU

    printf("\n\n\nTravaux Pratiques de Programmation en langage C/C++ :\n");
    printf("-------------------------------------------------------\n");
    printf("          Modélisation d'un billard américain\n");
    printf("-------------------------------------------------------\n");
    printf("\nThomas Gomes & Loïc Simon / 137e Promotion,  Groupe 1\n");
    printf("© ESPCI Paris, mars-avril 2019\n\n\n");
    printf("Le terrain va s'afficher. Pour jouer, redonner le focus à ce terminal.\n\n");

    char joueur = 1;      // Le joueur 1 commence

    gw.draw(par,T,p,FPS,0,0, joueur);   // Affichage
    sleep(1);                           // Pause d'1s avant affichage de la flèche

    int mabite = 2;

    while (1) {
        joueur = coupsuivant(FPS, Np, diameter, L, d, lmin, lmax, &gw, p, par, T, joueur);    // Réalisation du coup
        if(p[8].enJeu == 0) {           // Si la noire est sortie...
          printf("\n\nVOUS ÊTES UNE SOUPIÈRE\n");
          printf("Cordialement\n");
          printf("La direction\n\n");
          break;
        }
    }
    return (0);
}

// ********** FIN DU PROGRAMME **********
