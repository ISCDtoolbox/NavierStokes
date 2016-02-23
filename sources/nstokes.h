#ifndef _NSTOKES_H
#define _NSTOKES_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "ns_calls.h"

#define NS_VER   "5.2a"
#define NS_REL   "Jan. 29, 2016"
#define NS_CPY   "(C) Copyright 2006- , ICS-SU"

#define NS_NU     1.0
#define NS_RHO    1.0
#define NS_MAT    50
#define NS_CL     50
#define NS_RES    1.e-6
#define NS_EPS    1.e-6
#define NS_EPS2   1.e-12
#define NS_MAXIT  10000
#define NS_TGV    1.e+30
#define NS_EPSD   1.e-200

#define NS_MAX(a,b) (((a) < (b)) ? (b) : (a))
#define NS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define NS_MIN3(a,b,c) ( (a) < (b) ? ((a)<(c) ? (a) : (c)) : ((b)<(c) ? (b) : (c)) )
#define NS_MAX3(a,b,c) ( (a) > (b) ? ((a)>(c) ? (a) : (c)) : ((b)>(c) ? (b) : (c)) )


/* data structures */
typedef struct {
  double  c[3];         /* coordinates */
  int     ref,s,new;    /* label and seed for simplex, new: packed value */
  char    flag;         /* used in advect to mark vertex */
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[3],ref;
} Edge;
typedef Edge * pEdge;
  
typedef struct {
  int     v[6],adj[3],ref,mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[10],adj[4],ref,mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
	int      dim,ver;
	int      np,np2,na,nt,ne,npi,nai,nti,nei;
  char     verb,typ,zip,mfree;
  mytime   ctim[TIMEMAX];
} Info;

typedef struct {
  int       mark;
  char     *name;
  Point    *point;
  Edge     *edge;
  Tria     *tria;
  Tetra    *tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  double  u[3];
  int     ref;
  char    typ,elt,att;
} Cl;
typedef Cl * pCl;

typedef struct {
  double  nu,rho;
  int     ref;
} Mat;
typedef Mat * pMat;

typedef struct {
	int      nit,nt,ts,nbcl,nmat;
  double  *u,*p,*F,*un,dt,mt,res,gr[3],tim;
  char    *namein,*nameout,*namepar,cltyp,clelt,sim;
  Cl      *cl;
  Mat     *mat;
} Sol;
typedef Sol * pSol;

struct _NSst{
  Mesh    mesh;
	Sol     sol;
	Info    info;
};

/* prototypes */
int  loadMesh(NSst *nsst);
int  loadSol(NSst *nsst);
int  saveSol(NSst *nsst,int it);
int  pack_2d(NSst *nsst);
int  pack_3d(NSst *nsst);
int  unpack(NSst *nsst);
int  hashel_2d(NSst *nsst);
int  hashel_3d(NSst *nsst);
int  addnod_2d(NSst *nsst);
int  addnod_3d(NSst *nsst);
pCl  getCl(pSol sol,int ref,int elt,char typ);
int  getMat(pSol sol,int ref,double *nu,double *rho);
int  advect_P1_2d(NSst *nsst);
int  advect_P2_2d(NSst *nsst);
int  advect_P1_3d(NSst *nsst);
int  advect_P2_3d(NSst *nsst);
int  nstokes1_2d(NSst *nsst);
int  nstokes1_3d(NSst *nsst);

double kappa_2d(pMesh mesh,int ip,double *n,double *len);


#endif


