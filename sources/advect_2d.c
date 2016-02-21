#include "nstokes.h"


/* barycentric coordinates of point c[] in iel=p0,p1,p2 */
static inline int bar_2d(pPoint pp[3],double *c,int iel,double *cb) {
  double    det;
  char      i,i1,i2;

  det = (pp[1]->c[0] - pp[0]->c[0]) * (pp[2]->c[1] - pp[0]->c[1])
      - (pp[1]->c[1] - pp[0]->c[1]) * (pp[2]->c[0] - pp[0]->c[0]);
  if ( fabs(det) < NS_EPSD )  return(0);
  det = 1.0 / det;

  /* barycentric coordinate, by use of Cramer's formulas */
  for (i=0; i<2; i++) {
    i1 = (i+1) % 3;
    i2 = (i+2) % 3;
    cb[i]  = (pp[i1]->c[0]-c[0])*(pp[i2]->c[1]-c[1]) - (pp[i1]->c[1]-c[1])*(pp[i2]->c[0]-c[0]);
    cb[i] *= det;
  }
  cb[2] = 1.0 - cb[0] - cb[1];

  return(1);
}


/* return velocity P1_interpolation in v in element iv[3] */
static inline double vecint_P1_2d(double *un,int *iv,double *cb,double *v) {
  double   *u0,*u1,*u2,dd;

  u0 = &un[2*(iv[0]-1)];
  u1 = &un[2*(iv[1]-1)];
  u2 = &un[2*(iv[2]-1)];

  /* P1 interpolate of the speed */   
  v[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];  
  v[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  dd   = sqrt(v[0]*v[0] + v[1]*v[1]);

  return(dd);
}


/* find element containing c, starting from nsdep, return baryc. coord */
static int locelt_2d(pMesh mesh,int nsd,double *c,double *cb) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    ax,ay,bx,by,cx,cy,eps;
  double    aire1,aire2,aire3,dd; 
  int       i,sgn,nsf,nsp;

  nsf  = nsd;
  nsp  = 0;
  mesh->mark = ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tria[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;

    /* area of triangle */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    sgn = dd > 0.0 ? 1 : -1;
    eps = sgn == 1 ? -NS_EPS*dd : NS_EPS*dd;

    /* barycentric */
    bx = p1->c[0] - c[0];
    by = p1->c[1] - c[1];
    cx = p2->c[0] - c[0];
    cy = p2->c[1] - c[1];

    /* p in half-plane lambda_0 > 0 */
    aire1 = sgn*(bx*cy - by*cx);
    if ( aire1 < eps ) {
      nsp = nsf;
      nsf = pt->adj[0] / 3;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    ax = p0->c[0] - c[0];
    ay = p0->c[1] - c[1];
    aire2 = sgn*(cx*ay - cy*ax);
    if ( aire2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 3;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    aire3 = sgn*dd - aire1 - aire2;
    if ( aire3 < eps ) {
      nsp = nsf;
      nsf = pt->adj[2] / 3;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    aire1 = NS_MAX(aire1,0.0);
    aire2 = NS_MAX(aire2,0.0);
    aire3 = NS_MAX(aire3,0.0);
    dd    = aire1 + aire2 + aire3;
    if ( dd > NS_EPSD ) {
      dd = 1.0 / dd;
      cb[0] = aire1 * dd;
      cb[1] = aire2 * dd;
      cb[2] = aire3 * dd;
    }
    return(nsf);
  }

  /* no need for exhaustive search */
  return(nsp);
}


/* computes the characteristic line emerging from point with barycentric coordinates cb
 in triangle iel and follows the characteristic on a length dt at most. Most often has to 
 cross the boundary of the current triangle, thus stores in it the new triangle, in cb
 the barycentric coordinates in the new triangle of the crossing point, and updates dt 
 with the remaining time to follow characteristic line */
static int travel_2d(NSst *nsst,double *cb,int *iel,double *dt) {
  pTria       pt;
  pPoint      p[3];
  double     *u0,*u1,*u2,m[3],ddt,tol,ux,uy,c[2],cb1[3];
  int         k;
  char        i,i0,i1,i2;

  tol = *dt;
  k   = *iel;
  pt  = &nsst->mesh.tria[k];

  p[0] = &nsst->mesh.point[pt->v[0]];
  p[1] = &nsst->mesh.point[pt->v[1]];
  p[2] = &nsst->mesh.point[pt->v[2]];

  /* velocity at each vertex of iel */
  u0 = &nsst->sol.un[2*(pt->v[0]-1)];
  u1 = &nsst->sol.un[2*(pt->v[1]-1)];
  u2 = &nsst->sol.un[2*(pt->v[2]-1)];

  /* u = P1 velocity at the point of barycentric coordinates cb */
  ux = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];
  uy = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  if ( ux*ux+uy*uy < NS_EPSD )  return(0);

  /* endpoint of the characteristic line starting from cb */
  /* segment is of norm sqrt(ux*ux+uy*uy), and is to be followed for time dt, 
     in the - direction */
  c[0] = (cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0]) - tol*ux;   
  c[1] = (cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1]) - tol*uy;

  /* barycentric coordinate of uu in the current triangle */
  if ( !bar_2d(p,c,k,cb1) )  return(0);

  /* check if endpoint is in k */
  for (i=0; i<3; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 3 ) {
    memcpy(cb,cb1,3*sizeof(double));
    *dt -= tol;
    return(*dt > 0.0);
  }

  /* argument of the smallest value among 3 */
  ddt = NS_TGV;
  i0  = -1;
  for (i=0; i<3; i++) {
    m[i] = cb[i] - cb1[i]; 
    if ( m[i] > 0.0 ) {
      if ( tol*cb[i]/m[i] < ddt ) {
        ddt = tol*cb[i]/m[i];
        i0  = i;
      } 
    }
  }
  /* case when advection stays in same triangle */
  if ( ddt > tol ) {
    memcpy(cb,cb1,3*sizeof(double));
    *dt -= tol;
  }
  /* advection goes out of triangle: advect a minimum value */
  if ( ddt < NS_EPS*tol ) {
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - NS_EPS*ux;
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - NS_EPS*uy;
    /* find the new triangle */
    k = locelt_2d(&nsst->mesh,k,c,cb); 
    if ( !k )  return(0);
    *iel = k;
    *dt -= NS_EPS;
  }
  else {
    /* barycentric coordinates of the exit point */   
    for (i=0; i<3; i++)
      cb1[i] = cb[i] - ddt* m[i]/tol;
    *dt -= ddt;

    /* find output triangle */
    if ( !pt->adj[i0] ) {
      memcpy(cb,cb1,3*sizeof(double));
      return(0);
    }
    else {
      i1   = (i0+1) % 3;
      i2   = (i0+2) % 3;
      *iel = pt->adj[i0] / 3;
      i0   = pt->adj[i0] % 3; 
      cb[i0]         = 0.0;
      cb[(i0+1) % 3] = cb1[i2];
      cb[(i0+2) % 3] = cb1[i1];
    }
  }

  return(*dt > 0.0);
}


static int nxtpt_2d(NSst *nsst,int *iel,double *c,double *cb,double step,double *v) {
  double  cc[2];
  int     k;

  k = *iel;
  cc[0] = c[0] - step*v[0];
  cc[1] = c[1] - step*v[1];
  k = locelt_2d(&nsst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */  
  c[0] = cc[0];
  c[1] = cc[1];
  vecint_P1_2d(nsst->sol.un,nsst->mesh.tria[k].v,cb,v);
  *iel = k;

  return(1);
}


/* advect P2 solution u0, result in u */
int advect_P2_2d(NSst *nsst) {
  return(1);
}


/* advect P1 solution u0, result in u */
int advect_P1_2d(NSst *nsst) {
  pTria    pt;
  pPoint   ppt;
  double  *u0,*u1,*u2,cb[3],v[2],c[2],dt,dte,norm,st;  
  int      i,j,k,iel,kp,ns;

  dt = nsst->sol.dt;
  ns = 10;
  st = dt / ns;

  ++nsst->mesh.mark;
  for (k=1; k<=nsst->info.np; k++) {
    ppt = &nsst->mesh.point[k];

    /* velocity at point p */
    v[0] = nsst->sol.un[2*(k-1)+0];
    v[1] = nsst->sol.un[2*(k-1)+1];
    norm = sqrt(v[0]*v[0] + v[1]*v[1]);
    if ( norm < NS_EPSD )  continue;
    
    pt = &nsst->mesh.tria[ppt->s];
    for (i=0; i<3; i++)  if ( pt->v[i] == k )  break;

    /* barycentric coordinates of point p in triangle iel */
    memset(cb,0,3*sizeof(double));
    cb[i]         = 1.0;

    /* next point = foot of the characteristic line */
    c[0] = ppt->c[0];
    c[1] = ppt->c[1];
    dte  = dt;

    for (iel=ppt->s,j=0; j<ns; j++) {
      kp = iel;
      if ( nxtpt_2d(nsst,&iel,c,cb,st,v) < 1 )  break;
      dte -= st;
    }
    if ( j < ns && !iel ) {
      iel = kp;
      while ( travel_2d(nsst,cb,&iel,&dte) );
    }

    /* interpolate value at foot */
    pt = &nsst->mesh.tria[iel];
    
    u0 = &nsst->sol.un[2*(pt->v[0]-1)];
    u1 = &nsst->sol.un[2*(pt->v[1]-1)];
    u2 = &nsst->sol.un[2*(pt->v[2]-1)];
    nsst->sol.u[2*(k-1)+0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];
    nsst->sol.u[2*(k-1)+1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  }

  return(1);
}

