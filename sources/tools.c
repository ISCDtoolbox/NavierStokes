#include "nstokes.h"


#define TMAX     128

/* find boundary conds in list for given ref and type */
pCl getCl(pSol sol,int ref,int elt,char typ) {
  pCl     pcl;
  int     i;

  for (i=0; i<sol->nbcl; i++) {
    pcl = &sol->cl[i];
    if ( (pcl->ref == ref) && (pcl->elt == elt) && (pcl->typ == typ) )  return(pcl);
  }
  return(0);
}

/* retrieve physical properties in list */
int getMat(pSol sol,int ref,double *nu,double *rho) {
  pMat   pm;
  int    i;

  *nu  = NS_NU;
  *rho = NS_RHO;
  if ( sol->nmat == 0 )  return(1);

  for (i=0; i<sol->nmat; i++) {
    pm = &sol->mat[i];
    if ( pm->ref == ref ) {
      *nu   = pm->nu;
      *rho  = pm->rho;
      return(1);
    }
  }

  return(0);
}


static inline int nortri(double *a,double *b,double *c,double *len,double *n) {
  double    det,abx,aby,abz,acx,acy,acz;

  /* 3D triangle area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  /* average length of 2 incident edges */
  *len  = 0.5 * sqrt(abx*abx + aby*aby + abz*abz);
  *len += 0.5 * sqrt(acx*acx + acy*acy + acz*acz);

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;

  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( det < NS_EPSD )  return(0);

  det = 1.0 / sqrt(det);
  n[0] *= det;
  n[1] *= det;
  n[2] *= det;
  return(1);
}


/* curvature approximation at point ip (use triangles around point ip) */
int kappa_3d(pMesh mesh,int ip,double *n,double *len,double *kappa) {
  pPoint   p0,p1,p2,p3;
  pTria    pt;
  pTetra   ptt;
  double   tlist[TMAX],ni[3],b1[3],b2[3],ma[6],mb[3];
  double   dd,l,dz,ux,uy,uz,x,y,z,k1,k2,x2,y2,z2,xy;
  int      k,adj,ref;
  char     i,i1,i2,j2,il,ier,cls;

  *kappa = 0.0;
  memset(n,0,3*sizeof(double));

  /* find entity on interface */
  p0 = &mesh->point[ip];
  k  = p0->s3;
  if ( k == 0 )  return(0);
  pt  = &mesh->tria[k];
  ref = pt->ref;

  for (i=0; i<3; i++)
    if ( pt->v[i] == ip )  break;
  i1  = (i+1) % 3;
  j2  = (i+2) % 3;

  /* loop over triangles around ip, store in tlist */
  il = 0;
  tlist[il] = 3*k + i;

  cls = 1;  
  adj = pt->adj[i1] / 3;
  i2  = pt->adj[i1] % 3;
  pt  = &mesh->tria[adj];
  while ( (adj > 0) && (adj != k) && (pt->ref == ref) ) {
    il++;
    tlist[il] = 3*adj + (i2+1) % 3;

    i1  = (i2+2) % 3;
    adj = pt->adj[i1] / 3;
    i2  = pt->adj[i1] % 3;
    pt  = &mesh->tria[adj];
  }

  /* reverse loop if open ball */
  if ( adj != k ) {
    cls = 0;
    pt  = &mesh->tria[k];
    i2  = j2;
    adj = pt->adj[i2] / 3;
    i1  = pt->adj[i2] % 3;
    pt  = &mesh->tria[adj];
    while ( (adj > 0) && (adj != k) && (pt->ref == ref) ) {
      il++;
      tlist[il] = 3*adj + (i1+2) % 3;

      i2  = (i1+1) % 3;
      adj = pt->adj[i2] / 3;
      i1  = pt->adj[i2] % 3;
      pt  = &mesh->tria[adj];
    }
  }

  if ( il < 3 )  return(0);

  /* compute average normal */
  *len = 0.0;
  for (k=0; k<=il; k++) {
    pt  = &mesh->tria[k / 3];
    p0  = &mesh->point[pt->v[0]];
    p1  = &mesh->point[pt->v[1]];
    p2  = &mesh->point[pt->v[2]];
    ier = nortri(p0->c,p1->c,p2->c,&l,ni);
    if ( !ier )  return(0);
    n[0] += ni[0];  n[1] += ni[1];  n[2] += ni[2];
    *len += l;
  }
  *len /= il+1;

  /* unit normal */
  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd < NS_EPSD )  return(0);

  dd = 1.0 / sqrt(dd);
  n[0] *= dd;  n[1] *= dd;  n[2] *= dd;

  /* local frame */
  if ( fabs(n[0]) > fabs(n[1]) ) {
    if ( fabs(n[0]) > fabs(n[2]) ) {
      b1[0] = -(n[1]+n[2]) / n[0];
      b1[1] = b1[2] = 1.0;
    }
    else {
      b1[2] = -(n[0]+n[1]) / n[2];
      b1[0] = b1[1] = 1.0;
    }
  }
  else {
    if ( fabs(n[1]) > fabs(n[2]) ) {
      b1[1] = -(n[0]+n[2]) / n[1];
      b1[0] = b1[2] = 1.0;
    }
    else {
      b1[2] = -(n[0]+n[1]) / n[2];
      b1[0] = b1[1] = 1.0;
    }
  }
  dd = b1[0]*b1[0] + b1[1]*b1[1] + b1[2]*b1[2];
  if ( dd < NS_EPSD ) return(0);
  dd = 1.0 / sqrt(dd);
  b1[0] *= dd;  b1[1] *= dd;  b1[2] *= dd;

  b2[0] = n[1]*b1[2] - n[2]*b1[1];
  b2[1] = n[2]*b1[0] - n[0]*b1[2];
  b2[2] = n[0]*b1[1] - n[1]*b1[0];

  /* curvature evaluation */
  p0 = &mesh->point[ip];
  dz = 0.0;
  for (k=0; k<=il; k++) {
    pt = &mesh->tria[k / 3];
    i  = k % 3;
    i1 = (i+1) % 3;
    i2 = (i+2) % 3;
    p0 = &mesh->point[pt->v[i]];
    p1 = &mesh->point[pt->v[i1]];
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    /* coordinates of point ip in local frame */
    x = ux*b1[0] + uy*b1[1] + uz*b1[2];
    y = ux*b2[0] + uy*b2[1] + uz*b2[2];
    z = ux*n[0]  + uy*n[1]  + uz*n[2];

    x2 = x*x;  y2 = y*y;  z2 = z*z;  xy = 2*x*y;
    if ( dz < z2 )  dz = z2;

    ma[0] += x2*x2;
    ma[1] += x2*xy;
    ma[2] += x2*y2;
    ma[3] += 4*x2*y2;
    ma[4] += xy*y2;
    ma[5] += y2*y2;

    mb[0] += x2*z;
    mb[1] += xy*z;
    mb[2] += y2*z;
  }
  // CLOSED ?

  /* flat surface */
  if ( dz < NS_EPSD )  return(1);

  return(1);
}


/* curvature approximation at point ip (via edges around ip) */
int kappa_2d(pMesh mesh,int ip,double *n,double *len,double *kappa) {
  pPoint   p0,p1,p2;
  pTria    pt;
  double   dd,ux,uy,vx,vy,k1,k2,d1,d2,nx,ny;
  int      k,adj,ip1,ip2,ref;
  char     i,i1,i2,j2;

  *kappa = 0.0;

  /* find entity on interface */
  p0 = &mesh->point[ip];
  k  = p0->s2;
  if ( k == 0 )  return(0);
  pt  = &mesh->tria[k];
  ref = pt->ref;

  for (i=0; i<3; i++)
    if ( pt->v[i] == ip )  break;
  i1  = (i+1) % 3;
  j2  = (i+2) % 3;
  ip1 = pt->v[i1];
  ip2 = pt->v[j2];

  /* loop over triangles around ip */
  adj = pt->adj[i1] / 3;
  i2  = pt->adj[i1] % 3;
  pt  = &mesh->tria[adj];
  while ( (adj > 0) && (adj != k) && (pt->ref == ref) ) {
    ip2 = pt->v[i2];
    i1  = (i2+2) % 3;
    adj = pt->adj[i1] / 3;
    i2  = pt->adj[i1] % 3;
    pt  = &mesh->tria[adj];
  }
  /* reverse loop if open ball */
  if ( adj != k ) {
    pt  = &mesh->tria[k];
    i2  = j2;
    adj = pt->adj[i2] / 3;
    i1  = pt->adj[i2] % 3;
    pt  = &mesh->tria[adj];
    while ( (adj > 0) && (adj != k) && (pt->ref == ref) ) {
      ip1 = pt->v[i1];
      i2  = (i1+1) % 3;
      adj = pt->adj[i2] / 3;
      i1  = pt->adj[i2] % 3;
      pt  = &mesh->tria[adj];
    }
  }

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  vx = p2->c[0] - p0->c[0];
  vy = p2->c[1] - p0->c[1];
  /* average length of 2 edges */
  *len  = sqrt(ux*ux + uy*uy);
  *len += sqrt(vx*vx + vy*vy);
  *len  = *len / 2.0;

  /* outer normal */
  n[0] = uy - vy;
  n[1] = vx - ux;
  dd   = sqrt(n[0]*n[0] + n[1]*n[1]);
  if ( dd > NS_EPSD ) {
    n[0] /= dd;
    n[1] /= dd;
  }
  /* free interface (not closed) */
  if ( p2->ref != p0->ref || p1->ref != p0->ref ) {
    return(1);
  }

  /* compute algebraic curvature */
  d1 = ux*n[0] + uy*n[1];
  d2 = vx*n[0] + vy*n[1];
  if ( fabs(d1) > NS_EPSD )
    k1 = (ux*ux + uy*uy) / d1;
  else
    return(1);
  if ( fabs(d2) > NS_EPSD )
    k2 = (vx*vx + vy*vy) / d2;
  else
    return(1);

  if ( fabs(k1+k2) > NS_EPSD )
    *kappa = 4.0 / (k1+k2);

  return(1);
}
