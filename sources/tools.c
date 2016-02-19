#include "nstokes.h"


/* curvature approximation at point ip */
double kappa_2d(pMesh mesh,int ip,double *n,double *len) {
  pPoint   p0,p1,p2;
  pTria    pt;
  double   dd,ux,uy,vx,vy,k1,k2,d1,d2,nx,ny,kappa;
  int      k,adj,ip1,ip2,ref;
  char     i,i1,i2,j2;

  kappa = 0.0;
  p0 = &mesh->point[ip];
  k  = p0->s;
  if ( k == 0 )  return(kappa);
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
puts("ICI");
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
  *len  = sqrt(ux*ux + uy*uy);
  *len += sqrt(vx*vx + vy*vy);

  /* outer normal */
  n[0] = uy - vy;
  n[1] = vx - ux;
  dd   = sqrt(n[0]*n[0] + n[1]*n[1]);
  n[0] /= dd;
  n[1] /= dd;

  /* free interface (not closed) */
  if ( p2->ref != p0->ref || p1->ref != p0->ref ) {
    return(kappa);
  }

  /* compute algebraic curvature */
  d1 = ux*n[0] + uy*n[1];
  d2 = vx*n[0] + vy*n[1];
  k1 = (ux*ux + uy*uy) / d1;
  k2 = (vx*vx + vy*vy) / d2;
  kappa = 4.0 / (k1+k2);

  return(kappa);
}
