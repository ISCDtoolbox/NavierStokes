#include "nstokes.h"


static double det_3d(double *cp0,double *cp1,double *cp2,double *cp3) {
  double  det,x01,x02,x03,y01,y02,y03,z01,z02,z03; 

  x01 = cp1[0] - cp0[0];  y01 = cp1[1] - cp0[1];  z01 = cp1[2] - cp0[2];
  x02 = cp2[0] - cp0[0];  y02 = cp2[1] - cp0[1];  z02 = cp2[2] - cp0[2];
  x03 = cp3[0] - cp0[0];  y03 = cp3[1] - cp0[1];  z03 = cp3[2] - cp0[2];

  det = x01*(y02*z03-y03*z02) + y01*(z02*x03-z03*x02) + z01*(x02*y03-x03*y02);
  return(det);  
}


/* barycentric coordinates of point c[] in iel=p0,p1,p2,p3 */
static inline int bar_3d(pPoint pp[4],double *c,int iel,double *cb) { 
  double    det,cp0[3],cp1[3],cp2[3],cp3[3];
  char      i;

  for (i=0; i<3; i++) {
    cp0[i] = pp[0]->c[i];
    cp1[i] = pp[1]->c[i];
    cp2[i] = pp[2]->c[i];
    cp3[i] = pp[3]->c[i];
  }
  det = det_3d(cp0,cp1,cp2,cp3);
  if ( fabs(det) < NS_EPSD )  return(0);
  det = 1.0 / det;

  /* barycentric coordinate */
  cb[0] =  det_3d(c,cp1,cp2,cp3) * det;
  cb[1] = -det_3d(c,cp2,cp3,cp0) * det;
  cb[2] =  det_3d(c,cp3,cp0,cp1) * det;
  cb[3] = 1.0 - cb[0] - cb[1] - cb[2];

  return(1);
}


/* return velocity interpolation in v in element iv[3] */
static double vecint_P1_3d(double *un,int *iv,double *cb,double *v) {
  double   *u0,*u1,*u2,*u3,dd;

  u0 = &un[3*(iv[0]-1)+1];
  u1 = &un[3*(iv[1]-1)+1];
  u2 = &un[3*(iv[2]-1)+1];
  u3 = &un[3*(iv[3]-1)+1];

  /* P1 interpolate of the speed */   
  v[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0] + cb[3]*u3[0];
  v[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1] + cb[3]*u3[1];
  v[2] = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2] + cb[3]*u3[2];
  dd   = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  return(dd);
}


/* find element containing c, starting from nsdep, return baryc. coord */
static int locelt_3d(pMesh mesh,int nsd,double *c,double *cb) {
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   eps,vto,vol1,vol2,vol3,vol4,dd; 
  int      i,nsf,nsp;

  nsf = nsd;
  nsp = nsd;
  ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tetra[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;
    
    /* measure of element */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    /* barycentric and volume */
    bx  = p1->c[0] - p0->c[0];
    by  = p1->c[1] - p0->c[1];
    bz  = p1->c[2] - p0->c[2];
    cx  = p2->c[0] - p0->c[0];
    cy  = p2->c[1] - p0->c[1];
    cz  = p2->c[2] - p0->c[2];
    dx  = p3->c[0] - p0->c[0];
    dy  = p3->c[1] - p0->c[1];
    dz  = p3->c[2] - p0->c[2];

    /* test volume */
    vx  = cy*dz - cz*dy;
    vy  = cz*dx - cx*dz;
    vz  = cx*dy - cy*dx;
    vto = bx*vx + by*vy + bz*vz;
		eps = NS_EPS*vto;

    /* barycentric */
    apx = c[0] - p0->c[0];
    apy = c[1] - p0->c[1];
    apz = c[2] - p0->c[2];

    /* p in 2 */
    vol2  = apx*vx + apy*vy + apz*vz;
    if ( vol2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 4;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 3 */
    vx  = by*apz - bz*apy;
    vy  = bz*apx - bx*apz;
    vz  = bx*apy - by*apx;
    vol3 = dx*vx + dy*vy + dz*vz;
    if ( vol3 < eps ) {
			nsp = nsf;
      nsf = pt->adj[2] / 4;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 4 */
    vol4 = -cx*vx - cy*vy - cz*vz;
    if ( vol4 < eps ) {
			nsp = nsf;
      nsf = pt->adj[3] / 4;
      if ( !nsf ) {
        cb[3] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 1 */
    vol1 = vto - vol2 - vol3 - vol4;
    if ( vol1 < eps ) {
			nsp = nsf;
      nsf = pt->adj[0] / 4;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    dd = fabs(vol1+vol2+vol3+vol4);
    if ( dd > NS_EPSD ) {
	    dd = 1.0 / dd;
      cb[0] = fabs(vol1) * dd;
      cb[1] = fabs(vol2) * dd;
      cb[2] = fabs(vol3) * dd;
      cb[3] = fabs(vol4) * dd;
    }
    return(nsf);
  }

  return(nsp);
}


/* computes the characteristic line emerging from point with barycentric coordinates cb
 in tetra iel and follows the characteristic on a length dt at most. Most often has to 
 cross the boundary of the current tetra, thus stores in it the new tetra, in cb
 the barycentric coordinates in the new tetra of the crossing point, and updates dt 
 with the remaining time to follow characteristic line */
static int travel_3d(NSst *nsst,double *cb,int *iel,double *dt) {
  pTetra      pt;
  pPoint      p[4];
  double     *u0,*u1,*u2,*u3,m[4],dd,ddt,tol,ux,uy,uz,c[3],cb1[4];
  int         k;
  char        i,i0,i1,i2,i3;

  tol  = *dt;
  k    = *iel;
  pt   = &nsst->mesh.tetra[k];
  
  p[0] = &nsst->mesh.point[pt->v[0]];
  p[1] = &nsst->mesh.point[pt->v[1]];
  p[2] = &nsst->mesh.point[pt->v[2]];
  p[3] = &nsst->mesh.point[pt->v[3]];

  /* the speed at each vertex of iel */
  u0 = &nsst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &nsst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &nsst->sol.u[3*(pt->v[2]-1)+1];
  u3 = &nsst->sol.u[3*(pt->v[3]-1)+1];

  /* u = P1 velocity at the point of barycentric coordinates cb */
  ux = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0] + cb[3]*u3[0];
  uy = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1] + cb[3]*u3[1];
  uz = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2] + cb[3]*u3[2];  
  if ( ux*ux+uy*uy+uz*uz < NS_EPSD )  return(0);
  
  /* endpoint of the characteristic line starting from cb */
  /* segment is of norm sqrt(ux*ux+uy*uy+uz*uz), and is to be followed for time dt, in the - direction */ 
  c[0] = (cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] + cb[3]*p[3]->c[0]) - ux;
  c[1] = (cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] + cb[3]*p[3]->c[1]) - uy;
  c[2] = (cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] + cb[3]*p[3]->c[2]) - uz;

  /* barycentric coordinate of uu in the current tetra */
  if ( !bar_3d(p,c,k,cb1) )  return(0);

  /* check if endpoint is in k */
  for (i=0; i<4; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 4 ) {
    memcpy(cb,cb1,4*sizeof(double));
    *dt -= tol;
    return(*dt > 0.0);
  }

  /* ddt = least value of cb[i]/m[i] = (weighted) distance to lambda_i0 = 0 */
  /* ddt = interval of time elapsed before going out of the current element */
  /* condition m[i] > 0 ensures going out of element in "marching direction"*/  
  ddt = NS_TGV;
  i0  = -1;
  for (i=0; i<4; i++) {
    m[i] = cb[i] - cb1[i];
    if ( m[i] > 0.0 ) {
      if ( tol*cb[i]/m[i] < ddt ) {
        ddt = tol*cb[i]/m[i];
        i0  = i;
      }
    }
  }

  /* incomplete advection remain in element */
  if ( ddt > tol ) {
    for (i=0; i<4; i++)  cb[i] = cb[i] - tol*m[i];
    *dt -= tol;
  }
  /* advection goes out of element: advect a minimum value */
  if ( ddt < NS_EPS*tol ) {
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - NS_EPS*ux;
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - NS_EPS*uy;
    c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - NS_EPS*uz;
    /* find the new element */
    k = locelt_3d(&nsst->mesh,k,c,cb); 
    if ( k < 1 )  return(0);
    *iel = k;
    *dt -= NS_EPS;
  }
  else {
    i1 = (i0+1) % 4;
    i2 = (i1+1) % 4;
    i3 = (i2+1) % 4;

    /* barycentric coordinates of the exit point */ 
    cb1[i0] = 0.0;
    cb1[i1] = cb[i1] - ddt * m[i1];
    cb1[i2] = cb[i2] - ddt * m[i2];
    cb1[i3] = 1.0 - cb1[i1] - cb1[i2];
    if ( cb1[i1] < NS_EPS2 ) {
	    cb1[i2] -= (NS_EPS2-cb1[i1]) / 2;
      cb1[i3] -= (NS_EPS2-cb1[i1]) / 2;
      cb1[i1]  =  NS_EPS2;
    }
    else if ( cb1[i1] > 1.0-NS_EPS2 ) {
	    cb1[i2] += (cb1[i1]-1.0+NS_EPS2) / 2;
      cb1[i3] += (cb1[i1]-1.0+NS_EPS2) / 2;
      cb1[i1]  = 1.0-NS_EPS2;
    }
    if ( cb1[i2] < NS_EPS2 ) {
	    cb1[i3] -= (NS_EPS2-cb1[i2]) / 2;
      cb1[i1] -= (NS_EPS2-cb1[i2]) / 2;
      cb1[i2]  = NS_EPS2;
    }
    else if ( cb1[i2] > 1.0-NS_EPS2 ) {
	    cb1[i3] += (cb1[i2]-1.0+NS_EPS2) / 2;
      cb1[i1] += (cb1[i2]-1.0+NS_EPS2) / 2;
      cb1[i2]  = 1.0-NS_EPS2;
    }
    *dt -= ddt;

    /* coordinates of the point of barycentric coordinates cb1 */
    c[0] = cb1[0]*p[0]->c[0] + cb1[1]*p[1]->c[0] + cb1[2]*p[2]->c[0] + cb1[3]*p[3]->c[0];
    c[1] = cb1[0]*p[0]->c[1] + cb1[1]*p[1]->c[1] + cb1[2]*p[2]->c[1] + cb1[3]*p[3]->c[1];
    c[2] = cb1[0]*p[0]->c[2] + cb1[1]*p[1]->c[2] + cb1[2]*p[2]->c[2] + cb1[3]*p[3]->c[2];

    if ( !pt->adj[i0] ) {
      memcpy(cb,cb1,4*sizeof(double));
      return(0);
    }
    else {
      *iel = pt->adj[i0] / 4;
      if ( !bar_3d(p,c,*iel,cb) )  return(0);
    }
  }

  return(*dt > 0.0);
}


static int nxtpt_3d(NSst *nsst,int *iel,double *c,double *cb,double step,double *v) {
  double  cc[3];
  int     k;

  k = *iel;
  cc[0] = c[0] - step*v[0];
  cc[1] = c[1] - step*v[1];
  cc[2] = c[2] - step*v[2];
  k = locelt_3d(&nsst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */  
  c[0] = cc[0];
  c[1] = cc[1];
  c[2] = cc[2];
  vecint_P1_3d(nsst->sol.un,nsst->mesh.tria[k].v,cb,v);
  *iel = k;

  return(1);
}


/* advect P2 solution u0, result in u */
int advect_P2_3d(NSst *nsst) {
  return(1);
}


/* solve advection, solution in rv */
int advect_P1_3d(NSst *nsst) {
  pTetra   pt,pt1;
  pPoint   ppt;
  double  *u0,*u1,*u2,*u3,cb[4],v[3],c[3],dt,dte,norm,st;  
  int      i,j,k,ip,iel,kp,ns;

  dt = nsst->sol.dt;
  ns = 10;
  st = dt / ns;

  ++nsst->mesh.mark;
  for (k=1; k<=nsst->info.np; k++) {
    ppt = &nsst->mesh.point[k];

    /* velocity at point p */
    v[0] = nsst->sol.un[3*(k-1)+0];
    v[1] = nsst->sol.un[3*(k-1)+1];
    v[2] = nsst->sol.un[3*(k-1)+2];
    norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if ( norm < NS_EPSD )  continue;
    
    pt = &nsst->mesh.tetra[ppt->s];
    for (i=0; i<4; i++)  if ( pt->v[i] == k )  break;

    /* barycentric coordinates of point p in triangle k */
    memset(cb,0,4*sizeof(double));
		cb[i] = 1.0;

    /* next point = foot of the characteristic line */
    c[0] = ppt->c[0];
    c[1] = ppt->c[1];
    c[2] = ppt->c[2];
    dte  = dt;
    for (iel=ppt->s,j=0; j<ns; j++) {
      kp = iel;
      if ( nxtpt_3d(nsst,&iel,c,cb,st,v) < 1 )  break;
      dte -= st;
    }
    if ( j < ns ) {
      iel = kp;
      while ( travel_3d(nsst,cb,&iel,&dte) );
    }

    /* interpolate value at foot  */
    pt1 = &nsst->mesh.tetra[iel];
    u0 = &nsst->sol.un[3*(pt->v[0]-1)];
    u1 = &nsst->sol.un[3*(pt->v[1]-1)];
    u2 = &nsst->sol.un[3*(pt->v[2]-1)];
    u3 = &nsst->sol.un[3*(pt->v[3]-1)];
    nsst->sol.u[3*(k-1)+0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0] + cb[3]*u3[0];
    nsst->sol.u[3*(k-1)+1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1] + cb[3]*u3[1];
    nsst->sol.u[3*(k-1)+2] = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2] + cb[3]*u3[2];
  }

  return(1);
}
