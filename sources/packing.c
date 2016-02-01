#include "nstokes.h"


/* compactify mesh structure */
int pack_3d(NSst *nsst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    nu,rho,w[3];
  int       i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=nsst->info.ne; k++) {
    pe = &nsst->mesh.tetra[k];
    if ( getMat(&nsst->sol,pe->ref,&nu,&rho) ) {
      nf++;
      for (i=0; i<4; i++)  nsst->mesh.point[pe->v[i]].new = pe->v[i];
    }
  }
  if ( nf == nsst->info.ne )  return(-1);

  /* store permutations */
  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  nsst->info.zip = 1;

  /* compress and renum vertices */
  nf = nsst->info.ne;
  k  = 0;
  while ( ++k <= nf ) {
    if ( nsst->mesh.point[k].new == 0 ) {
      while ( (nsst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&nsst->mesh.point[0],&nsst->mesh.point[nf],sizeof(Point));
        memcpy(&nsst->mesh.point[nf],&nsst->mesh.point[k],sizeof(Point));
        memcpy(&nsst->mesh.point[k],&nsst->mesh.point[0],sizeof(Point));
        nsst->mesh.point[k].new  = nf;
        nsst->mesh.point[nf].new = k;
        
        if ( nsst->sol.u ) {
          memcpy(&w,&nsst->sol.u[3*(nf-1)],3*sizeof(double));
          memcpy(&nsst->sol.u[3*(nf-1)],&nsst->sol.u[3*(k-1)],3*sizeof(double));
          memcpy(&nsst->sol.u[3*(k-1)],&w,3*sizeof(double));
        }
      }
      nf++;
    }
  }
  nsst->info.np = nf;

  /* compress and renum tetrahedra */
  nf = nsst->info.ne;
  k  = 0;
  while ( ++k <= nf ) {
    pe = &nsst->mesh.tetra[k];
    if ( !getMat(&nsst->sol,pe->ref,&nu,&rho) ) {
      while ( !getMat(&nsst->sol,nsst->mesh.tetra[nf].ref,&nu,&rho) && (k < nf) )  nf --;
      /* put nf into k */
      memcpy(&nsst->mesh.tetra[k],&nsst->mesh.tetra[nf],sizeof(Tetra));
      nf--;
    }
    for (i=0; i<4; i++)  pe->v[i] = nsst->mesh.point[pe->v[i]].new;
  }
  nsst->info.ne = nf;

  /* renum triangles */
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    pt->v[0] = nsst->mesh.point[pt->v[0]].new;
    pt->v[1] = nsst->mesh.point[pt->v[1]].new;
    pt->v[2] = nsst->mesh.point[pt->v[2]].new;
  }
  
  /* renum edges */
  for (k=1; k<=nsst->info.na; k++) {
    pa = &nsst->mesh.edge[k];
    pa->v[0] = nsst->mesh.point[pa->v[0]].new;
    pa->v[1] = nsst->mesh.point[pa->v[1]].new;
  }

  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",nsst->info.np);
    if ( nsst->info.na )  fprintf(stdout,", %d edges",nsst->info.na);
    if ( nsst->info.nt )  fprintf(stdout,", %d triangles",nsst->info.nt);
    if ( nsst->info.ne )  fprintf(stdout,", %d tetrahedra",nsst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* mesh renumbering and packing */
int pack_2d(NSst *nsst) {
  pTria     pt;
  pEdge     pa,pb;
  double    l,m,w[2];
  int       i,k,nf,id,dof;

  /* check if compression needed */
  nf = 0;
  dof = nsst->info.typ == P1 ? 3 : 6;
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    if ( getMat(&nsst->sol,pt->ref,&l,&m) ) {
      nf++;
      for (i=0; i<dof; i++)  nsst->mesh.point[pt->v[i]].new = pt->v[i];
    }
  }
  if ( nf == nsst->info.nt )  return(-1);

  /* store permutations */
  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  nsst->info.zip = 1;

  /* compress and realloc vertices+solution/data */
  nf = nsst->info.np;
  k  = 0;
  while ( ++k <= nf ) {
    if ( nsst->mesh.point[k].new == 0 ) {
      while ( (nsst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&nsst->mesh.point[0],&nsst->mesh.point[nf],sizeof(Point));
        memcpy(&nsst->mesh.point[nf],&nsst->mesh.point[k],sizeof(Point));
        memcpy(&nsst->mesh.point[k],&nsst->mesh.point[0],sizeof(Point));
        nsst->mesh.point[k].new  = nf;
        nsst->mesh.point[nf].new = k;

        if ( nsst->sol.u ) {
          memcpy(&w,&nsst->sol.u[2*(nf-1)],2*sizeof(double));
          memcpy(&nsst->sol.u[2*(nf-1)],&nsst->sol.u[2*(k-1)],2*sizeof(double));
          memcpy(&nsst->sol.u[2*(k-1)],&w,2*sizeof(double));
        }
        if ( nsst->sol.p ) {
          w[0] = nsst->sol.p[nf-1];
          nsst->sol.p[nf-1] = nsst->sol.p[k-1];
          nsst->sol.p[k-1]  = w[0];
        }
      }
      nf--;
    }
  }
  nsst->info.np = nf;

  /* compress and renum triangles */
  nf = nsst->info.nt;
  k  = 0;
  while ( ++k <= nf ) {
    pt = &nsst->mesh.tria[k];
    if ( !getMat(&nsst->sol,pt->ref,&l,&m) ) {
      while ( !getMat(&nsst->sol,nsst->mesh.tria[nf].ref,&l,&m) && (k < nf) )  nf --;
      /* put nf into k */
      memcpy(&nsst->mesh.tria[k],&nsst->mesh.tria[nf],sizeof(Tria));
      nf--;
    }
    for (i=0; i<dof; i++)  pt->v[i] = nsst->mesh.point[pt->v[i]].new;
  }
  nsst->info.nt = nf;

  /* compress and renum edges */
  for (k=1; k<=nsst->info.na; k++) {
    pa = &nsst->mesh.edge[k];
    for (i=0; i<3; i++)  pa->v[i] = nsst->mesh.point[pa->v[i]].new;
  }
  nf = nsst->info.na;
  k  = 0;
  while ( ++k <= nf ) {
    pa = &nsst->mesh.edge[k];
    if ( (pa->v[0] > nsst->info.np || pa->v[0] == 0) || 
         (pa->v[1] > nsst->info.np || pa->v[1] == 0) ) {
      pb = &nsst->mesh.edge[nf];
      while ( ((pb->v[0] > nsst->info.np || pb->v[0] == 0) || 
               (pb->v[1] > nsst->info.np || pb->v[1] == 0)) && (k < nf) ) {
        nf--;
        pb = &nsst->mesh.edge[nf];
      }
      memcpy(&nsst->mesh.edge[k],&nsst->mesh.edge[nf],sizeof(Edge));
      nf--;
    }
  }
  nsst->info.na = nf;

  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",nsst->info.np);
    if ( nsst->info.na )  fprintf(stdout,", %d edges",nsst->info.na);
    if ( nsst->info.nt )  fprintf(stdout,", %d triangles",nsst->info.nt);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* restore solution at initial vertices */
int unpack(NSst *nsst) {
  pPoint  ppt;
	double  w[3];
  int     k,dim;
	char    i;

  if ( nsst->info.np == nsst->info.npi )  return(1);
  else if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"    Uncompressing data: ");
    fflush(stdout);
  }

  dim = nsst->info.dim;
  for (k=1; k<=nsst->info.np; k++) {
    ppt = &nsst->mesh.point[k];
    if ( ppt->new != k ) {
      memcpy(&w,&nsst->sol.u[dim*(k-1)+0],dim*sizeof(double));
      memcpy(&nsst->sol.u[dim*(k-1)+0],&nsst->sol.u[dim*(ppt->new-1)+0],dim*sizeof(double));
      memcpy(&nsst->sol.u[dim*(ppt->new-1)+0],&w,dim*sizeof(double));
    }
  }

  nsst->info.np = nsst->info.npi;
  nsst->info.na = nsst->info.nai;
  nsst->info.nt = nsst->info.nti;
  nsst->info.ne = nsst->info.nei;

  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",nsst->info.np);
  }

  return(1);
}

