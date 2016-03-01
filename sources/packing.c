#include "nstokes.h"


/* compactify mesh structure */
int pack_3d(NSst *nsst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  pPoint    ppt;
  double    nu,rho,w[3];
  int      *prm,i,k,nf,id;

  /* check if compression needed */
  nf  = 0;
  for (k=1; k<=nsst->info.nei; k++) {
    pe = &nsst->mesh.tetra[k];
    if ( getMat(&nsst->sol,pe->ref,&nu,&rho) ) {
      nf++;
      for (i=0; i<4; i++)  nsst->mesh.point[pe->v[i]].new = pe->v[i];
    }
  }
  if ( nf == nsst->info.nei )  return(-1);

  /* store permutations */
  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  nsst->info.zip = 1;

  /* compress and renum vertices */
  nf = nsst->info.npi;
  k  = 1;
  while ( k <= nf ) {
    if ( nsst->mesh.point[k].new == 0 ) {
      while ( (nsst->mesh.point[nf].new == 0) && (k <= nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&nsst->mesh.point[0],&nsst->mesh.point[nf],sizeof(Point));
        memcpy(&nsst->mesh.point[nf],&nsst->mesh.point[k],sizeof(Point));
        memcpy(&nsst->mesh.point[k],&nsst->mesh.point[0],sizeof(Point));
        /* swap solution too */
        if ( nsst->sol.u ) {
          memcpy(&w,&nsst->sol.u[3*(nf-1)],3*sizeof(double));
          memcpy(&nsst->sol.u[3*(nf-1)],&nsst->sol.u[3*(k-1)],3*sizeof(double));
          memcpy(&nsst->sol.u[3*(k-1)],&w,3*sizeof(double));
        }
        if ( nsst->sol.p ) {
          w[0] = nsst->sol.p[nf-1];
          nsst->sol.p[nf-1] = nsst->sol.p[k-1];
          nsst->sol.p[k-1]  = w[0];
        }
        nsst->mesh.point[k].new  = nf;
        nsst->mesh.point[nf].new = k;
        nf--;
      }
    }
    k++;
  }
  nsst->info.np = nf;

  /* compress and renum tetrahedra */
  prm = (int*)calloc(nsst->info.nei,sizeof(int));
  assert(prm);
  for (k=1; k<=nsst->info.nei; k++) {
    pe = &nsst->mesh.tetra[k];
    for (i=0; i<4; i++)  pe->v[i] = nsst->mesh.point[pe->v[i]].new;
    prm[k] = k;
  }
  nf = nsst->info.nei;
  k  = 1;
  while ( k <= nf ) {
    pe = &nsst->mesh.tetra[k];
    if ( !getMat(&nsst->sol,pe->ref,&nu,&rho) ) {
      do {
        pe = &nsst->mesh.tetra[nf];
        if ( getMat(&nsst->sol,pe->ref,&nu,&rho) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf into k */
      if ( k < nf ) {
        memcpy(&nsst->mesh.tetra[0],&nsst->mesh.tetra[nf],sizeof(Tetra));
        memcpy(&nsst->mesh.tetra[nf],&nsst->mesh.tetra[k],sizeof(Tetra));
        memcpy(&nsst->mesh.tetra[k],&nsst->mesh.tetra[0],sizeof(Tetra));
        prm[k]  = nf;
        prm[nf] = k;
        nf--;
      }
    }
    k++;
  }
  nsst->info.ne = nf;

  /* update adjacency */
  for (k=1; k<=nsst->info.nei; k++) {
    pe = &nsst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      if ( pe->adj[i] > 0 )
        pe->adj[i] = 4*prm[pe->adj[i]/4] + pe->adj[i] % 4;
    }
  }
  /* update simplices */
  for (k=1; k<=nsst->info.npi; k++) {
    ppt = &nsst->mesh.point[k];
    if ( ppt->s > 0 )  ppt->s = prm[ppt->s];
  }
  free(prm);

  /* renum triangles */
  for (k=1; k<=nsst->info.nti; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = nsst->mesh.point[pt->v[i]].new;
  }
  nf = nsst->info.nti;
  k  = 1;
  while ( k <= nf ) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] > nsst->info.np )  break;
    if ( i < 3 ) {
      do {
        pt = &nsst->mesh.tria[nf];
        for (i=0; i<3; i++)
          if ( pt->v[i] > nsst->info.np )  break;
        if ( i == 3 )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&nsst->mesh.tria[0],&nsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&nsst->mesh.tria[nf],&nsst->mesh.tria[k],sizeof(Tria));
        memcpy(&nsst->mesh.tria[k],&nsst->mesh.tria[0],sizeof(Tria));
        nf--;
      }
    }
    k++;
  }
  nsst->info.nt = nf;

  /* renum edges */
  for (k=1; k<=nsst->info.na; k++) {
    pa = &nsst->mesh.edge[k];
    for (i=0; i<2; i++)  pa->v[i] = nsst->mesh.point[pa->v[i]].new;
  }
  nf = nsst->info.nai;
  k  = 1;
  while ( k <= nf ) {
    pa = &nsst->mesh.edge[k];
    if ( (pa->v[0] > nsst->info.np) || (pa->v[1] > nsst->info.np) ) {
      do {
        pa = &nsst->mesh.edge[nf];
        if ( (pa->v[0] > nsst->info.np) || (pa->v[1] > nsst->info.np) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&nsst->mesh.edge[k],&nsst->mesh.edge[nf],sizeof(Edge));
        nf--;
      }
    }
    k++;
  }
  nsst->info.na = nf;

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
  pEdge     pa;
  pPoint    ppt;
  double    nu,rho,w[2];
  int      *prm,i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=nsst->info.nti; k++) {
    pt = &nsst->mesh.tria[k];
    if ( getMat(&nsst->sol,pt->ref,&nu,&rho) ) {
      nf++;
      for (i=0; i<3; i++)  nsst->mesh.point[pt->v[i]].new = pt->v[i];
    }
  }
  if ( nf == nsst->info.nti )  return(-1);

  /* store permutations */
  if ( nsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  nsst->info.zip = 1;

  /* compress and realloc vertices+solution/data */
  nf = nsst->info.npi;
  k  = 1;
  while ( k <= nf ) {
    if ( nsst->mesh.point[k].new == 0 ) {
      while ( (nsst->mesh.point[nf].new == 0) && (k <= nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&nsst->mesh.point[0],&nsst->mesh.point[nf],sizeof(Point));
        memcpy(&nsst->mesh.point[nf],&nsst->mesh.point[k],sizeof(Point));
        memcpy(&nsst->mesh.point[k],&nsst->mesh.point[0],sizeof(Point));
        /* swap solution too */
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
        nsst->mesh.point[k].new  = nf;
        nsst->mesh.point[nf].new = k;
        nf--;
      }
    }
    k++;
  }
  nsst->info.np = nf;

  /* compress and renum triangles (!need to check with getMat) */
  prm = (int*)calloc(nsst->info.nti,sizeof(int));
  assert(prm);
  for (k=1; k<=nsst->info.nti; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = nsst->mesh.point[pt->v[i]].new;
    prm[k] = k;
  }
  nf = nsst->info.nti;
  k  = 1;
  while ( k <= nf ) {
    pt = &nsst->mesh.tria[k];
    if ( !getMat(&nsst->sol,pt->ref,&nu,&rho) ) {
      do {
        pt = &nsst->mesh.tria[nf];
        if ( getMat(&nsst->sol,pt->ref,&nu,&rho) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf into k */
      if ( k < nf ) {
        memcpy(&nsst->mesh.tria[0],&nsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&nsst->mesh.tria[nf],&nsst->mesh.tria[k],sizeof(Tria));
        memcpy(&nsst->mesh.tria[k],&nsst->mesh.tria[0],sizeof(Tria));
        prm[k]  = nf;
        prm[nf] = k;
        nf--;
      }
    }
    k++;
  }
  nsst->info.nt = nf;

  /* update adjacency */
  for (k=1; k<=nsst->info.nti; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++) {
      if ( pt->adj[i] > 0 )
        pt->adj[i] = 3*prm[pt->adj[i]/3] + pt->adj[i] % 3;
    }
  }
  /* update simplices */
  for (k=1; k<=nsst->info.npi; k++) {
    ppt = &nsst->mesh.point[k];
    if ( ppt->s > 0 )  ppt->s = prm[ppt->s];
  }
  free(prm);

  /* compress and renum edges */
  for (k=1; k<=nsst->info.na; k++) {
    pa = &nsst->mesh.edge[k];
    for (i=0; i<2; i++)  pa->v[i] = nsst->mesh.point[pa->v[i]].new;
  }
  nf = nsst->info.nai;
  k  = 1;
  while ( k <= nf ) {
    pa = &nsst->mesh.edge[k];
    if ( (pa->v[0] > nsst->info.np) || (pa->v[1] > nsst->info.np) ) {
      do {
        pa = &nsst->mesh.edge[nf];
        if ( (pa->v[0] > nsst->info.np) || (pa->v[1] > nsst->info.np) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&nsst->mesh.edge[k],&nsst->mesh.edge[nf],sizeof(Edge));
        nf--;
      }
    }
    k++;
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
  for (k=nsst->info.np+1; k<=nsst->info.npi; k++)
    memset(&nsst->sol.u[dim*(k-1)+0],0,dim*sizeof(double));
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

