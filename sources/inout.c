#include "nstokes.h"
#include "ns_calls.h"
#include "libmesh5.h"


/* read mesh */
int loadMesh(NSst *nsst) {
  pMat       pm;
  pPoint     ppt;
  pEdge      pa;
  pTria      pt1;
  pTetra     pt;
  double    *a,*b,*c,*d,aire,vol;
  float      fp1,fp2,fp3;
  int        i,k,dof,inm,nf,nc,tmp;
  char      *ptr,data[256];

  strcpy(data,nsst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = GmfOpenMesh(data,GmfRead,&nsst->info.ver,&nsst->info.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&nsst->info.ver,&nsst->info.dim)) ) {
        fprintf(stderr," # %s: file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&nsst->info.ver,&nsst->info.dim)) ) {
    fprintf(stderr," # %s: file not found.\n",data);
    return(0);
  }

  if ( nsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  nsst->info.np = GmfStatKwd(inm,GmfVertices);
  nsst->info.na = GmfStatKwd(inm,GmfEdges);
  nsst->info.nt = GmfStatKwd(inm,GmfTriangles);
  nsst->info.ne = GmfStatKwd(inm,GmfTetrahedra);

  if ( !nsst->info.np ) {
    if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }
	nsst->info.npi = nsst->info.np;
	nsst->info.nai = nsst->info.na;
	nsst->info.nti = nsst->info.nt;
	nsst->info.nei = nsst->info.ne;

  /* memory allocation */
  dof = nsst->info.typ == P2 ? 10 : 1;  /* bound on number of nodes */
  nsst->mesh.point = (pPoint)calloc(dof*nsst->info.np+1,sizeof(Point));
  assert(nsst->mesh.point);
  GmfGotoKwd(inm,GmfVertices);
  if ( nsst->info.dim == 2 ) {
    /* 2d mesh vertices */
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      if ( nsst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->ref);
    }
  }
  else {
    /* 3d mesh vertices */
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      if ( nsst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
  }
  /* read mesh corners */
  nc = GmfStatKwd(inm,GmfCorners);
  if ( nc ) {
    GmfGotoKwd(inm,GmfCorners);
    for (k=1; k<=nc; k++) {
      GmfGetLin(inm,GmfCorners,&tmp);
      if ( tmp > 0 && tmp <= nsst->info.npi ) {
        ppt = &nsst->mesh.point[tmp];
        ppt->tag |= Corner;
      }
    }
  }
  /* read mesh edges */
  if ( nsst->info.na > 0 ) {
    nsst->mesh.edge = (pEdge)calloc(nsst->info.na+1,sizeof(Edge));
    assert(nsst->mesh.edge);
    GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=nsst->info.na; k++) {
      pa = &nsst->mesh.edge[k];
      GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
    }
  }
  
  /* read mesh triangles */
  if ( nsst->info.nt > 0 ) {
    nsst->mesh.tria  = (pTria)calloc(nsst->info.nt+1,sizeof(Tria));
    assert(nsst->mesh.tria);
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=nsst->info.nt; k++) {
      pt1 = &nsst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
  }
  /* read mesh tetrahedra */
  if ( nsst->info.ne > 0 ) {
    nsst->mesh.tetra  = (pTetra)calloc(nsst->info.ne+1,sizeof(Tetra));
    assert(nsst->mesh.tetra);
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=nsst->info.ne; k++) {
      pt = &nsst->mesh.tetra[k];
      GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
    }
  }
  GmfCloseMesh(inm);

  if ( nsst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",nsst->info.npi);
    if ( nsst->info.na )  fprintf(stdout,", %d edges",nsst->info.nai);
    if ( nsst->info.nt )  fprintf(stdout,", %d triangles",nsst->info.nti);
    if ( nsst->info.ne )  fprintf(stdout,", %d tetrahedra",nsst->info.nei);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load initial solution */
int loadSol(NSst *nsst) {
  double      bufd[GmfMaxTyp],pmin,pmax,ddp;
  float       fp,buf[GmfMaxTyp];
  int         i,k,dim,ver,np,type,inm,typtab[GmfMaxTyp],offset;
  char       *ptr,data[128];

  if ( !nsst->sol.namein )  return(-1);
  strcpy(data,nsst->sol.namein);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != nsst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
  if ( !np || typtab[0] != GmfVec || np != nsst->info.np )  return(-1);

  if ( nsst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read mesh solutions u_1,..,u_d,p */
  pmin =  DBL_MAX;
  pmax = -DBL_MAX;
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( ver == GmfFloat ) {
    for (k=0; k<nsst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,buf);
      for (i=0; i<nsst->info.dim; i++)
        nsst->sol.u[nsst->info.dim*k+i] = buf[i];
      nsst->sol.p[k] = buf[i];
      if ( nsst->sol.p[k] < pmin )  pmin = nsst->sol.p[k];
      else if ( nsst->sol.p[k] > pmax )  pmax = nsst->sol.p[k];
    }
  }
  else {
    for (k=0; k<nsst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
      for (i=0; i<nsst->info.dim; i++)
        nsst->sol.u[nsst->info.dim*k+i] = bufd[i];
      nsst->sol.p[k] = bufd[i];
      if ( nsst->sol.p[k] < pmin )  pmin = nsst->sol.p[k];
      else if ( nsst->sol.p[k] > pmax )  pmax = nsst->sol.p[k];
    }
  }
  /* rescale pressure field */
  ddp = fabs(pmax - pmin);
  for (k=0; k<nsst->info.np; k++)
    nsst->sol.p[k] /= ddp;

  if ( GmfStatKwd(inm,GmfIterations) ) {
    GmfGotoKwd(inm,GmfIterations);
    GmfGetLin(inm,GmfTime,&np);
  }
  if ( GmfStatKwd(inm,GmfTime) ) {
    GmfGotoKwd(inm,GmfTime);
    if ( nsst->info.ver == GmfFloat ) {
      GmfGetLin(inm,GmfTime,buf);
      nsst->sol.tim = (double)buf[0];
    }
    else {
      GmfGetLin(inm,GmfTime,bufd);
      nsst->sol.tim = bufd[0];
    }
  }
  GmfCloseMesh(inm); 

  if ( nsst->info.verb != '0' ) {
    fprintf(stdout," %d data vectors\n",nsst->info.np);
  }

  return(1);
}


/* save solution */
int saveSol(NSst *nsst,int it) {
  double       dbuf[GmfMaxTyp],pmin,pmax,ddp;
  float        fbuf[GmfMaxTyp],tmpf;
  int          i,k,outm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128],buf[64];

  strcpy(data,nsst->sol.nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    if ( it > 0 ) {
      sprintf(buf,".%d",it);
      strcat(data,buf);
    }
    strcat(data,nsst->info.ver == 1 ? ".solb" : ".sol");
  }
  else {
    ptr = strstr(data,".sol");
    if ( ptr && it > 0 ) {
      *ptr = '\0';
      sprintf(buf,".%d",it);
      strcat(data,buf);
    }
    else if ( !ptr ) {
      if ( it > 0 ) {
        sprintf(buf,".%d",it);
        strcat(data,buf);
      }
      strcat(data,".sol");
    }
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,nsst->info.ver,nsst->info.dim)) ) {
    fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);
  type = 2;
  typtab[0] = GmfVec;
  typtab[1] = GmfSca;

  /* rescale pressure field */
  pmin =  DBL_MAX;
  pmax = -DBL_MAX;
  for (k=0; k<nsst->info.np+nsst->info.np2; k++) {
    if ( nsst->sol.p[k] < pmin )  pmin = nsst->sol.p[k];
    else if ( nsst->sol.p[k] > pmax )  pmax = nsst->sol.p[k];
  }
  ddp = fabs(pmax - pmin);

  /* write P1 sol */
  GmfSetKwd(outm,GmfSolAtVertices,nsst->info.np+nsst->info.np2,type,typtab);
  if ( nsst->info.ver == GmfFloat ) {
    for (k=0; k<nsst->info.np+nsst->info.np2; k++) {
      for (i=0; i<nsst->info.dim; i++)      
        fbuf[i] = nsst->sol.u[nsst->info.dim*k+i];
      fbuf[i] = nsst->sol.p[k] / ddp;
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
  }
  else {
    for (k=0; k<nsst->info.np+nsst->info.np2; k++) {
      for (i=0; i<nsst->info.dim; i++)      
        dbuf[i] = nsst->sol.u[nsst->info.dim*k+i];
      dbuf[i] = nsst->sol.p[k] / ddp;
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }

  /* unsteady case */
  if ( nsst->sol.nt > 1 ) {
    GmfSetKwd(outm,GmfIterations);
    GmfSetLin(outm,GmfIterations,nsst->sol.nt); 
  }
  if ( nsst->sol.dt > 0.0 || nsst->sol.tim > 0.0 ) {
    GmfSetKwd(outm,GmfTime);
    if ( nsst->info.ver == GmfFloat ) {
      tmpf = nsst->sol.tim;
      GmfSetLin(outm,GmfTime,tmpf);
    }
    else
      GmfSetLin(outm,GmfTime,nsst->sol.tim);
  }
  GmfCloseMesh(outm);

  if ( nsst->info.verb != '0' )  fprintf(stdout," %d data vectors\n",nsst->info.np+nsst->info.np2);

  return(1);
}


/* save vorticity solution */
int saveVor(NSst *nsst,int it) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,outm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128],buf[64];

  if ( !nsst->sol.w )  return(0);

  strcpy(data,nsst->sol.nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    if ( it > 0 ) {
      sprintf(buf,".vor.%d",it);
      strcat(data,buf);
    }
    strcat(data,nsst->info.ver == 1 ? ".solb" : ".sol");
  }
  else {
    ptr = strstr(data,".sol");
    if ( ptr && it > 0 ) {
      *ptr = '\0';
      sprintf(buf,".vor.%d",it);
      strcat(data,buf);
    }
    else if ( !ptr ) {
      if ( it > 0 ) {
        sprintf(buf,".vor.%d",it);
        strcat(data,buf);
      }
      strcat(data,".vor.sol");
    }
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,nsst->info.ver,nsst->info.dim)) ) {
    fprintf(stderr," # unable to open vorticity file\n");
    return(0);
  }

  type = 1;
  typtab[0] = GmfSca;
  GmfSetKwd(outm,GmfSolAtVertices,nsst->info.np,type,typtab);

  if ( nsst->info.ver == GmfFloat ) {
    for (k=0; k<nsst->info.np; k++) {
      fbuf[0] = nsst->sol.w[k];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
  }
  else {
    for (k=0; k<nsst->info.np; k++) {
      dbuf[0] = nsst->sol.w[k];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }    
  }
  GmfCloseMesh(outm);

  /* release memory */
  free(nsst->sol.w);
  nsst->sol.w = NULL;

  return(1);
}


