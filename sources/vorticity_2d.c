#include "nstokes.h"
#include "sparse.h"


/* right hand side for vorticity */
static int rhsFv_2d(NSst *nsst,double *Fv) {
  pTria    pt;
  pPoint   ppt;
  double  *a,*b,*c,*u0,*u1,*u2,dd,x[3],y[3];
  int      k;

  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];

    a = &nsst->mesh.point[pt->v[0]].c[0];
    b = &nsst->mesh.point[pt->v[1]].c[0];
    c = &nsst->mesh.point[pt->v[2]].c[0];

    x[0] = c[0] - b[0]; 
    x[1] = a[0] - c[0]; 
    x[2] = b[0] - a[0]; 

    y[0] = b[1] - c[1];
    y[1] = c[1] - a[1];
    y[2] = a[1] - b[1];

    u0 = &nsst->sol.u[2*(pt->v[0]-1)+0];
    u1 = &nsst->sol.u[2*(pt->v[1]-1)+0];
    u2 = &nsst->sol.u[2*(pt->v[2]-1)+0];

    dd = u0[0]*x[0]+u1[0]*x[1]+u2[0]*x[2] - (u0[1]*y[0]+u1[1]*y[1]+u2[1]*y[2]);
    dd = dd / 6.0;
    
    Fv[pt->v[0]-1] += dd;
    Fv[pt->v[1]-1] += dd;
    Fv[pt->v[2]-1] += dd;
  }

  /* set homogeneous boundary condition */
  for (k=1; k<=nsst->info.np; k++) {
    ppt = &nsst->mesh.point[k];
    if ( ppt->ref > 0 )  Fv[k-1] = 0.0;
  }

  return(1);
}


/* set TGV to diagonal coefficient when Dirichlet */
static void setTGVv_2d(NSst *nsst,pCsr A) {
  pPoint   ppt;
  int      k;

  for (k=1; k<=nsst->info.np; k++) {
    ppt = &nsst->mesh.point[k];
    if ( ppt->ref ) {
      csrSet(A,k-1,k-1,NS_TGV);
    }
  }
}


/* Laplacian matrix (for vorticity computation) */
static int matAv_2d(NSst *nsst,pCsr A) {
  pTria    pt;
  double  *a,*b,*c,val,x[3],y[3],vol,ivo;
  int      i,j,k,ig,jg,nr,nc,nbe;

  /* memory allocation (estimate) */
  nr  = nsst->info.np;
  nc  = nr;
  nbe = 8 * nsst->info.np;
  csrAlloc(A,nr,nc,nbe,CS_UT+CS_SYM);

  /* store values in A */
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];

    /* measure of K */
    a = &nsst->mesh.point[pt->v[0]].c[0]; 
    b = &nsst->mesh.point[pt->v[1]].c[0]; 
    c = &nsst->mesh.point[pt->v[2]].c[0]; 

    x[0] = c[0] - b[0];
    x[1] = a[0] - c[0];
    x[2] = b[0] - a[0];

    y[0] = b[1] - c[1];
    y[1] = c[1] - a[1];
    y[2] = a[1] - b[1];

    vol  = 0.5 * (x[2]*y[1] - y[2]*x[1]);
    if ( vol < NS_EPSD )  vol = NS_EPSD;
    ivo = 1.0 / (4.0*vol);

    /* vertices */
    for (i=0; i<3; i++) {
      ig = pt->v[i]-1;
      for (j=0; j<3; j++) {
        jg  = pt->v[j]-1;

        val = ivo * (y[i]*y[j] + x[i]*x[j]);
	      if ( ig <= jg )  csrPut(A,ig,jg,val);
      }
    }
  }
  setTGVv_2d(nsst,A);
  csrPack(A);

  if ( nsst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/(nr*nc));

  return(1);
}


/* compute vorticity: -Delta psi = \nabla \times u */ 
int vorticity_2d(NSst *nsst) {
  Csr      Av;
  double  *Fv,res;
  int      ier,nit;

  if ( nsst->info.verb != '0' )  fprintf(stdout,"    Computing vorticity\n");

  /* allocating memory (for dylib) */
  nsst->sol.w = (double*)calloc(nsst->info.np,sizeof(double));
  assert(nsst->sol.w);  
  Fv = (double*)calloc(nsst->info.np,sizeof(double));
  assert(Fv);

  ier = matAv_2d(nsst,&Av);
  if ( !ier ) {
    free(nsst->sol.w);
    free(Fv);
    return(ier > 0);
  }

  ier = rhsFv_2d(nsst,Fv);
  if ( !ier ) {
    free(nsst->sol.w);
    free(Fv);
    return(ier > 0);
  }

  /* solve vorticity problem */
  res = nsst->sol.res;
  nit = nsst->sol.nit;
  ier = csrPrecondGrad(&Av,nsst->sol.w,Fv,&res,&nit,0);

  free(Fv);

  return(ier > 0);
}


