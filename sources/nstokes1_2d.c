#include "nstokes.h"
#include "sparse.h"


/* triangle area */
static inline double area_2d(double *a,double *b,double *c) {
  return(0.5 * ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])));
}


/* set TGV to diagonal coefficient when Dirichlet */
static int setTGV_2d(NSst *nsst,pCsr A) {
  pCl      pcl;
  pEdge    pa;
  pPoint   ppt;
  pTria    pt;
  double   *a,*b,*c,val,vol,ivo,x[3],y[3];
  int      k,i,dof;

  /* at vertices */
  if ( nsst->sol.clelt & NS_ver ) {
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
      if ( pcl ) {
        csrSet(A,2*(k-1)+0,2*(k-1)+0,NS_TGV);
        csrSet(A,2*(k-1)+1,2*(k-1)+1,NS_TGV);
      }
	  }
    if ( nsst->info.typ == P2 && nsst->info.np2 ) {
      for (k=nsst->info.np+1; k<=nsst->info.np2; k++) {
        ppt = &nsst->mesh.point[k];
        pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
        if ( pcl ) {
          csrSet(A,2*(k-1)+0,2*(k-1)+0,NS_TGV);
          csrSet(A,2*(k-1)+1,2*(k-1)+1,NS_TGV);
        }
      }
    }
  }
  /* at edge nodes */
  if ( nsst->sol.clelt & NS_edg ) {
    dof = nsst->info.typ == P1 ? 2 : 3;
    for (k=1; k<=nsst->info.na; k++) {
      pa = &nsst->mesh.edge[k];
      pcl = getCl(&nsst->sol,pa->ref,NS_edg,Dirichlet);
	    if ( pcl ) {
        for (i=0; i<dof; i++) {
          csrSet(A,2*(pa->v[i]-1)+0,2*(pa->v[i]-1)+0,NS_TGV);
          csrSet(A,2*(pa->v[i]-1)+1,2*(pa->v[i]-1)+1,NS_TGV);
        }
      }
    }  
  }
  
  return(1);
}


/* local stiffness matrix */
static int matAe_P1(double *a,double *b,double *c,double dt,double nu,double rho,double Ae[12][12]) {
  double   det,idet,m[2][2],mm[4][8],Dp[2][4],ph[4][8];
  char     i,j,p,s;
  static double w[3] = { 1./3., 1./3., 1./3.};
  static double q[3][2] = { {1./6.,1./6.}, {1./6.,2./3.}, {2./3.,1./6.} };

  Dp[0][0] = 1.0;  Dp[0][1] = 0.0;  Dp[0][2] = -1.0;  
  Dp[1][0] = 0.0;  Dp[1][1] = 1.0;  Dp[1][2] = -1.0;

  /* m = tBT^-1 */
  det  = area_2d(a,b,c);
  idet = 0.5 / det;
	m[0][0] = idet*(b[1]-c[1]);    m[0][1] = idet*(c[1]-a[1]);
	m[1][0] = idet*(c[0]-b[0]);    m[1][1] = idet*(a[0]-c[0]);

  /* Ae:loop on 3 quadrature points */
  memset(Ae,0,12*12*sizeof(double));
  for (p=0; p<3; p++) {
	  /* Dp for bubble */
    Dp[0][3] = 27.0 * q[p][1]*(1.0-2.0*q[p][0]-q[p][1]); 
    Dp[1][3] = 27.0 * q[p][0]*(1.0-q[p][0]-2.0*q[p][1]);

    /* unsteady case */
    if ( dt > 0.0 ) {
      /* ph_basis function */
      memset(ph,0,32*sizeof(double));
      ph[0][0] = ph[2][4] = q[p][0];  
      ph[0][1] = ph[2][5] = q[p][1];  
      ph[0][2] = ph[2][6] = 1.0-q[p][0]-q[p][1];  
      ph[0][3] = ph[2][7] = 27.0 * q[p][0]*q[p][1]*(1.0-q[p][0]-q[p][1]);   
    }
    /* mm = (tBt^-1) Dp */
    memset(mm,0,4*8*sizeof(double));
    for (i=0; i<2; i++) {
      for (j=0; j<4; j++) {
        for (s=0; s<2; s++)
          mm[i][j]  += m[i][s] * Dp[s][j];
        mm[i+2][j+4] = mm[i][j];
      }
    }
	  /* Ae = area*wp*tmm N mm */
	  for (i=0; i<8; i++) {
	    for (j=i; j<8; j++) {
	      for (s=0; s<4; s++) {
	        Ae[i][j] += w[p]*det*nu * mm[s][i] * mm[s][j];
	        if ( dt > 0)
	          Ae[i][j] += (w[p]*det*rho * ph[s][i]*ph[s][j]) / dt;
        }
	    }
	  }
  }

  return(1);
}


/* local mass matrix */
static int matBe_P1(double *a,double *b,double *c,double Be[6][6]) {
  double   Dp[2][4],m[2][2],ph[2][6],mm[2][4];
  double   det,idet;
  char     i,j,p,s;
  static double w[3] = { 1./3., 1./3., 1./3.};
  static double q[3][2] = { {1./6.,1./6.}, {1./6.,2./3.}, {2./3.,1./6.} };

  /* m = tBT^-1 */
	det  = area_2d(a,b,c);
  idet = 0.5 / det;
	m[0][0] = idet*(b[1]-c[1]);    m[0][1] = idet*(c[1]-a[1]);
	m[1][0] = idet*(c[0]-b[0]);    m[1][1] = idet*(a[0]-c[0]);

  Dp[0][0] = 1.0;  Dp[0][1] = 0.0;  Dp[0][2] = -1.0;  
  Dp[1][0] = 0.0;  Dp[1][1] = 1.0;  Dp[1][2] = -1.0;

  memset(Be,0,6*6*sizeof(double));
  for (p=0; p<3; p++) {
	  /* Dp for bubble */
    Dp[0][3] = 27.0 * q[p][1]*(1.0-2.0*q[p][0]-q[p][1]); 
    Dp[1][3] = 27.0 * q[p][0]*(1.0-q[p][0]-2.0*q[p][1]);

    /* ph_basis function */
    memset(ph,0,2*6*sizeof(double));
    ph[0][0] = ph[1][3] = -q[p][0];  
    ph[0][1] = ph[1][4] = -q[p][1];  
    ph[0][2] = ph[1][5] =  q[p][0]+q[p][1]-1.0;  

    /* mm = (tBt^-1) Dp */
    memset(mm,0,2*4*sizeof(double));
    for (i=0; i<2; i++) {
      for (j=0; j<4; j++) {
        for (s=0; s<2; s++)
          mm[i][j] += m[i][s]*Dp[s][j];
      }
    }
	  /* Be = area*wp*tmm N mm */
	  for (i=0; i<6; i++) {
	    for (j=0; j<4; j++) {
	      for (s=0; s<2; s++)
	        Be[i][j] += w[p] * det * ph[s][i] * mm[s][j];
	    }
	  }
  }

  return(1);
}


/* local stiffness matrix */
static int matAe_P2(double *a,double *b,double *c,double dt,double nu, double rho, double Ae[12][12]) {
  double   m[2][2],mi[2][2],mm[4][12],Dp[2][6],ph[4][12];
  double   det,dd,area,x[6],y[6];
  char     i,j,p,s;
  static double w[3] = { 1./3., 1./3., 1./3.};
  static double q[3][2] = { {1./6.,1./6.}, {1./6.,2./3.}, {2./3.,1./6.} };

  /* m = tBT^-1 */
	x[0]=a[0]; x[1]=b[0]; x[2]=c[0]; x[3]=(a[0]+b[0])/2; x[4]=(b[0]+c[0])/2; x[5]=(c[0]+a[0])/2;
	y[0]=a[1]; y[1]=b[1]; y[2]=c[1]; y[3]=(a[1]+b[1])/2; y[4]=(b[1]+c[1])/2; y[5]=(c[1]+a[1])/2;

  /* Ae:loop on 3 quadrature points */
	det = area_2d(a,b,c);
  memset(Ae,0,12*12*sizeof(double));
  for (p=0; p<3; p++) {
	  /* Dp for P2 */           
	  Dp[0][0] = 4*(q[p][0]+q[p][1])-3;    Dp[0][1] = 4*q[p][0]-1;            
		Dp[0][2] = 0;                        Dp[0][3] = 4*(1-2*q[p][0]-q[p][1]);
	  Dp[0][4] = 4*q[p][1];                Dp[0][5] =  -4*q[p][1];
	  Dp[1][0] = 4*(q[p][0]+q[p][1])-3;    Dp[1][1] = 0;
	  Dp[1][2] = 4*q[p][1]-1;              Dp[1][3] = -4*q[p][0];
	  Dp[1][4] = 4*q[p][0];                Dp[1][5] = 4*(1-q[p][0]-2*q[p][1]);
	 
	  m[0][0]  = m[1][0] = m[0][1] = m[1][1] = 0.;
		for (i=0; i<6; i++) {
			m[0][0] += x[i]*Dp[0][i];          m[0][1]  +=x[i]*Dp[1][i];
			m[1][0] += y[i]*Dp[0][i];          
			m[1][1]  +=y[i]*Dp[1][i];
    }
		dd = m[0][0]*m[1][1]-m[0][1]*m[1][0];
		mi[0][0] =  m[1][1]/dd;             mi[0][1] = -m[1][0]/dd;
	 	mi[1][0] = -m[0][1]/dd;             mi[1][1] =  m[0][0]/dd;

    /* unsteady case */
    if ( dt > 0.0 ) {
      /* ph_basis function */
      memset(ph,0,48*sizeof(double));
      ph[0][0] = ph[2][6]  = (1.0-q[p][0]-q[p][1])*(1.0-2.0*q[p][0]-2.0*q[p][1]); 
      ph[0][1] = ph[2][7]  = q[p][0]*(2.0*q[p][0]-1);
      ph[0][2] = ph[2][8]  = q[p][1]*(2.0*q[p][1]-1);
      ph[0][3] = ph[2][9]  = 4.0*q[p][0]*(1.0-q[p][0]-q[p][1]);
      ph[0][4] = ph[2][10] = 4.0*q[p][0]*q[p][1];
      ph[0][5] = ph[2][11] = 4.0*q[p][1]*(1.0-q[p][0]-q[p][1]);
    }
 
    /* mm = (tBt^-1) Dp */
    memset(mm,0,4*12*sizeof(double));
    for (i=0; i<2; i++) {
      for (j=0; j<6; j++) {
        for (s=0; s<2; s++)
          mm[i][j]  += mi[i][s]*Dp[s][j];
        mm[i+2][j+6] = mm[i][j];
      }
    }
	  /* Ae = area*wp*tmm N mm */
	  area = w[p] * det;
	  for (i=0; i<12; i++) {
	    for (j=i; j<12; j++) {
	      for (s=0; s<4; s++) {
	        Ae[i][j] += area * nu * mm[s][i] * mm[s][j];
	        if ( dt > 0.0 )
	          Ae[i][j] += (area*rho * ph[s][i]*ph[s][j]) / dt;;
	      }
	    }
	  }
  }

  return(1);
}


/* local mass matrix */
static int matBe_P2(double *a,double *b,double *c,double Be[6][6]) {
  double   m[2][2],mi[2][2],mm[2][6],Dp[2][6],ph[2][6];
  double   det,dd,area,x[6],y[6];
  char     i,j,p,s;
  static double w[3] = { 1./3., 1./3., 1./3.};
  static double q[3][2] = { {1./6.,1./6.}, {1./6.,2./3.}, {2./3.,1./6.} };

  x[0]=a[0]; x[1]=b[0]; x[2]=c[0]; x[3]=(a[0]+b[0])/2; x[4]=(b[0]+c[0])/2; x[5]=(c[0]+a[0])/2;
  y[0]=a[1]; y[1]=b[1]; y[2]=c[1]; y[3]=(a[1]+b[1])/2; y[4]=(b[1]+c[1])/2; y[5]=(c[1]+a[1])/2;

  /* Be:loop on 3 quadrature points */
	det = area_2d(a,b,c);
  memset(Be,0,6*6*sizeof(double));
  for (p=0; p<3; p++) {
	  /* Dp for P2 */
    Dp[0][0] = 4*(q[p][0]+q[p][1])-3;    Dp[0][1] = 4*q[p][0]-1;            
    Dp[0][2] = 0;                        Dp[0][3] = 4*(1-2*q[p][0]-q[p][1]);
    Dp[0][4] = 4*q[p][1];                Dp[0][5] =  -4*q[p][1];
    Dp[1][0] = 4*(q[p][0]+q[p][1])-3;    Dp[1][1] = 0;
    Dp[1][2] = 4*q[p][1]-1;              Dp[1][3] = -4*q[p][0];
    Dp[1][4] = 4*q[p][0];                Dp[1][5] = 4*(1-q[p][0]-2*q[p][1]);

    m[0][0]  = m[1][0] = m[0][1] = m[1][1] = 0.;
		for (i=0; i<6; i++) {
			m[0][0] += x[i]*Dp[0][i];          m[0][1]  += x[i]*Dp[1][i];
			m[1][0] += y[i]*Dp[0][i];          
			m[1][1]  +=y[i]*Dp[1][i];
    }
		dd = m[0][0]*m[1][1]-m[0][1]*m[1][0];
		mi[0][0] =  m[1][1]/dd;             mi[0][1] = -m[1][0]/dd;
	 	mi[1][0] = -m[0][1]/dd;             mi[1][1] =  m[0][0]/dd;
	
    /* P1 ph_basis function */
    memset(ph,0,2*6*sizeof(double));
    ph[0][0] = ph[1][3] = q[p][0]+q[p][1]-1.0;
	  ph[0][1] = ph[1][4] = -q[p][0];  
	  ph[0][2] = ph[1][5] =  -q[p][1];
	
    /* mm = (tBt^-1) Dp */
    memset(mm,0,2*6*sizeof(double));
    for (i=0; i<2; i++) {
      for (j=0; j<6; j++) {
        for (s=0; s<2; s++)
          mm[i][j] += mi[i][s]*Dp[s][j];
      }
    }
	  /* Be = area*wp*tmm N mm */
	  area = det * w[p];
	  for (i=0; i<6; i++) {
	    for (j=0; j<6; j++) {
	      for (s=0; s<2; s++)
	        Be[i][j] += area * ph[s][i] * mm[s][j];
	    }
	  }
  }

  return(1);
}


/* local slip penalization matrix*/
static int matPNe_P1(double *a,double *b,double PNe[4][4]) {
  double        n[2],dx,dy,d,eps;
  char          i;
  static double w[2] ={1./3., 1./6.};

  /* lenght(a,b) */
  dx  = b[0] - a[0];
  dy  = b[1] - a[1];
  d   = sqrt(dx*dx+dy*dy);
  eps = 1.0 / NS_EPS2;
  /* normal: orthogonal to the edge a,b oriented according to the vector ab */
  n[0] = -dy / d;
  n[1] =  dx / d;

  for (i=0; i<4; i++)
    PNe[i][3-i] = eps*d * w[1] * n[0]*n[1];
  PNe[0][1] = eps*d * w[0] * n[0]*n[1];
  PNe[1][0] = PNe[0][1];
  PNe[2][3] = PNe[3][2] = PNe[0][1];
  for (i=0; i<2; i++)
    PNe[i][i]= eps*d * w[0] * n[i]*n[i];

  for (i=2; i<4; i++)
    PNe[i][i]= eps*d * w[0] * n[i-2]*n[i-2];
  PNe[0][2] = PNe[2][0] = eps*d * w[1] * n[0]*n[0];
  PNe[1][3] = PNe[3][1] = eps*d * w[1] * n[1]*n[1];

  return(1);
}


static int matPNe_P2(double *a,double *b,double PNe[4][4]) {
  return(0);
}


/*Local friction matrix*/
static int matPFe_P1(double *a,double *b,double PFe[4][4]) {
  double        t[2], dx, dy,d;
  static double w[2] ={1./3., 1./6.};
  char          i;

  /* Distance between a and b */
  dx = b[0] - a[0];
  dy = b[1] - a[1];
  d  = sqrt(dx*dx+dy*dy);
  /* Tangent to the edge a, b oriented according to the vector ab */
  t[0] = dx / d;
  t[1] = dy / d;

  for (i=0; i<4; i++)
    PFe[i][3-i]= d*w[1]*t[0]*t[1];
  PFe[0][1]= d*w[0]*t[0]*t[1];
  PFe[1][0]= PFe[0][1];
  PFe[2][3]= PFe[3][2]= PFe[0][1];
  for (i=0; i<2; i++)
    PFe[i][i]= d*w[0]*t[i]*t[i];
  for (i=2; i<4; i++)
    PFe[i][i]= d*w[0]*t[i-2]*t[i-2];
  PFe[0][2]= PFe[2][0] = d*w[1]*t[0]*t[0];
  PFe[1][3]= PFe[3][1] = d*w[1]*t[1]*t[1];

  return(1);
}


static int matPFe_P2(double *a,double *b,double PFe[4][4]) {
  return(0);
}


/* set slip condition */
static int slipon_2d(NSst *nsst,pCsr A) {
  pEdge    pa;
  pCl      pcl;
  double  *a,*b,PNe[4][4],PFe[4][4];
  int      j,k,l,dof,iga,igb;

  dof = nsst->info.typ == P1 ? 2 : 3;

  for (k=1; k<=nsst->info.na; k++) {
    pa  = &nsst->mesh.edge[k];
    pcl = getCl(&nsst->sol,pa->ref,NS_edg,Slip);
    if ( !pcl )  continue;

    if ( pa->v[0] > nsst->info.np || pa->v[1] > nsst->info.np ) continue;
    iga = 2*(pa->v[0]-1);
    igb = 2*(pa->v[1]-1);
    a = &nsst->mesh.point[pa->v[0]].c[0];
    b = &nsst->mesh.point[pa->v[1]].c[0];

    if ( nsst->info.typ == P1 ) {
      matPNe_P1(a,b,PNe);
      matPFe_P1(a,b,PFe);
    }
    else {
      matPNe_P2(a,b,PNe);
      matPFe_P2(a,b,PFe);
    }

    /* P1b P1 */
    for (j=0; j<2; j++) {
      for (l=0; l<2; l++) {
        csrPut(A,iga+j,iga+l,PNe[j][l]);
        csrPut(A,iga+j,igb+l,PNe[j][dof+l]);

        csrPut(A,igb+j,igb+l,PNe[dof+j][dof+l]);
        csrPut(A,igb+j,iga+l,PNe[dof+j][l]);
      }
    }

    /* global friction term */
    if ( fabs(pcl->u[0]) < NS_EPSD )  continue;

    for (j=0; j<2; j++) {
      for (l=0; l<2; l++) {
        csrPut(A,iga+j,iga+l,pcl->u[0] * PFe[j][l]);
        csrPut(A,iga+j,igb+l,pcl->u[0] * PFe[j][dof+l]);

        csrPut(A,igb+j,igb+l,pcl->u[0] * PFe[dof+j][dof+l]);
        csrPut(A,igb+j,iga+l,pcl->u[0] * PFe[dof+j][l]);
      }
    }
  }

  return(1);
}


/* build stiffness matrix */
static int matAB_2d(NSst *nsst,pCsr A,pCsr B) {
  pTria    pt;
  pCl      pcl;
  double  *a,*b,*c,Ae[12][12],Be[6][6];
  double   rho,nu;
  int      i,j,k,dof,ier,nra,nca,nrb,ncb,nbe,ig,jg;

  /* memory allocation (estimate) */
  if ( nsst->info.typ == P1 ) {
    nra = nca = nsst->info.dim * (nsst->info.np+nsst->info.nt);
    nbe = 13 * nsst->info.nt;
    csrAlloc(A,nra,nca,nbe,CS_UT+CS_SYM);

    nrb = nsst->info.np;
    ncb = nsst->info.dim * (nsst->info.np+nsst->info.nt);
    nbe = 14 * nsst->info.nt;
    csrAlloc(B,nrb,ncb,nbe,0);
    dof = 4;
  }
  else {
    nra = nca = nsst->info.dim * (nsst->info.np+nsst->info.np2);
    nbe = 28 * nsst->info.nt;
    csrAlloc(A,nra,nca,nbe,CS_UT+CS_SYM);

    nrb = nsst->info.np;
    ncb = nsst->info.dim * (nsst->info.np+nsst->info.np2);
    nbe = 22 * nsst->info.nt;
    csrAlloc(B,nrb,ncb,nbe,0);
    dof = 6;
  }

  /* stiffness and mass matrices assembly */
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];

    /* viscosity / density of fluid domain */
    if ( !getMat(&nsst->sol,pt->ref,&nu,&rho) )  continue;
 
    /* measure of K */
    a = &nsst->mesh.point[pt->v[0]].c[0]; 
    b = &nsst->mesh.point[pt->v[1]].c[0]; 
    c = &nsst->mesh.point[pt->v[2]].c[0]; 

    /* local stiffness matrix A */
    if ( nsst->info.typ == P1 ) {
      ier  = matAe_P1(a,b,c,nsst->sol.dt,nu,rho,Ae);
      ier += matBe_P1(a,b,c,Be);
    }
    else {
      ier  = matAe_P2(a,b,c,nsst->sol.dt,nu,rho,Ae);
      ier += matBe_P2(a,b,c,Be);
    }
		if ( ier < 2 )  continue;

    /* put a(i,j) into global matrix */
	  for (i=0; i<2*dof; i++) {
	    ig = pt->v[i % dof];
	    ig = 2*(ig-1) + (i / dof);
	    for (j=i; j<2*dof; j++) {
        if ( fabs(Ae[i][j]) < NS_EPSD )  continue;
	      jg = pt->v[j % dof];
	      jg = 2*(jg-1) + (j / dof);
	      if ( ig < jg )
	        csrPut(A,ig,jg,Ae[i][j]);
	      else
	        csrPut(A,jg,ig,Ae[i][j]);
	    }
	  }

    /* local mass matrix B */
  	for (i=0; i<6; i++) {
  	  ig = pt->v[i % 3]-1;
  	  for (j=0; j<dof; j++) {
        if ( fabs(Be[i][j]) < NS_EPSD )  continue;
  	    jg = 2*(pt->v[j]-1) + (i / 3);
        csrPut(B,ig,jg,Be[i][j]);
      }
    }
  }
  
  /* slip boundary condition */
  if ( (nsst->sol.cltyp & Slip) && (nsst->info.na > 0) ) {
    ier = slipon_2d(nsst,A);
  }
  setTGV_2d(nsst,A);
  csrPack(A);
  csrPack(B);

  if ( nsst->info.verb == '+' ) {
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nra,nca,100.0*A->nbe/(nra*nca));
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nrb,ncb,100.0*B->nbe/(nrb*ncb));
  }

  return(1);
}


/* build right hand side vector and set boundary conds. */
static int rhsF_P1_2d(NSst *nsst,double *F) {
  pTria    pt;
  pEdge    pa;
  pPoint   ppt;
  pCl      pcl;
  double  *vp,area,len,n[2],w[2],*a,*b,*c,gamma,kappa,nu,rho;
  int      i,k,ier,nc;
  char     okkap;

  if ( nsst->info.verb == '+' )  fprintf(stdout,"     gravity and body forces\n");

  /* gravity as external force */
  if ( nsst->sol.cltyp & Gravity ) {
    nc = 0;
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      getMat(&nsst->sol,pt->ref,&nu,&rho);

      a = &nsst->mesh.point[pt->v[0]].c[0];
      b = &nsst->mesh.point[pt->v[1]].c[0];
      c = &nsst->mesh.point[pt->v[2]].c[0];
      area = area_2d(a,b,c);
      for (i=0; i<3; i++) {
        F[2*(pt->v[i]-1)+0] += rho*area * nsst->sol.gr[0] / 3.0;
        F[2*(pt->v[i]-1)+1] += rho*area * nsst->sol.gr[1] / 3.0;
      }
      /* bubble part */
      F[2*(k-1)+0] += 9.0*rho*area * nsst->sol.gr[0] / 20.0;
      F[2*(k-1)+1] += 9.0*rho*area * nsst->sol.gr[1] / 20.0;
      nc++;
    }
    if ( nsst->info.verb == '+' )  fprintf(stdout,"     %d gravity values assigned\n",nc);
  }

  /* check for surface tension or atmosph. pressure */
  if ( (nsst->sol.cltyp & Tension) || (nsst->sol.cltyp & AtmPres) ) {
    nc  = 0;
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      if ( ppt->tag & Corner )  continue;

      /* surface tension (check for sign) */
      okkap = 0;
      ier   = 0;
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Tension);
      if ( pcl ) {
        gamma = pcl->u[0];
        ier   = kappa_2d(&nsst->mesh,k,n,&len,&kappa);
        okkap = 1;
        if ( ier ) {
          F[2*(k-1)+0] += gamma * kappa * n[0] * len/2.0;
          F[2*(k-1)+1] += gamma * kappa * n[1] * len/2.0;
          nc++;
        }
      }

      /* atmospheric pressure */
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,AtmPres);
      if ( pcl ) {
        if ( !okkap )  ier = kappa_2d(&nsst->mesh,k,n,&len,&kappa);
        if ( ier ) {
          F[2*(k-1)+0] += pcl->u[0] * n[0] * len/2.0;
          F[2*(k-1)+1] += pcl->u[1] * n[1] * len/2.0;
          nc++;
        }
      }
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d values at interface\n",nc);
  }

  /* nodal boundary conditions */
  if ( nsst->sol.clelt & NS_ver ) {
    nc = 0;
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];

      /* Dirichlet conditions */
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
			if ( pcl ) {
        vp = pcl->att == 'f' ? &nsst->sol.u[2*(k-1)] : &pcl->u[0];
        F[2*(k-1)+0] = NS_TGV * vp[0];
        F[2*(k-1)+1] = NS_TGV * vp[1];
      }
      nc++;
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d nodal values\n",nc);
  }

  return(1);
}


/* build right hand side vector and set boundary conds. */
static int rhsF_P2_2d(NSst *nsst,double *F) {
  return(1);
}


/* update RHS for unsteady NS: F^n = F + 1/dt u^n(X^n(x)) */
static int rhsFu_2d(NSst *nsst,double *Fk) {
  pTria    pt;
  double  *a,*b,*c,area,d1,d2,idt,rho,nu;
  int      i,k,dof;

  dof = nsst->info.typ == P1 ? 3 : 6;
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    getMat(&nsst->sol,pt->ref,&nu,&rho);

    a = &nsst->mesh.point[pt->v[0]].c[0];
    b = &nsst->mesh.point[pt->v[1]].c[0];
    c = &nsst->mesh.point[pt->v[2]].c[0];
    area = area_2d(a,b,c);
    d1 = rho * area / (3.0 * nsst->sol.dt);
    d2 = rho * area * 9.0 / (20.0 * nsst->sol.dt);
    for (i=0; i<dof; i++) {
      Fk[2*(pt->v[i]-1)+0] += d1*nsst->sol.u[2*(pt->v[i]-1)+0];
      Fk[2*(pt->v[i]-1)+1] += d1*nsst->sol.u[2*(pt->v[i]-1)+1];
    }
    if ( nsst->info.typ == P1 ) {
      /* bubble part */
      Fk[2*(pt->v[3]-1)+0] += d2*nsst->sol.u[2*(pt->v[3]-1)+0];
      Fk[2*(pt->v[3]-1)+1] += d2*nsst->sol.u[2*(pt->v[3]-1)+1];
    }
  }

  return(1);
}


/* 2d Navier-Stokes */
int nstokes1_2d(NSst *nsst) {
  Csr      A,B;
  double  *F,*Fk,res;
  int      ier,it,jt,nit,sz;
  char     verb,stim[32];

  /* -- Part I: matrix assembly */
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");

  /* counting P2 nodes (for dylib) */
	if ( nsst->info.typ == P2 && !nsst->info.np2 ) {
		nsst->info.np2 = hashel_2d(nsst);
		if ( nsst->info.np2 == 0 ) {
			fprintf(stdout," # error: no P2 node added.\n");
			return(0);
		}
	}

  /* allocating memory (for dylib) */
  sz = (nsst->info.typ == P1) ? nsst->info.npi+nsst->info.nti : nsst->info.npi+nsst->info.np2;

  if ( !nsst->sol.u ) {
    nsst->sol.u = (double*)calloc(nsst->info.dim*sz,sizeof(double));
    assert(nsst->sol.u);
	  nsst->sol.p = (double*)calloc(nsst->info.npi,sizeof(double));
    assert(nsst->sol.p);
  }

  /* build matrices */
  ier = matAB_2d(nsst,&A,&B);
  if ( !ier )  return(0);
  F = (double*)calloc(nsst->info.dim*sz,sizeof(double));
  assert(F);
  ier = nsst->info.typ == P1 ? rhsF_P1_2d(nsst,F) : rhsF_P2_2d(nsst,F);
  if ( ! ier ) {
    free(nsst->sol.u);
    free(nsst->sol.p);
    free(F);
    return(ier > 0);
  }

  /* -- Part II: solver */
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    Solving linear system:\n");

  /* steady-state */
  if ( nsst->sol.dt < 0.0 ) {
    if ( nsst->info.mfree ) {
		  free(nsst->mesh.tria);
      if ( nsst->info.na )  free(nsst->mesh.edge);
      if ( !nsst->info.zip )  free(nsst->mesh.point);
    }
    /* Uzawa solver */
    ier = csrUzawa(&A,&B,nsst->sol.u,nsst->sol.p,F,&nsst->sol.res,&nsst->sol.nit,nsst->info.verb);
  }
  /* unsteady problem */
  else {
	  /* memory allocation */
    nsst->sol.un = (double*)calloc(nsst->info.dim*sz,sizeof(double));
	  assert(nsst->sol.un);
    Fk = (double*)calloc(nsst->info.dim*sz,sizeof(double));
    assert(Fk);

    it = jt = 1;
    do {
      nsst->sol.tim += nsst->sol.dt;
      /* copy solution at time u^n */
      memcpy(nsst->sol.un,nsst->sol.u,nsst->info.dim*sz*sizeof(double));

      /* non-linear term: u_t + u\nabla u */
      if ( nsst->sol.sim == Navier ) {
        ier = nsst->info.typ == P1 ? advect_P1_2d(nsst) : advect_P2_2d(nsst);
        if ( !ier )  break;
      }
      /* right-hand side */
      memcpy(Fk,F,nsst->info.dim*sz*sizeof(double));
      ier = rhsFu_2d(nsst,Fk);

      /* Uzawa solver */
      res = nsst->sol.res;
      nit = nsst->sol.nit;
      ier = csrUzawa(&A,&B,nsst->sol.u,nsst->sol.p,Fk,&res,&nit,'0');
      if ( nsst->info.verb != '0' ) {
        fprintf(stdout,"     iteration %d: res=%E, nit=%5d, time=%8.4f s\r",it,res,nit,nsst->sol.tim);
        fflush(stdout);
      }
      if ( ier < 1 )  break;

      /* save solution */
      if ( nsst->sol.ts > 0 && it % nsst->sol.ts == 0 ) {
        verb = nsst->info.verb;
        nsst->info.verb = '0';
        saveSol(nsst,jt);
        /* compute vorticity */
        if ( nsst->info.vort > 0 ) {
          ier = NS_vorticity(nsst);
          if ( ier > 0 )  ier = saveVor(nsst,jt);
          if ( ier < 1 )  break;
        }
        nsst->info.verb = verb;
        jt++;
      }
    }
    while ( ++it <= nsst->sol.nt );
    if ( nsst->info.verb != '0' )  fprintf(stdout,"\n");

    /* free mesh structure */
    if ( nsst->info.mfree ) {
		  free(nsst->mesh.tria);
      if ( nsst->info.na )  free(nsst->mesh.edge);
      if ( !nsst->info.zip )  free(nsst->mesh.point);
    }
    free(Fk);
    free(nsst->sol.un);
  }
  free(F);

  return(ier > 0);
}


