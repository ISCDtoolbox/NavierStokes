#include "nstokes.h"
#include "sparse.h"


/* tetrahedron volume */
static inline double volu_3d(double *a,double *b,double *c,double *d) {
  double    bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,dd;

  bx  = b[0] - a[0];
  by  = b[1] - a[1];
  bz  = b[2] - a[2];
  cx  = c[0] - a[0];
  cy  = c[1] - a[1];
  cz  = c[2] - a[2];
  dx  = d[0] - a[0];
  dy  = d[1] - a[1];
  dz  = d[2] - a[2];

  /* test volume */
  vx  = cy*dz - cz*dy;
  vy  = cz*dx - cx*dz;
  vz  = cx*dy - cy*dx; 
  dd  = (bx*vx + by*vy + bz*vz) / 6.0;

  return(dd);
}


/* invert 3x3 non-symmetric matrix */
static int invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < NS_EPSD )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}


/* set TGV to diagonal coefficient when Dirichlet */
static int setTGV_3d(NSst *nsst,pCsr A) {
  pCl      pcl;
  pEdge    pa;
  pPoint   ppt;
  pTria    pt;
  double   *a,*b,*c,val,vol,ivo,x[3],y[3];
  int      k,i,dof,nc;

  if ( nsst->info.verb == '+' )  fprintf(stdout,"     set diagonal coefficients\n");

  /* at vertices */
  if ( nsst->sol.clelt & NS_ver ) {
    nc = 0;
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
      if ( pcl ) {
        csrSet(A,3*(k-1)+0,3*(k-1)+0,NS_TGV);
        csrSet(A,3*(k-1)+1,3*(k-1)+1,NS_TGV);
        csrSet(A,3*(k-1)+2,3*(k-1)+2,NS_TGV);
        nc++;
      }
	  }
    if ( nsst->info.typ == P2 && nsst->info.np2 ) {
      for (k=nsst->info.np+1; k<=nsst->info.np2; k++) {
        ppt = &nsst->mesh.point[k];
        pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
        if ( pcl ) {
          csrSet(A,3*(k-1)+0,3*(k-1)+0,NS_TGV);
          csrSet(A,3*(k-1)+1,3*(k-1)+1,NS_TGV);
          csrSet(A,3*(k-1)+2,3*(k-1)+2,NS_TGV);
          nc++;
        }
      }
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d nodal values\n",nc);
  }
  /* at edge nodes (not meaningful) */
  if ( nsst->sol.clelt & NS_edg ) {
    dof = nsst->info.typ == P1 ? 2 : 3;
    for (k=1; k<=nsst->info.na; k++) {
      pa = &nsst->mesh.edge[k];
      pcl = getCl(&nsst->sol,pa->ref,NS_edg,Dirichlet);
	    if ( pcl ) {
        for (i=0; i<dof; i++) {
          csrSet(A,3*(pa->v[i]-1)+0,3*(pa->v[i]-1)+0,NS_TGV);
          csrSet(A,3*(pa->v[i]-1)+1,3*(pa->v[i]-1)+1,NS_TGV);
          csrSet(A,3*(pa->v[i]-1)+2,3*(pa->v[i]-1)+2,NS_TGV);
        }
      }
    }  
  }
  /* at triangle nodes */
  if ( nsst->sol.clelt & NS_tri ) {
    dof = nsst->info.typ == P1 ? 3 : 6;
    nc  = 0;
    for (k=1; k<=nsst->info.nt; k++) {
      pt  = &nsst->mesh.tria[k];
      pcl = getCl(&nsst->sol,pt->ref,NS_tri,Dirichlet);
      if ( pcl ) {
        for (i=0; i<dof; i++) {
          csrSet(A,3*(pt->v[i]-1)+0,3*(pt->v[i]-1)+0,NS_TGV);
          csrSet(A,3*(pt->v[i]-1)+1,3*(pt->v[i]-1)+1,NS_TGV);
          csrSet(A,3*(pt->v[i]-1)+2,3*(pt->v[i]-1)+2,NS_TGV);
        }
        nc++;
      }
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d element values\n",nc);
  }

  return(1);
}


/* local stiffness matrix */
static int matAe_P1(double *a,double *b,double *c,double *d,double dt,double nu, double rho, double Ae[30][30]) {
  double   vol,m[9],im[9],mm[9][15],Dp[3][5],ph[9][15];
  char     i,j,p,s;
  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                            {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

  Dp[0][0] = 1;  Dp[0][1] = 0;  Dp[0][2] = 0;  Dp[0][3] = -1;  
  Dp[1][0] = 0;  Dp[1][1] = 1;  Dp[1][2] = 0;  Dp[1][3] = -1; 
  Dp[2][0] = 0;  Dp[2][1] = 0;  Dp[2][2] = 1;  Dp[2][3] = -1;

  /* measure of K(a,b,c,d) */
  vol = volu_3d(a,b,c,d);

  /* mm = tB^-1 */
  for (i=0; i<3; i++) {
    m[i+0] = a[i] - d[i];
    m[i+3] = b[i] - d[i];
    m[i+6] = c[i] - d[i];
  }
	if ( !invmatg(m,im) )  return(0);
    
  /* Ae: boucle sur les 5 points de quadrature */
  memset(Ae,0,30*30*sizeof(double));
  for (p=0; p<5; p++) {
    /* Dp for bubble*/
    Dp[0][4] = 256.0*q[p][1]*q[p][2]*(1.0-2.0*q[p][0]-q[p][1]-q[p][2]);
    Dp[1][4] = 256.0*q[p][0]*q[p][2]*(1.0-q[p][0]-2.0*q[p][1]-q[p][2]);
    Dp[2][4] = 256.0*q[p][1]*q[p][0]*(1.0-q[p][0]-q[p][1]-2.0*q[p][2]);
  
    /* unsteady case */
    if ( dt > 0 ) {
      /* ph_basis function */
      memset(ph,0,9*15*sizeof(double));
      ph[0][0] = ph[3][5] = ph[6][10] = q[p][0];  
      ph[0][1] = ph[3][6] = ph[6][11] = q[p][1];  
      ph[0][2] = ph[3][7] = ph[6][12] = q[p][2];  
      ph[0][3] = ph[3][8] = ph[6][13] = 1.0 - q[p][0]-q[p][1]-q[p][2];   
      ph[0][4] = ph[3][9] = ph[6][14] = 256.0 * q[p][0]*q[p][1]*q[p][2]*(1.0-q[p][0]-q[p][1]-q[p][2]);  
    }

    /* mm = (tBt^-1) Dp */
    memset(mm,0,9*15*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<5; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
        mm[i+3][j+5]  = mm[i][j];
        mm[i+6][j+10] = mm[i][j];
      }
    } 
    /* Ae = vol tmm nn */
    for (i=0; i<15; i++) {
      for (j=i; j<15; j++) {
        for (s=0; s<9; s++) {
          Ae[i][j] += w[p]*vol*nu * mm[s][i] * mm[s][j];
				  if ( dt > 0 ) 
						Ae[i][j] += (w[p]*vol*rho * ph[s][i]*ph[s][j]) / dt;
			  }
      }
    }
  }

  return(1);
}


/* local mass matrix */
static int matBe_P1(double *a,double *b,double *c,double *d,double Be[12][10]) {
  double   m[9],im[9],mm[3][5],ph[3][12],Dp[3][5];
  double   vol;
  char     i,j,p,s;
  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                            {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

  Dp[0][0]=1;  Dp[0][1]=0;  Dp[0][2]=0;  Dp[0][3]=-1;  
  Dp[1][0]=0;  Dp[1][1]=1;  Dp[1][2]=0;  Dp[1][3]=-1; 
  Dp[2][0]=0;  Dp[2][1]=0;  Dp[2][2]=1;  Dp[2][3]=-1;

  /* measure of K(a,b,c,d) */
  vol = volu_3d(a,b,c,d);

  /* mm = tB^-1 */
  for (i=0; i<3; i++) {
    m[i+0] = a[i] - d[i];
    m[i+3] = b[i] - d[i];
    m[i+6] = c[i] - d[i];
  }
	if ( !invmatg(m,im) )  return(0);

  memset(Be,0,12*10*sizeof(double));
  for (p=0; p<5; p++) {
	  /* Dp for bubble */
    Dp[0][4] = 256.0*q[p][1]*q[p][2]*(1-2*q[p][0]-q[p][1]-q[p][2]);
    Dp[1][4] = 256.0*q[p][0]*q[p][2]*(1-q[p][0]-2*q[p][1]-q[p][2]);
    Dp[2][4] = 256.0*q[p][1]*q[p][0]*(1-q[p][0]-q[p][1]-2*q[p][2]);

    /* ph_basis function */
    memset(ph,0,3*12*sizeof(double));
    ph[0][0] = ph[1][4] = ph[2][8] = -q[p][0];  
    ph[0][1] = ph[1][5] = ph[2][9] = -q[p][1];  
    ph[0][2] = ph[1][6] = ph[2][10] = -q[p][2];  
    ph[0][3] = ph[1][7] = ph[2][11] = q[p][0]+q[p][1]+q[p][2]-1.0;
    /* mm = (tBt^-1) Dp */
    memset(mm,0,3*5*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<5; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
      }
    }
	  /* Be = area*wp*tmm N mm */
	  for (i=0; i<12; i++) {
	    for (j=0; j<5; j++) {
	      for (s=0; s<3; s++)
	        Be[i][j] += w[p]*vol * ph[s][i] * mm[s][j];
	    }
	  }
  }

	return(1);
}


static int matAe_P2(double *a,double *b,double *c,double *d,double dt,double nu, double rho, double Ae[30][30]) {
	double   m[9],im[9],mm[9][30],Dp[3][10],ph[9][30];
  double   vol;
  char     i,j,p,s;
  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                            {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

  /* mm = tB^-1 */
  for (i=0; i<3; i++) {
    m[i+0] = a[i] - d[i];
    m[i+3] = b[i] - d[i];
    m[i+6] = c[i] - d[i];
  }
	if ( !invmatg(m,im) )  return(0);

  /* measure of K(a,b,c,d) */
  vol = volu_3d(a,b,c,d);

  /* Ae: boucle sur les 5 points de quadrature */
  memset(Ae,0,30*30*sizeof(double));
  for (p=0; p<5; p++) {
    /* Dp for P2*/
    Dp[0][0]=4*q[p][0]-1; Dp[0][1]=0;            Dp[0][2]=0;           Dp[0][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[1][0]=0;           Dp[1][1]=4*q[p][1]-1;  Dp[1][2]=0;           Dp[1][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[2][0]=0;           Dp[2][1]=0;            Dp[2][2]=4*q[p][2]-1; Dp[2][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    
    Dp[0][4]=4*q[p][1];   Dp[0][5]=4*q[p][2];    Dp[0][6]=4*(1-2*q[p][0]-q[p][1]-q[p][2]); 
    Dp[1][4]=4*q[p][0];   Dp[1][5]=0;            Dp[1][6]=-4*q[p][0]; 
    Dp[2][4]=0;           Dp[2][5]=4*q[p][0];    Dp[2][6]=-4*q[p][0];   
  
    Dp[0][7]=0;           Dp[0][8]=-4*q[p][1];                         Dp[0][9]=-4*q[p][2];
    Dp[1][7]=4*q[p][2];   Dp[1][8]=4*(1-q[p][0]-2*q[p][1]-q[p][2]);    Dp[1][9]=-4*q[p][2];
    Dp[2][7]=4*q[p][1];   Dp[2][8]=-4*q[p][1];                         Dp[2][9]=4*(1-q[p][0]-q[p][1]-2*q[p][2]);
		
		/* unsteady case */
    if ( dt > 0.0 ) {
      /* ph_basis function */
      memset(ph,0,9*30*sizeof(double));
      ph[0][0] = ph[3][10] = ph[6][20] = q[p][0]*(2*q[p][0]-1);  
      ph[0][1] = ph[3][11] = ph[6][21] = q[p][1]*(2*q[p][1]-1);  
      ph[0][2] = ph[3][12] = ph[6][22] = q[p][2]*(2*q[p][2]-1);  
      ph[0][3] = ph[3][13] = ph[6][23] = (1.0 - q[p][0]-q[p][1]-q[p][2])*(1.0 - 2*q[p][0]-2*q[p][1]-2*q[p][2]);   
      ph[0][4] = ph[3][14] = ph[6][24] = 4 * q[p][0]*q[p][1];
      ph[0][5] = ph[3][15] = ph[6][25] = 4 * q[p][0]*q[p][2];
      ph[0][6] = ph[3][16] = ph[6][26] = 4 * q[p][0]*(1.0 - q[p][0]-q[p][1]-q[p][2]);
      ph[0][7] = ph[3][17] = ph[6][27] = 4 * q[p][1]*q[p][2]; 
      ph[0][8] = ph[3][18] = ph[6][28] = 4 * q[p][1]*(1.0 - q[p][0]-q[p][1]-q[p][2]);
      ph[0][9] = ph[3][19] = ph[6][29] = 4 * q[p][2]*(1.0 - q[p][0]-q[p][1]-q[p][2]); 
    }  
    /* mm = (tBt^-1) Dp */
    memset(mm,0,9*30*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<10; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
        mm[i+3][j+10] = mm[i][j];
        mm[i+6][j+20] = mm[i][j];
      }
    } 
    /* Ae = vol tmm nn */
    for (i=0; i<30; i++) {
      for (j=i; j<30; j++) {
        for (s=0; s<9; s++)
          Ae[i][j] += w[p]*vol * nu * mm[s][i] * mm[s][j];
          if ( dt > 0.0 ) 
						Ae[i][j] += (w[p]*vol * rho * ph[s][i] * ph[s][j]) / dt;
      }
    }
  }

  return(1);
}


static int matBe_P2(double *a,double *b,double *c,double *d,double Be[12][10]) {
	double   m[9],im[9],mm[3][10],ph[3][12],Dp[3][10];
  double   vol;
  char     i,j,p,s;
  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                            {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

  /* mm = tB^-1 */
  for (i=0; i<3; i++) {
    m[i+0] = a[i] - d[i];
    m[i+3] = b[i] - d[i];
    m[i+6] = c[i] - d[i];
  }
	if ( !invmatg(m,im) )  return(0);

  /* measure of K(a,b,c,d) */
  vol = volu_3d(a,b,c,d);

  memset(Be,0,12*10*sizeof(double));
  for (p=0; p<5; p++) {
	  /* Dp for bubble */
    Dp[0][0]=4*q[p][0]-1; Dp[0][1]=0;            Dp[0][2]=0;           Dp[0][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[1][0]=0;           Dp[1][1]=4*q[p][1]-1;  Dp[1][2]=0;           Dp[1][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[2][0]=0;           Dp[2][1]=0;            Dp[2][2]=4*q[p][2]-1; Dp[2][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    
    Dp[0][4]=4*q[p][1];   Dp[0][5]=4*q[p][2];    Dp[0][6]=4*(1-2*q[p][0]-q[p][1]-q[p][2]); 
    Dp[1][4]=4*q[p][0];   Dp[1][5]=0;            Dp[1][6]=-4*q[p][0]; 
    Dp[2][4]=0;           Dp[2][5]=4*q[p][0];    Dp[2][6]=-4*q[p][0];   
  
    Dp[0][7]=0;           Dp[0][8]=-4*q[p][1];                         Dp[0][9]=-4*q[p][2];
    Dp[1][7]=4*q[p][2];   Dp[1][8]=4*(1-q[p][0]-2*q[p][1]-q[p][2]);    Dp[1][9]=-4*q[p][2];
    Dp[2][7]=4*q[p][1];   Dp[2][8]=-4*q[p][1];                         Dp[2][9]=4*(1-q[p][0]-q[p][1]-2*q[p][2]);
		
    /* ph_basis function */
    memset(ph,0,3*12*sizeof(double));
    ph[0][0] = ph[1][4] = ph[2][8] = -q[p][0];  
    ph[0][1] = ph[1][5] = ph[2][9] = -q[p][1];  
    ph[0][2] = ph[1][6] = ph[2][10] = -q[p][2];  
    ph[0][3] = ph[1][7] = ph[2][11] = q[p][0]+q[p][1]+q[p][2]-1.0;
    /* mm = (tBt^-1) Dp */
    memset(mm,0,3*10*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<10; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
      }
    }
	  /* Be = area*wp*tmm N mm */
	  for (i=0; i<12; i++) {
	    for (j=0; j<10; j++) {
	      for (s=0; s<3; s++)
	        Be[i][j] += w[p]*vol * ph[s][i] * mm[s][j];
	    }
	  }
  }

  return(1);
}


static int matPNe_P1(double *a,double *b,double PNe[4][4]) {
  return(1);
}


static int matPNe_P2(double *a,double *b,double PNe[4][4]) {
  return(1);
}


/*Local friction matrix*/
static int matPFe_P1(double *a,double *b,double PFe[4][4]) {
  return(1);
}


static int matPFe_P2(double *a,double *b,double PFe[4][4]) {
  return(1);
}


/* set slip condition */
static int slipon_3d(NSst *nsst,pCsr A) {
  pTria    pt;
  pCl      pcl;
  double  *a,*b,*c,PNe[4][4],PFe[4][4];
  int      j,k,l,dof,iga,igb,igc,ier;

  dof = nsst->info.typ == P1 ? 3 : 6;

  for (k=1; k<=nsst->info.nt; k++) {
    pt  = &nsst->mesh.tria[k];
    pcl = getCl(&nsst->sol,pt->ref,NS_tri,Slip);
    if ( !pcl )  continue;
    iga = 3*(pt->v[0]-1);
    igb = 3*(pt->v[1]-1);
    igc = 3*(pt->v[2]-1);
    a = &nsst->mesh.point[pt->v[0]].c[0];
    b = &nsst->mesh.point[pt->v[1]].c[0];
    c = &nsst->mesh.point[pt->v[2]].c[0];
  }

  return(1);
}


/* build stiffness matrix */
static int matAB_3d(NSst *nsst,pCsr A,pCsr B) {
  pTetra   pt;
  pCl      pcl;
  double  *a,*b,*c,*d,Ae[30][30],Be[12][10];
  double   rho,nu;
  int      i,j,k,dof,ier,nra,nca,nrb,ncb,nbe,ig,jg;

  /* memory allocation (estimate) */
  if ( nsst->info.typ == P1 ) {
    nra = nca = nsst->info.dim * (nsst->info.np+nsst->info.ne);
    nbe = 20 * nsst->info.ne;
    csrAlloc(A,nra,nca,nbe,CS_UT+CS_SYM);

    nrb = nsst->info.np;
    ncb = nsst->info.dim * (nsst->info.np+nsst->info.ne);
    nbe = 21 * nsst->info.ne;
    csrAlloc(B,nrb,ncb,nbe,0);
    dof = 5;
  }
  else {
    nra = nca = nsst->info.dim * (nsst->info.np+nsst->info.np2);
    nbe = 60 * (nsst->info.np+nsst->info.np2);
    csrAlloc(A,nra,nca,nbe,CS_UT+CS_SYM);

    nrb = nsst->info.np;
    ncb = nsst->info.dim * (nsst->info.np+nsst->info.np2);
    nbe = 180 * nsst->info.np;
    csrAlloc(B,nrb,ncb,nbe,0);
    dof = 10;
  }

  /* stiffness and mass matrices assembly */
  for (k=1; k<=nsst->info.ne; k++) {
    pt = &nsst->mesh.tetra[k];

    /* viscosity / density of fluid domain */
    if ( !getMat(&nsst->sol,pt->ref,&nu,&rho) )  continue;

    /* measure of K */
    a = &nsst->mesh.point[pt->v[0]].c[0]; 
    b = &nsst->mesh.point[pt->v[1]].c[0]; 
    c = &nsst->mesh.point[pt->v[2]].c[0]; 
    d = &nsst->mesh.point[pt->v[3]].c[0];

    /* local stiffness matrix A */
    if ( nsst->info.typ == P1 ) {
      ier  = matAe_P1(a,b,c,d,nsst->sol.dt,nu,rho,Ae);
      ier += matBe_P1(a,b,c,d,Be);
    }
    else {
      ier  = matAe_P2(a,b,c,d,nsst->sol.dt,nu,rho,Ae);
      ier += matBe_P2(a,b,c,d,Be);
    }
		if ( ier < 2 )  continue;

    /* put a(i,j) into global matrix */
    for (i=0; i<3*dof; i++) {
      ig = pt->v[i % dof];
      ig = 3*(ig-1) + (i / dof);
      for (j=i; j<3*dof; j++) {
        if ( fabs(Ae[i][j]) < NS_EPSD )  continue;
        jg = pt->v[j % dof];
        jg = 3*(jg-1) + (j / dof);
        if ( ig < jg )
          csrPut(A,ig,jg,Ae[i][j]);
        else
          csrPut(A,jg,ig,Ae[i][j]);
      }
    }

    /* local mass matrix B */
    for (i=0; i<12; i++) {    /*(12=mesh->dim*4) */
      ig = pt->v[i % 4]-1;
      for (j=0; j<dof; j++) {
        if ( fabs(Be[i][j]) < NS_EPSD )  continue;
        jg = 3*(pt->v[j]-1) + (i / 4);
        csrPut(B,ig,jg,Be[i][j]);
      }
    }
  }

  /* slip boundary condition */
  if ( (nsst->sol.cltyp & Slip) && (nsst->info.nt > 0) ) {
    ier = slipon_3d(nsst,A);
  }
  setTGV_3d(nsst,A);
  csrPack(A);
  csrPack(B);

  if ( nsst->info.verb == '+' ) {
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nra,nca,100.0*A->nbe/(nra*nca));
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nrb,ncb,100.0*B->nbe/(nrb*ncb));
  }

  return(1);
}


/* build right hand side vector and set boundary conds. */
static int rhsF_P1_3d(NSst *nsst,double *F) {
  pTetra   pt;
  pTria    ptt;
  pEdge    pa;
  pPoint   ppt;
  pCl      pcl;
  double  *vp,vol,len,n[3],w[3],*a,*b,*c,*d,kappa,nu,rho;
  int      i,k,nc;

  if ( nsst->info.verb == '+' )  fprintf(stdout,"     gravity and body forces\n");

  /* gravity as external force */
  if ( nsst->sol.cltyp & Gravity ) {
    nc = 0;
    for (k=1; k<=nsst->info.ne; k++) {
      pt = &nsst->mesh.tetra[k];
      getMat(&nsst->sol,pt->ref,&nu,&rho);

      a = &nsst->mesh.point[pt->v[0]].c[0];
      b = &nsst->mesh.point[pt->v[1]].c[0];
      c = &nsst->mesh.point[pt->v[2]].c[0];
      d = &nsst->mesh.point[pt->v[3]].c[0];
      vol = volu_3d(a,b,c,d);
      for (i=0; i<4; i++) {
        F[3*(pt->v[i]-1)+0] += rho*vol * nsst->sol.gr[0] / 4.0;
        F[3*(pt->v[i]-1)+1] += rho*vol * nsst->sol.gr[1] / 4.0;
        F[3*(pt->v[i]-1)+2] += rho*vol * nsst->sol.gr[2] / 4.0;
      }
      /* bubble part */
      F[3*(k-1)+0] += 32.0*rho*vol * nsst->sol.gr[0] / 105.0;
      F[3*(k-1)+1] += 32.0*rho*vol * nsst->sol.gr[1] / 105.0;
      F[3*(k-1)+2] += 32.0*rho*vol * nsst->sol.gr[2] / 105.0;
      nc++;
    }
    if ( nsst->info.verb == '+' )  fprintf(stdout,"     %d gravity values assigned\n",nc);
  }

  /* check for surface tension or atmosph. pressure */
  if ( (nsst->sol.cltyp & Tension) || (nsst->sol.cltyp & AtmPres) ) {
  }
  
  /* nodal boundary conditions */
  if ( nsst->sol.clelt & NS_ver ) {
    nc = 0;
    for (k=1; k<=nsst->info.np; k++) {
      ppt = &nsst->mesh.point[k];
      
      /* Dirichlet conditions */
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Dirichlet);
			if ( !pcl )  continue;
      vp = pcl->att == 'f' ? &nsst->sol.u[3*(k-1)] : &pcl->u[0];
      F[3*(k-1)+0] = NS_TGV * vp[0];
      F[3*(k-1)+1] = NS_TGV * vp[1];
      F[3*(k-1)+2] = NS_TGV * vp[2];
      
      /* slip conditions */
      pcl = getCl(&nsst->sol,ppt->ref,NS_ver,Slip);
      if ( pcl ) {
        fprintf(stdout," # NOT IMPLEMENTED\n");
        exit(1);
      }
      nc++;
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d nodal values\n",nc);
  }

  if ( nsst->sol.clelt & NS_tri ) {
    nc = 0;
    for (k=1; k<=nsst->info.nt; k++) {
      ptt = &nsst->mesh.tria[k];
      /* Dirichlet conditions */
      pcl = getCl(&nsst->sol,ptt->ref,NS_tri,Dirichlet);
			if ( !pcl )  continue;
      for (i=0; i<3; i++) {
        vp = pcl->att == 'f' ? &nsst->sol.u[3*(ptt->v[i]-1)] : &pcl->u[0];
        F[3*(ptt->v[i]-1)+0] = NS_TGV * vp[0];
        F[3*(ptt->v[i]-1)+1] = NS_TGV * vp[1];
        F[3*(ptt->v[i]-1)+2] = NS_TGV * vp[2];
      }
      nc++;
    }
    if ( nsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d element values\n",nc);
  }

  return(1);
}


/* build right hand side vector and set boundary conds. */
static int rhsF_P2_3d(NSst *nsst,double *F) {
  return(1);
}


/* update RHS for unsteady NS: F^n = F + 1/dt u^n(X^n(x)) */
static int rhsFu_3d(NSst *nsst,double *Fk) {
  pTetra   pt;
  double  *a,*b,*c,*d,vol,d1,d2,idt,rho,nu;
  int      i,k,dof;

  dof = nsst->info.typ == P1 ? 4 : 10;
  for (k=1; k<=nsst->info.ne; k++) {
    pt = &nsst->mesh.tetra[k];
    getMat(&nsst->sol,pt->ref,&nu,&rho);

    a = &nsst->mesh.point[pt->v[0]].c[0];
    b = &nsst->mesh.point[pt->v[1]].c[0];
    c = &nsst->mesh.point[pt->v[2]].c[0];
    d = &nsst->mesh.point[pt->v[3]].c[0];
    vol = volu_3d(a,b,c,d);
    d1 = rho * vol / (4.0 * nsst->sol.dt);
    d2 = rho * vol * 32.0 / (105.0 * nsst->sol.dt);
    for (i=0; i<dof; i++) {
      Fk[3*(pt->v[i]-1)+0] += d1*nsst->sol.u[3*(pt->v[i]-1)+0];
      Fk[3*(pt->v[i]-1)+1] += d1*nsst->sol.u[3*(pt->v[i]-1)+1];
      Fk[3*(pt->v[i]-1)+2] += d1*nsst->sol.u[3*(pt->v[i]-1)+2];
    }
    if ( nsst->info.typ == P1 ) {
      /* bubble part */
      Fk[3*(pt->v[4]-1)+0] += d2*nsst->sol.u[3*(pt->v[4]-1)+0];
      Fk[3*(pt->v[4]-1)+1] += d2*nsst->sol.u[3*(pt->v[4]-1)+1];
      Fk[3*(pt->v[4]-1)+2] += d2*nsst->sol.u[3*(pt->v[4]-1)+1];
    }
  }

  return(1);
}


/* 3d Navier-Stokes */
int nstokes1_3d(NSst *nsst) {
  Csr      A,B;
  double  *F,*Fk,res;
  int      ier,it,jt,nit,sz;
  char     verb,stim[32];
  
  /* -- Part I: matrix assembly */
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");

  /* counting P2 nodes (for dylib) */
	if ( nsst->info.typ == P2 && !nsst->info.np2 ) {
		nsst->info.np2 = hashel_3d(nsst);
		if ( nsst->info.np2 == 0 ) {
			fprintf(stdout," # error: no P2 node added.\n");
			return(0);
		}
	}

  /* allocating memory (for dylib) */
  sz = (nsst->info.typ == P1) ? nsst->info.npi+nsst->info.nei : nsst->info.npi+nsst->info.np2;

  if ( !nsst->sol.u ) {
    nsst->sol.u = (double*)calloc(nsst->info.dim*sz,sizeof(double));
    assert(nsst->sol.u);
	  nsst->sol.p = (double*)calloc(nsst->info.npi,sizeof(double));
    assert(nsst->sol.p);
  }

  /* build matrices */
  ier = matAB_3d(nsst,&A,&B);
  if ( !ier )  return(0);
  F = (double*)calloc(nsst->info.dim*sz,sizeof(double));
  assert(F);
  ier = nsst->info.typ == P1 ? rhsF_P1_3d(nsst,F) : rhsF_P2_3d(nsst,F);

  /* -- Part II: solver */
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    Solving linear system:\n");

  /* steady-state */
  if ( nsst->sol.dt < 0.0 ) {
    if ( nsst->info.mfree ) {
		  free(nsst->mesh.tetra);
      if ( nsst->info.nt )  free(nsst->mesh.tria);
      if ( nsst->info.na )  free(nsst->mesh.edge);
      if ( !nsst->info.zip )  free(nsst->mesh.point);
    }
    /* Uzawa solver */
    ier = csrUzawa(&A,&B,nsst->sol.u,nsst->sol.p,F,&nsst->sol.res,&nsst->sol.nit,nsst->info.verb);
  }
  /* unsteady problem */
  else {
	  nsst->sol.un = (double*)calloc(nsst->info.dim*sz,sizeof(double));
	  assert(nsst->sol.un);
    Fk = (double*)calloc(nsst->info.dim*sz,sizeof(double));
    assert(F);
    it = jt = 1;
    do {
      nsst->sol.tim += nsst->sol.dt;
      /* copy solution at time u^n */
      memcpy(nsst->sol.un,nsst->sol.u,nsst->info.dim*sz*sizeof(double));
      
      /* non-linear term: u_t + u\nabla u */
      if ( nsst->sol.sim == Navier ) {
        ier = nsst->info.typ == P1 ? advect_P1_3d(nsst) : advect_P2_3d(nsst);
        if ( !ier )  break;
      }

      /* right-hand side */
      memcpy(Fk,F,nsst->info.dim*sz*sizeof(double));
      ier = rhsFu_3d(nsst,Fk);

      /* Uzawa solver */
      res = nsst->sol.res;
      nit = nsst->sol.nit;
      ier = csrUzawa(&A,&B,nsst->sol.u,nsst->sol.p,Fk,&res,&nit,nsst->info.verb);
      if ( ier < 1 )  break;
      if ( nsst->info.verb != '0' ) {
        fprintf(stdout,"     iteration %d: res=%E, nit=%5d\r",it,res,nit);
        fflush(stdout);
      }
      /* save solution */
      if ( nsst->sol.ts > 0 && it % nsst->sol.ts == 0 ) {
        verb = nsst->info.verb;
        nsst->info.verb = '0';
        saveSol(nsst,jt);
        nsst->info.verb = verb;
        jt++;
      }
    }
    while ( ++it <= nsst->sol.nt );

    /* free mesh structure */
    if ( nsst->info.mfree ) {
		  free(nsst->mesh.tetra);
		  if ( nsst->info.nt )  free(nsst->mesh.tria);
      if ( nsst->info.na )  free(nsst->mesh.edge);
      if ( !nsst->info.zip )  free(nsst->mesh.point);
    }
    free(Fk);
  }
  free(F);

  return(ier > 0);
}


