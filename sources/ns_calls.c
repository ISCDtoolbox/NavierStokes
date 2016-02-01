#include "nstokes.h"
#include "ns_calls.h"


NSst *NS_init(int dim,int ver,char typ,char mfree) {
	NSst   *nsst;

  /* default values */
  nsst = (NSst*)calloc(1,sizeof(NSst));
  memset(&nsst->mesh,0,sizeof(Mesh));

  /* solution structure */
  memset(&nsst->sol,0,sizeof(Sol));
  nsst->sol.cl  = (Cl*)calloc(NS_CL,sizeof(Cl));
  nsst->sol.mat = (Mat*)calloc(NS_MAT,sizeof(Mat));
  nsst->sol.res = NS_RES;
  nsst->sol.nit = NS_MAXIT;
  nsst->sol.nbcl = 0;
  nsst->sol.nmat = 0;

  /* global parameters */
  nsst->info.dim    = dim;
  nsst->info.ver    = ver;
  nsst->info.verb   = '1';
  nsst->info.zip    = 0;
  nsst->info.typ    = typ;
  nsst->info.mfree  = mfree;

  /* init timer */
  tminit(nsst->info.ctim,TIMEMAX);
  chrono(ON,&nsst->info.ctim[0]);

  return(nsst);
}


/* free global data structure */
int NS_stop(NSst *nsst) {
	char   stim[32];

	/* release memory */
  free(nsst->sol.u);
	free(nsst->sol.cl);
	free(nsst->sol.mat);

  chrono(OFF,&nsst->info.ctim[0]);
  if ( nsst->info.verb != '0' ) {
	  printim(nsst->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s sec.\n",stim);
  }

	return(1);
}

/* set params (facultative): verb= '-|0|+',  zip = 0|1 */
void NS_setPar(NSst *nsst,char verb,int zip) {
	nsst->info.verb = verb;
	nsst->info.zip  = zip;
}

/* handle boundary conditions:
  typ= Dirichlet, Load
  ref= integer
  att= char 'v', 'f', 'n'
  elt= enum NS_ver, NS_edg, NS_tri, NS_tet */
int NS_setBC(NSst *nsst,int typ,int ref,char att,int elt,double *u) {
  Cl    *pcl;
  int    i;

  pcl = &nsst->sol.cl[nsst->sol.nbcl];
  pcl->typ = typ;
  pcl->ref = ref;
  pcl->att = att;
  pcl->elt = elt;

  if ( pcl->typ == Dirichlet ) {
    if ( !strchr("fv",pcl->att) ) {
      fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
  }
  else if ( pcl->typ == Neumann ) {
    if ( !strchr("fnv",pcl->att) ) {
      if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
    if ( (pcl->elt == NS_ver) && (pcl->att == 'n') ) {
      if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: %c\n",pcl->att);
      return(0);
    }
  }

  if ( pcl->att == 'v' ) {
    for (i=0; i<nsst->info.dim; i++)  pcl->u[i] = u[i];
  }
  else if ( pcl->att == 'n' ) {
    pcl->u[0] = u[0];
  }

  if ( nsst->sol.nbcl == NS_CL-1 )  return(0);
  nsst->sol.nbcl++;

  return(1);
}


/* specify gravity value */
void NS_setGra(NSst *nsst, double *gr) {
  int   i;

	nsst->sol.cltyp |= Gravity;
  for (i=0; i<nsst->info.dim; i++)
    nsst->sol.gr[i] = gr[i];
}


/* specify viscosity+density coefficients */
int NS_setCoef(NSst *nsst,int ref,double nu,double rho) {
  Mat    *pm;

  if ( nsst->sol.nmat == NS_MAT-1 )  return(0);
  
  pm = &nsst->sol.mat[nsst->sol.nmat];
  pm->ref  = ref;
  pm->nu   = nu;
  pm->rho  = rho;
  
  nsst->sol.nmat++;

  return(1);
}

/* Construct solution */
int NS_newSol(NSst *nsst) {
  
  nsst->sol.u = (double*)calloc(nsst->info.dim * (nsst->info.np+nsst->info.np2),sizeof(double));
  assert(nsst->sol.u);
  return(1);
}

/* Add element u[dim] to solution at ip */
int NS_addSol(NSst *nsst,int ip,double *u) {
  memcpy(&nsst->sol.u[nsst->info.dim*(ip-1)],u,nsst->info.dim*sizeof(double));
  return(1);
}

/* construct mesh */
int NS_mesh(NSst *nsst,int np,int na,int nt,int ne) {
	int   dof;

	if ( !nsst )  return(0);

	nsst->info.np = np;
	nsst->info.na = na;
	nsst->info.nt = nt;
	nsst->info.ne = ne;

  /* bound on number of nodes */
  dof = nsst->info.typ == P2 ? 10 : 1;
  nsst->mesh.point = (pPoint)calloc(dof*nsst->info.np+1,sizeof(Point));
	assert(nsst->mesh.point);

  if ( nsst->info.na ) {
    nsst->mesh.edge  = (pEdge)calloc(nsst->info.na+1,sizeof(Edge));
    assert(nsst->mesh.edge);
  }
  if ( nsst->info.nt ) {
    nsst->mesh.tria  = (pTria)calloc(nsst->info.nt+1,sizeof(Tria));
    assert(nsst->mesh.tria);
  }
  if ( nsst->info.ne ) {
    nsst->mesh.tetra  = (pTetra)calloc(nsst->info.ne+1,sizeof(Tetra));
    assert(nsst->mesh.tetra);
  }

  return(1);
}

/* insert mesh elements into structure */
int NS_addVer(NSst *nsst,int idx,double *c,int ref) {
  pPoint   ppt;
	int      i;

	assert(idx > 0 && idx <= nsst->info.np);
	ppt = &nsst->mesh.point[idx];
	for (i=0; i<nsst->info.dim; i++)
    ppt->c[i] = c[i];
	ppt->ref = ref;

	return(1);
}

int NS_addEdg(NSst *nsst,int idx,int *v,int ref) {
	pEdge   pe;
	
	assert(idx > 0 && idx <= nsst->info.na);
	pe = &nsst->mesh.edge[idx];
	memcpy(&pe->v[0],&v[0],2*sizeof(int));
  pe->ref = ref;

	return(1);
}

int NS_addTri(NSst *nsst,int idx,int *v,int ref) {
	pTria   pt;

	assert(idx > 0 && idx <= nsst->info.nt);
	pt = &nsst->mesh.tria[idx];
	memcpy(&pt->v[0],&v[0],3*sizeof(int));
  pt->ref = ref;

	return(1);
}

int NS_addTet(NSst *nsst,int idx,int *v,int ref) {
	pTetra   pt;

	assert(idx > 0 && idx <= nsst->info.ne);
	pt = &nsst->mesh.tetra[idx];
	memcpy(&pt->v[0],&v[0],4*sizeof(int));
  pt->ref = ref;

	return(1);
}

/* return mesh header */
void NS_headMesh(NSst *nsst,int *np,int *na,int *nt,int *ne) {
  int k;
	*np = nsst->info.np;
	*na = nsst->info.na;
	*nt = nsst->info.nt;
	*ne = nsst->info.ne;
}


/* initialize solution vector or Dirichlet conditions 
   return: 1 if completion
           0 if no vertex array allocated
          -1 if previous data stored in struct. */
int NS_iniSol(NSst *nsst,double *u) {
  if ( !nsst->info.np )  return(0);

  /* no data already allocated */
  if ( !nsst->sol.u ) {
    nsst->sol.u  = (double*)u;
		return(1);
  }
	/* resolve potential conflict */
	else {
		free(nsst->sol.u);
    nsst->sol.u  = (double*)u;
		return(-1);
	}
}

/* return pointer to solution (Warning: starts at address 0) */
double *NS_getSol(NSst *nsst) {
	return(nsst->sol.u);
}


int NS_nstokes(NSst *nsst) {
  int   i,ier;
  Cl    *pcl;

  for (i=0; i<nsst->sol.nbcl; i++) {
    pcl = &nsst->sol.cl[i];
		nsst->sol.cltyp |= pcl->typ;
    nsst->sol.clelt |= pcl->elt;
  }

  if ( nsst->info.dim == 2)
		ier = nstokes1_2d(nsst);
	else
		ier = nstokes1_3d(nsst);

	return(ier);	
}


