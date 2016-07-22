/*
 * main program file for nstokes
 * (C) Copyright 2006 - 2015, ICS-SU
 *
 * This file is part of nstokes.
 *
 * nstokes is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * nstokes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SUstokes.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "nstokes.h"
#include "ns_calls.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n # unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout," abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout," code error...\n");  break;
    case SIGFPE:
      fprintf(stdout," floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout," illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout," segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout," programm killed.\n");  break;
  }
  fprintf(stdout," # no data file saved.\n"); 
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s [+/-v | -h | -N] [-dt step] [-mt val] [-nt n] [-n nit] [-r res] [-t typ] [-ts n] source[.mesh] [-p param[.elas]] [-s data[.sol]] [-o output[.sol]]\n",prog);

  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -N           Navier-Stokes solver\n\
  -dt step     time step (time units)\n\
  -mt val      max time (time units)\n\
  -nt n        number of time steps\n\
  -n nit       number of iterations max for convergence\n\
  -r res       value of the residual (Krylov space) for convergence (default: %e)\n\
  -t typ       specify the type of FE space: 1: P1bP1(*), 2: P2P1\n\
  -ts n        save solution every n time steps\n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\n\
  source.mesh    name of the mesh file\n\
  param.elas     name of file containing elasticity parameters\n\
  data.sol       name of file containing the initial solution or boundary conditions\n\
  output.sol     name of the output file (displacement field)\n",NS_RES);
  exit(1);
}


/* parsing arguments on command line */
static int parsar(int argc,char *argv[],NSst *nsst) {
  int      i;
  char    *ptr,*data;

  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case '-':
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],NS_VER,NS_REL);
          exit(1);
        }
        break;
      case 'h':  /* on-line help */
      case '?':
        usage(argv[0]);
        break;
      case 'N':
        nsst->sol.sim = Navier;
        break;
      case 'd':
        if ( !strcmp(argv[i],"-dt") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            nsst->sol.dt = strtod(argv[i],NULL);
          else {
            fprintf(stderr,"%s: missing argument option\n",argv[0]);
            usage(argv[0]);
          }
        }
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]); 
        }
        break;
      case 'i':
        if ( ++i < argc ) {
          nsst->mesh.name = argv[i];
          ptr = strstr(nsst->mesh.name,".mesh");
          if ( !ptr )  strcat(nsst->mesh.name,".mesh");
        }
        else {
          fprintf(stdout,"%s: missing input file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-mt") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            nsst->sol.mt = strtod(argv[i],NULL);
          else {
            fprintf(stderr,"%s: missing argument option\n",argv[0]);
            usage(argv[0]);
          }
        }
      break;
      case 'n':
        if ( !strcmp(argv[i],"-nt") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            nsst->sol.nt = atoi(argv[i]);
          else {
            fprintf(stderr,"%s: missing argument option\n",argv[0]);
            usage(argv[0]);
          }
        }
        else if ( ++i < argc && isdigit(argv[i][0]) )
          nsst->sol.nit = atoi(argv[i]);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 'o':
        if ( ++i < argc ) {
          nsst->sol.nameout = argv[i];
          ptr = strstr(nsst->sol.nameout,".sol");
          if ( !ptr )  strcat(nsst->sol.nameout,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'p':
        if ( ++i < argc ) {
          nsst->sol.namepar = argv[i];
          ptr = strstr(nsst->sol.namepar,".elas");
          if ( !ptr )  strcat(nsst->sol.namepar,".elas");
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'r':
        if ( ++i < argc && isdigit(argv[i][0]) )
          nsst->sol.res = strtod(argv[i],NULL);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
          break;
      case 's':
        if ( ++i < argc ) {
          nsst->sol.namein = argv[i];
          ptr = strstr(nsst->sol.namein,".sol");
          if ( !ptr )  strcat(nsst->sol.namein,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 't':
        if ( !strcmp(argv[i],"-ts") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            nsst->sol.ts = atoi(argv[i]);
          else {
            fprintf(stderr,"%s: missing argument option\n",argv[0]);
            usage(argv[0]);
          }
        }
        else if ( ++i < argc && isdigit(argv[i][0]) )
          nsst->info.typ = atoi(argv[i]);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          nsst->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          nsst->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( nsst->mesh.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        nsst->mesh.name = data;
      }
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( nsst->mesh.name == NULL ) {
    if ( nsst->info.verb != '0' )  fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }

  /* set time stepping */
  if ( (nsst->sol.dt > 0.0) && (nsst->sol.nt == 0) ) {
    if ( nsst->info.verb != '0' )
      fprintf(stdout," # incorrect time stepping: %d it, %f dt\n",nsst->sol.nt,nsst->sol.dt);
    usage(argv[0]);
  }

  return(1);
}


/* parsing boundary conditions */
static int parsop(NSst *nsst) {
  Cl         *pcl;
  Mat        *pm;
  float       fp1,fp2;
  int         i,j,ncld,npar,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  if ( !nsst->sol.namepar ) {
    strcpy(data,nsst->mesh.name);
    ptr = strstr(data,".mesh");
    if ( ptr )  *ptr = '\0';
    strcat(data,".nstokes");
    in = fopen(data,"r");
    if ( !in ) {
      sprintf(data,"%s","DEFAULT.nstokes");
      in = fopen(data,"r");
    }
  }
  else {
    strcpy(data,nsst->sol.namepar);
    ptr = strstr(data,".nstokes");
    if ( !ptr )  strcat(data,".nstokes");
    in = fopen(data,"r");
  }
  if ( !in ) {
    if ( nsst->info.verb != '0' )  fprintf(stdout," # parameter file %s not found\n",data); 
    return(0);
  }
  if ( nsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* read flow parameters */
  nsst->sol.nbcl = 0;
	npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"dirichlet") ){
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=nsst->sol.nbcl; i<nsst->sol.nbcl+ncld; i++) {
        pcl = &nsst->sol.cl[i];
        pcl->typ = Dirichlet;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);

        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fv",pcl->att) ) {
          if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<nsst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = NS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = NS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = NS_tri;
      }
      nsst->sol.nbcl += ncld;
    }
    /* load forces */
    else if ( !strcmp(data,"neumann") ) {
      fscanf(in,"%d",&ncld);
			npar++;
      for (i=nsst->sol.nbcl; i<nsst->sol.nbcl+ncld; i++) {
        pcl = &nsst->sol.cl[i];
        pcl->typ = Neumann;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);
        
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fnv",pcl->att) ) {
          if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<nsst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        else if ( pcl->att == 'n' )  fscanf(in,"%lf ",&pcl->u[0]);

        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = NS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = NS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = NS_tri;

        /* for the time being: no normal at vertices known */
        if ( (pcl->elt == NS_ver) && (pcl->att == 'n') ) {
          if ( nsst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: [%s] %c\n",buf,pcl->att);
          return(0);
        }
      }
      nsst->sol.nbcl += ncld;
    }
    /* slip condition: normal component vanishes */
    else if ( !strcmp(data,"slip") ) {
      npar++;
			fscanf(in,"%d",&ncld);
      for (i=nsst->sol.nbcl; i<nsst->sol.nbcl+ncld; i++) {
        pcl = &nsst->sol.cl[i];
        pcl->typ = Slip;
        fscanf(in,"%d %s %lf",&pcl->ref,buf,&pcl->u[0]);

        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = NS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = NS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = NS_tri;
      }
      nsst->sol.nbcl += ncld;
    }
    /* surface tension: ref+gamma; atmosph.pressure (only on vertices) */
    else if ( !strcmp(data,"tension") || !strcmp(data,"atmpres") ) {
      npar++;
			fscanf(in,"%d",&ncld);
      for (i=nsst->sol.nbcl; i<nsst->sol.nbcl+ncld; i++) {
        pcl = &nsst->sol.cl[i];
        pcl->typ = !strcmp(data,"tension") ? Tension : AtmPres;
        pcl->elt = NS_ver;
        fscanf(in,"%d %lf",&pcl->ref,&pcl->u[0]);
      }
      nsst->sol.nbcl += ncld;
    }
    /* gravity or body force */
    else if ( !strcmp(data,"gravity") ) {
      npar++;
      nsst->sol.cltyp |= Gravity;
      for (j=0; j<nsst->info.dim; j++)
        fscanf(in,"%lf",&nsst->sol.gr[j]);
    }
    /* flow coefficients */
    else if ( !strcmp(data,"domain") ) {
			npar++;
      fscanf(in,"%d",&ncld);
      assert(ncld <= NS_MAT);
      nsst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &nsst->sol.mat[i];
        fscanf(in,"%d %lf %lf",&pm->ref,&pm->nu,&pm->rho);
      }
    }
  }
  fclose(in);
  
  for (i=0; i<nsst->sol.nbcl; i++) {
    pcl = &nsst->sol.cl[i];
		nsst->sol.cltyp |= pcl->typ;
    nsst->sol.clelt |= pcl->elt;
  }

  if ( (npar > 0) && (nsst->info.verb != '0') )  fprintf(stdout," %d parameters\n",npar);

  return(1);
}


/* parsing boundary conditions */
static int parsdt(NSst *nsst) {
  double   dt;
  FILE    *in;

  /* check for parameter file */
  in = fopen(".dt","r");
  if ( !in )  return(1);
  fscanf(in,"%lf",&dt);
  fclose(in);
 
  /* check value */
  if ( nsst->sol.dt > 0.0 )
    nsst->sol.dt = NS_MIN(nsst->sol.dt,dt);

  return(1);
}


int main(int argc,char **argv) {
  NSst     nsst;
  int      ier,nel;
  char     stim[32];

  memset(&nsst,0,sizeof(NSst));
  tminit(nsst.info.ctim,TIMEMAX);
  chrono(ON,&nsst.info.ctim[0]);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&nsst.mesh,0,sizeof(Mesh));
  memset(&nsst.sol,0,sizeof(Sol));
  nsst.sol.cl  = (Cl*)calloc(NS_CL,sizeof(Cl));
  nsst.sol.mat = (Mat*)calloc(NS_MAT,sizeof(Mat));
  nsst.sol.sim   = Stokes;   /* Navier-Stokes */
  nsst.sol.dt    = -1.0;     /* time stepping */
  nsst.sol.mt    = -1.0;     /* max time      */
  nsst.sol.nt    = 0;        /* steady-state  */
  nsst.sol.ts    = 0;        /* no saving */
  nsst.sol.res   = NS_RES;
  nsst.sol.nit   = NS_MAXIT;

  /* global parameters */
	nsst.info.dim    = 3;
	nsst.info.ver    = 1;	
  nsst.info.verb   = '1';
  nsst.info.zip    = 0;
  nsst.info.typ    = P1;    /* mini element */
  nsst.info.mfree  = 1;

  /* parse command line */
  if ( !parsar(argc,argv,&nsst) )  return(1);

  /* loading data */
	chrono(ON,&nsst.info.ctim[1]);

  if ( nsst.info.verb != '0' ) {
    fprintf(stdout," - NSTOKES, Release %s, %s\n   %s\n\n",NS_VER,NS_REL,NS_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading mesh */
	ier = loadMesh(&nsst);
	if ( ier <= 0 )  return(1);
  nel = nsst.info.dim == 2 ? nsst.info.nt : nsst.info.ne;

  /* parse parameters in file */
  if ( !parsop(&nsst) )  return(1);
  if ( !parsdt(&nsst) )  return(1);

  /* allocating memory */
  if ( !nsst.sol.u ) {
    if ( nsst.info.typ == P1 )
      nsst.sol.u = (double*)calloc(nsst.info.dim*(nsst.info.np+nel),sizeof(double));
	  else {
	    nsst.info.np2 = 3*nsst.info.np;   /* for now an upper bound */
		  nsst.sol.u = (double*)calloc(nsst.info.dim * (nsst.info.np+nsst.info.np2),sizeof(double));
    }
    assert(nsst.sol.u);
    nsst.sol.p = (double*)calloc(nsst.info.np,sizeof(double));
    assert(nsst.sol.p);
  }

  /* loading solution */
  if ( nsst.sol.namein ) {
    ier = loadSol(&nsst);
    if ( !ier )  return(1);
  }

  /* packing mesh if needed */
  if ( nsst.sol.nmat ) {
    ier = nsst.info.dim == 2 ? pack_2d(&nsst) : pack_3d(&nsst);
		if ( ier == 0 ) {
			if ( nsst.info.verb != '0' )  fprintf(stdout," # Packing error.\n");
		  return(1);
		}
	}

  /* setting adjacency */
  if ( (nsst.sol.sim == Navier) || (nsst.sol.cltyp & Tension) )    
	  nsst.info.dim == 2 ? hashel_2d(&nsst) : hashel_3d(&nsst);

  /* counting P1b|P2 nodes */
  nsst.info.dim == 2 ? addnod_2d(&nsst) : addnod_3d(&nsst);

	chrono(OFF,&nsst.info.ctim[1]);
	printim(nsst.info.ctim[1].gdif,stim);
  if ( nsst.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  if ( !nsst.sol.nameout ) {
    nsst.sol.nameout = (char *)calloc(128,sizeof(char));
    assert(nsst.sol.nameout);
    strcpy(nsst.sol.nameout,nsst.mesh.name);
  }

  /* solve */
  chrono(ON,&nsst.info.ctim[2]);
  if ( nsst.info.verb != '0' )
    fprintf(stdout,"\n ** MODULE NSTOKES: %s (%s)\n",NS_VER,nsst.sol.sim == Navier ? "Navier-Stokes" : "Stokes");
	ier = NS_nstokes(&nsst);
  chrono(OFF,&nsst.info.ctim[2]);
  if ( nsst.info.verb != '0' ) {
		printim(nsst.info.ctim[2].gdif,stim);
    if ( ier )
      fprintf(stdout," ** COMPLETED: %s\n\n",stim);
    else
      fprintf(stdout," ** NOT COMPLETED!: %s\n\n",stim);
	}

  /* save file */
  if ( nsst.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&nsst.info.ctim[3]);
  if ( nsst.info.zip && !unpack(&nsst) )  return(1);

  ier = saveSol(&nsst,0);
	if ( !ier )   return(1);
  chrono(OFF,&nsst.info.ctim[3]);
  if ( nsst.info.verb != '0' ) {
    printim(nsst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
	free(nsst.sol.u);
  free(nsst.sol.p);
  free(nsst.sol.cl);
  free(nsst.sol.mat);

  chrono(OFF,&nsst.info.ctim[0]);
  if ( nsst.info.verb != '0' ) {
	  printim(nsst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}

