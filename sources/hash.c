#include "nstokes.h"

#define KTA     7
#define KTB    11
#define KA     31
#define KB     57
#define KC     79

/* data structures */
typedef struct {
  int    min,ind,nxt,elt;
} Cell;

typedef struct {
  Cell  *cell;
  int    nmax,hsiz,hnxt;
} Htab;


/* hash mesh edges for creating P2 nodes */
int hashel_3d(NSst *nsst) {
  return(1);
}


int addnod_3d(NSst *nsst) {
  
  return(1);
}


static int hcode_2d(Mesh *mesh,Htab *ht,int a,int b,int k,int i) {
  Cell     *pc;
  pTria     pt,pt1;
  int       abmin,adj,sum;

  sum = a+b;
  if ( sum >= ht->nmax )  return(0);

  /* check if edge ab stored */
  pc    = &ht->cell[sum];
  abmin = NS_MIN(a,b);
  if ( !pc->min ) {
    pc->min = abmin;
    pc->elt = k;
    pc->ind = i;
    return(1);
  }

  /* analyze linked list */
  pt  = &mesh->tria[k];
  do {
    pt1 = &mesh->tria[pc->elt];
    if ( pc->min == abmin ) {
      adj = pt1->adj[pc->ind];
      if ( !adj ) {
        pt->adj[i]        = 3*pc->elt+pc->ind;
        pt1->adj[pc->ind] = 3*k+i;
      }
      return(1);
    }
    else if ( !pc->nxt ) {
      pc->nxt = ht->hnxt;
      pc      = &ht->cell[ht->hnxt];
      if ( !pc )  return(0);
      pc->min  = abmin;
      pc->elt  = k;
      pc->ind  = i;
      ht->hnxt = pc->nxt;
      pc->nxt  = 0;

      /* check for size overflow */
      if ( !ht->hnxt )  return(0);
      return(1);
    }
    pc = &ht->cell[pc->nxt];
  } while (1);

  return(0);  
}


/* build adjacency table */
int hashel_2d(NSst *nsst) {
  Htab     ht;
  pTria    pt;
  pPoint   ppt;
  int      k,na;
  char     i,i1,i2;

  if ( nsst->info.verb == '+' )  fprintf(stdout,"    Adjacency table: ");

  /* alloc hash */
  ht.nmax = (int)(3.71 * nsst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = 2 * nsst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;

  /* update adjacency */
  na = 0;
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++) {
      i1 = (i+1) % 3;
      i2 = (i+2) % 3;
      if ( !hcode_2d(&nsst->mesh,&ht,pt->v[i1],pt->v[i2],k,i) )  return(0);
      na++;
    }
  }

  /* add seed with point */
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++) {
      if ( !pt->adj[i] )  nsst->mesh.point[pt->v[1]].s = k;
    }
  }
  for (k=1; k<=nsst->info.nt; k++) {
    pt = &nsst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ppt = &nsst->mesh.point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k; 
    }
  }
  free(ht.cell);

  if ( nsst->info.verb == '+' )  fprintf(stdout," %d updated\n",na);

  return(1);  
}


/* add nodes P1b or P2 to tables */
int addnod_2d(NSst *nsst) {
  pTria    pt;
  int      k,na;

  if ( nsst->info.verb == '+' )  fprintf(stdout,"    Adding nodes: ");

  /* store P1b nodes */
  if ( nsst->info.typ == P1 ) {
    na = nsst->info.nt;
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      pt->v[3] = nsst->info.np + k;
    }
  }
  /* create P2 nodes */
  else {
    
  }

  if ( nsst->info.verb == '+' )  fprintf(stdout," %d\n",na);

  return(1);
}