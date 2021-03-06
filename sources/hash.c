#include "nstokes.h"


/* data structures */
typedef struct {
  int    min,max,nxt,elt,ind;
} Cell;

typedef struct {
  Cell  *cell;
  int    nmax,hsiz,hnxt;
} Htab;


/* check for node na along edge a,b (create if needed) */
static int hpush_2d(Htab *ht,int a,int b,int *na) {
  Cell     *pc;
  pTria     pt;
  int       k,abmin,ia,ib,sum;
  
  sum = a+b;
  if ( sum >= ht->nmax )  return(0);

  /* check if edge ab stored */
  pc = &ht->cell[sum];
  ia = NS_MIN(a,b);
  ib = NS_MAX(a,b);

  /* existing edge: return node */
  if ( pc->min == ia && pc->elt == ib) 
    return(pc->ind);
  else if ( pc->min ) {
    while ( pc->nxt && pc->nxt < ht->nmax ) {
      pc = &ht->cell[pc->nxt];
      if ( pc->min == ia && pc->elt == ib )
        return(pc->ind);
    }
    /* check for new node */
    pc->nxt = ht->hnxt;
    pc = &ht->cell[pc->nxt];
    if ( ++ht->hnxt == ht->nmax ) {
      ht->nmax *= 1.2;
      ht->cell  = (Cell*)realloc(ht->cell,ht->nmax*sizeof(Cell));
      assert(ht->cell);
      for (k=ht->hnxt; k<ht->nmax; k++)  ht->cell[k].nxt = k+1;
    }
  }
  /* create new node */  
  *na = *na + 1;
  pc->min = ia;
  pc->elt = ib;
  pc->ind = *na;
  pc->nxt = 0;

  return(*na);
}


static int hcode_3d(Tetra *tetra,Htab *ht,int a,int b,int c,int k,int i) {
  Cell     *pc;
  pTetra    pt,pt1;
  int       adj,sum,min,max;

  sum = a+b+c;
  if ( sum >= ht->nmax )  return(0);

  /* check if edge ab stored */
  sum = sum % ht->hsiz;
  pc  = &ht->cell[sum];
  min = NS_MIN3(a,b,c);
  max = NS_MAX3(a,b,c);
  if ( !pc->min ) {
    pc->min = min;
    pc->max = max;
    pc->elt = k;
    pc->ind = i;
    return(1);
  }

  /* analyze linked list */
  pt  = &tetra[k];
  do {
    pt1 = &tetra[pc->elt];
    if ( pc->min == min && pc->max == max ) {
      adj = pt1->adj[pc->ind] / 4;
      if ( !adj ) {
        pt->adj[i]        = 4*pc->elt+pc->ind;
        pt1->adj[pc->ind] = 4*k+i;
      }
      return(-1);
    }
    else if ( !pc->nxt ) {
      pc->nxt = ht->hnxt;
      pc      = &ht->cell[ht->hnxt];
      if ( !pc )  return(0);
      pc->min  = min;
      pc->max  = max;
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


/* insert edge a,b and update adjacent triangle */
static int hcode_2d(Tria *tria,Htab *ht,int a,int b,int k,int i) {
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
  pt  = &tria[k];
  do {
    pt1 = &tria[pc->elt];
    if ( pc->min == abmin ) {
      adj = pt1->adj[pc->ind] / 3;
      if ( !adj ) {
        pt->adj[i]        = 3*pc->elt+pc->ind;
        pt1->adj[pc->ind] = 3*k+i;
      }
      return(-1);
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


/* hash mesh edges for creating P2 nodes */
int hashel_3d(NSst *nsst) {
  Htab     ht;
  pTetra   pt;
  pPoint   ppt;
  int      k,nt;
  char     i,i1,i2,i3,ier;

  /* alloc hash (all nodes considered) */
  ht.nmax = (int)(12.71 * nsst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = nsst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;

  /* update */
  nt = 0;
  for (k=1; k<=nsst->info.ne; k++) {
    pt = &nsst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      i1 = (i+1) % 4;
      i2 = (i+2) % 4;
      i3 = (i+3) % 4;
      ier = hcode_3d(nsst->mesh.tetra,&ht,pt->v[i1],pt->v[i2],pt->v[i3],k,i);
      if ( ier == 0 )
        return(0);
      else if ( ier > 0 )
        nt++;
    }
  }

  /* add seed with point */
  for (k=1; k<=nsst->info.ne; k++) {
    pt   = &nsst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      if ( !pt->adj[i] )  nsst->mesh.point[pt->v[i]].s3 = k;
    }
  }
  for (k=1; k<=nsst->info.ne; k++) {
    pt = &nsst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      ppt = &nsst->mesh.point[pt->v[i]];
      if ( !ppt->s3 )  ppt->s3 = k; 
    }
  }
  free(ht.cell);

  return(nt);
}


/* build adjacency table */
int hashel_2d(NSst *nsst) {
  Htab     ht;
  pTria    pt;
  pPoint   ppt;
  int      k,na;
  char     i,i1,i2,ier;

  /* alloc hash (all nodes considered) */
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
      i1  = (i+1) % 3;
      i2  = (i+2) % 3;
      ier = hcode_2d(nsst->mesh.tria,&ht,pt->v[i1],pt->v[i2],k,i);
      if ( ier == 0 )
        return(0);
      else if ( ier > 0 )
        na++;
    }
  }

  /* add seed with point */
  if ( nsst->info.dim == 2 ) {
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      for (i=0; i<3; i++) {
        if ( !pt->adj[i] )  nsst->mesh.point[pt->v[(i+1)%3]].s2 = k;
      }
    }
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      for (i=0; i<3; i++) {
        ppt = &nsst->mesh.point[pt->v[i]];
        if ( !ppt->s2 )  ppt->s2 = k; 
      }
    }
  }
  else {
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      for (i=0; i<3; i++) {
        if ( !pt->adj[i] )  nsst->mesh.point[pt->v[(i+1)%3]].s3 = k;
      }
    }
    for (k=1; k<=nsst->info.nt; k++) {
      pt = &nsst->mesh.tria[k];
      for (i=0; i<3; i++) {
        ppt = &nsst->mesh.point[pt->v[i]];
        if ( !ppt->s3 )  ppt->s3 = k; 
      }
    }
  }
  free(ht.cell);

  return(na);  
}


/* add nodes P1b or P2 to tables */
int addnod_3d(NSst *nsst) {
  Htab     ht;
  pTetra   pt;
  int      i,k,na;
  static int edg[6][2] = {0,1, 0,2, 0,3, 1,2, 1,3, 2,3};

  /* store P1b nodes */
  if ( nsst->info.typ == P1 ) {
    na = nsst->info.ne;
    for (k=1; k<=nsst->info.ne; k++) {
      pt = &nsst->mesh.tetra[k];
      pt->v[4] = nsst->info.np + k;
    }
  }
  /* create P2 nodes */
  else {
    /* alloc hash */
    ht.nmax = (int)(8.2 * nsst->info.np);
    ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
    assert(ht.cell);

    ht.hsiz = 2 * nsst->info.np;
    ht.hnxt = ht.hsiz;
    for (k=ht.hsiz; k<ht.nmax; k++)
      ht.cell[k].nxt = k+1;

    for (k=1; k<=nsst->info.nei; k++) {
      pt = &nsst->mesh.tetra[k];
      for (i=0; i<6; i++) {
        na = hpush_2d(&ht,pt->v[edg[i][0]],pt->v[edg[i][1]],&nsst->info.na);
        pt->v[i+3] = nsst->info.np + na;
      }
    }
    na = nsst->info.na;

    free(ht.cell);
  }

  return(na);
}


/* add nodes P1b or P2 to tables */
int addnod_2d(NSst *nsst) {
  Htab     ht;
  pTria    pt;
  int      i,k,na;

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
    /* alloc hash */
    ht.nmax = (int)(3.2 * nsst->info.np);
    ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
    assert(ht.cell);

    ht.hsiz = 2 * nsst->info.np;
    ht.hnxt = ht.hsiz;
    for (k=ht.hsiz; k<ht.nmax; k++)
      ht.cell[k].nxt = k+1;

    for (k=1; k<=nsst->info.nti; k++) {
      pt = &nsst->mesh.tria[k];
      for (i=0; i<3; i++) {
        na = hpush_2d(&ht,pt->v[i],pt->v[(i+1) % 3],&nsst->info.na);
        pt->v[i+3] = nsst->info.np + na;
      }
    }
    na = nsst->info.na;

    free(ht.cell);
  }

  return(na);
}