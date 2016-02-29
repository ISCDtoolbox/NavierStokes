#ifndef __NS_CALLS_H
#define __NS_CALLS_H

enum {Dirichlet=1, Neumann=2, Tension=4, Slip=8, AtmPres=16, Gravity=32};
enum {P1=1, P2};
enum {NS_ver=1,NS_edg=2,NS_tri=4,NS_tet=8};
enum {Stokes=1, Navier=2};

/* data structure */
typedef struct _NSst NSst;

/* prototypes */
NSst *NS_init(int dim, int ver, char typ,char mfree);
int   NS_stop(NSst *nsst);

int   NS_nstokes(NSst *nsst);


#endif
