/* command-line help and arguments passing adapted from vcolg.c version 2.0; B D McKay, May 11, 2017 */

#define USAGE \
"gencrn [-n#] [-q] [-z] [-c] [-t|-l|-m] [-x] [-s#] [-b] \n"

#define HELPTEXT \
"  Read digraphs and produce all non-isomorphic CRNs, where the digraphs specify the reaction network. The output format is LBS.\n\
   If the input digraphs are non-isomorphic then the output CRNs are also.\n\
   If the input graph is an undirected, it is considered a reversible network (each undirected edge is turned into two opposite directed edges).\n\
\n\
    -b  encode CRN reactions into 1 byte (#species <= 4) or 2 bytes, and print them as a stream of bytes\n\
    -c  do not produce CRNs containing independent sub-CRNs (e.g. | A->B | C->D) \n\
    -l  only produce CRNs with at least a conservation law\n\
    -m  only produce mass-conserving CRNs (all species are under some conservation law)\n\
    -n# number of species\n\
    -q  only count the number of CRNs\n\
    -s# partition species into classes. For example, \"-s1;2;3\" means that A is in class 1, B and C are in class 2, and E F are in class 3. Classes cardinality must be in non-decreasing order.\n\
    -t  only produce \"dynamically non-trivial\" CRNs (see https://arxiv.org/abs/1705.10820) \n\
    -x  only produce CRNs with no conservation laws\n\
    -z  do not include the naught complex in any reaction (e.g. | -> A)\n"

/*************************************************************************/

#include <math.h>  // for tgamma in fact
#include "gtools.h"
#include "naugroup.h"
#include "nautinv.h"
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#else
#include <math.h>
#endif

nauty_counter vc_nin, vc_nout;
FILE *outfile;

#define MAXNV 128 

static int col[MAXNV];
static booleann first;
static booleann first2;
static int lastreject[MAXNV];
static booleann lastrejok;
static double groupsize;
static unsigned long long newgroupsize;
static int graphSize;
static booleann Tswitch;

static int fail_level;

#define GROUPTEST_NOT 
#ifdef GROUPTEST
static long long totallab;
#endif

/* CS: changes */
#ifndef MAXNN
#define MAXNN 32         /* not more than max(32,WORDSIZE) */
#endif

#if MAXNN > 32
#error "Can't have MAXNN greater than 32"
#endif

#if MAXNN < 32
typedef int xword;   /* Must be as large as MAXNN bits, and
                    must be unsigned if equal to MAXNN bits */
#else
typedef unsigned int xword;
#endif

static unsigned long long int counter = 0;
static unsigned long long int blanksCounter = 0;
static void(*outproc)(FILE*, graph*, int, int, int, xword*);

// outfile, gx, nx, graphSize, spCount, colorMask
static int mindeg, maxdeg, maxn, mine, maxe, mod, res;
#if MAXNN <= 16
static xword xbit[] = { 0x0001,0x0002,0x0004,0x0008,
                   0x0010,0x0020,0x0040,0x0080,
                   0x0100,0x0200,0x0400,0x0800,
                   0x1000,0x2000,0x4000,0x8000 };

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : 15-leftbit[((x)>>8)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] + bytecount[(x)&0xFF])
#elif MAXNN <= 24
static xword xbit[] = { 0x000001,0x000002,0x000004,0x000008,
                   0x000010,0x000020,0x000040,0x000080,
                   0x000100,0x000200,0x000400,0x000800,
                   0x001000,0x002000,0x004000,0x008000,
                   0x010000,0x020000,0x040000,0x080000,
                   0x100000,0x200000,0x400000,0x800000 };

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : \
      (x)&0xFF00 ? 15-leftbit[((x)>>8)&0xFF] : 23-leftbit[((x)>>16)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] \
    + bytecount[((x)>>16)&0xFF] + bytecount[(x)&0xFF])
#else
static xword xbit[] = { 0x00000001,0x00000002,0x00000004,0x00000008,
                   0x00000010,0x00000020,0x00000040,0x00000080,
                   0x00000100,0x00000200,0x00000400,0x00000800,
                   0x00001000,0x00002000,0x00004000,0x00008000,
                   0x00010000,0x00020000,0x00040000,0x00080000,
                   0x00100000,0x00200000,0x00400000,0x00800000,
                   0x01000000,0x02000000,0x04000000,0x08000000,
                   0x10000000,0x20000000,0x40000000,0x80000000 };

#define XNEXTBIT(x) \
    ((x)&0xFF ? 7-leftbit[(x)&0xFF] : \
      (x)&0xFF00 ? 15-leftbit[((x)>>8)&0xFF] : \
    (x)&0xFF0000 ? 23-leftbit[((x)>>16)&0xFF] : \
                      31-leftbit[((x)>>24)&0xFF])
#define XPOPCOUNT(x) (bytecount[((x)>>8)&0xFF] \
    + bytecount[((x)>>16)&0xFF] + \
    + bytecount[((x)>>24)&0xFF] + bytecount[(x)&0xFF])
#endif
static int connec = 0;              /* 1 for -c, 2 for -C, 0 for neither */
statsblk nauty_stats;

static int leastNBitsOn = 0;
static int leastNBitsOnClass = 0;
static int maxColorEdges = 0; // total number of edges necessary to color a CRN

typedef struct
{
  int ne, dmax;            /* values used for xlb,xub calculation */
  int xlb, xub;            /* saved bounds on extension degree */
  xword lo, hi;            /* work purposes for orbit calculation */
  xword xstart[MAXNN + 1]; /* index into xset[] for each cardinality */
  xword *xset;             /* array of all x-sets in card order */
  xword *xcard;            /* cardinalities of all x-sets */
  xword *xinv;             /* map from x-set to index in xset */
  xword *xorb;             /* min orbit representative */
  xword *xx;               /* (-b, -t, -s, -m) candidate x-sets */
                           /*   note: can be the same as xcard */
  xword xlim;              /* number of x-sets in xx[] */
} leveldata;

static leveldata data;      /* data[n] is data for n -> n+1 */
static leveldata classData;

void crnextend(graph* g, int n, int ne, booleann isClassEnumeration); // declaration of crnextend to avoid cyclic dependencies in accept2

#ifndef  MAXN  /* maximum allowed n value; use 0 for dynamic sizing. */
#define MAXN 0
#define MAXMM 0
#else
#define MAXMM ((MAXNN+WORDSIZE-1)/WORDSIZE)  /* max setwords in a set */
#endif  /* MAXN */


static booleann connectedSwitch   = FALSE
       , nonTrivialSwitch         = FALSE
       , printCrnSwitch           = FALSE
       , conservationLawSwitch    = FALSE
       , nonConservationLawSwitch = FALSE
       , noZeroNodeSwitch         = FALSE
       , massConservingSwitch     = FALSE
       , byteEncodingSwitch       = FALSE
       , speciesClassesSwitch     = FALSE;




static int speciesCount;
static int maxComplexes;
static int* colorLab;
static int* colorPtn;
static int  colorMask[4];
static int currMol[MAXNN];
static int maxMol[MAXNN];
static int* speciesClass = NULL;
static int classCount = 0;
static int classCardTotal = 0;
static long speciesClasses2[50];
static int  speciesClasses[50];

/*************************************/
booleann canonise = 0;
static char* speciesNames = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static graph gcan[MAXNN];
static xword pxx[MAXNN];         // TODO: is 12 the max number of species supported?
static xword pxsetStreak[MAXNN]; // holds the current species assignments with the same degree;
                                 // the xsets are stored in Little-Endian order, as opposed to Nauty's graph format,
                                 // which is similar to Big-Endian. 
                                 // This variable helps to ensure that edges are added in canonical order in myFun.
                                 // The encoding must be Little-Endian, because data[n].xorb is encoded in Little-Endian too.

static xword imax;
static xword classMax;


/* copied from naugroup.c, in order to allow multiple copies of *group */
static grouprec *mygroup = NULL;
static int mygroup_depth = 0;
DYNALLSTAT(cosetrec, mycoset, mycoset_sz);

static permrec *mygens;
static permrec *myfreelist = NULL;
static int myfreelist_n = 0;
DYNALLSTAT(int, myallp, myallp_sz);
DYNALLSTAT(int, myid, myid_sz);

struct naugroupState {
  grouprec* group;
  int group_depth;

  // coset representatives; must be populated by makecosetreps(group);
  cosetrec* coset;
  size_t coset_sz;

  permrec *gens;
  permrec *freelist;
  int     freelist_n;
  int     *allp;
  size_t  allp_sz;
  int     *id;
  size_t  id_sz;
};

typedef struct naugroupState naugroupState;





struct crnextendState {
  xword assigned;
  xword assignedHeterodimers;
  int streak;
  booleann isRigid;
  int currentCardinality;
  int currentXSetIndex;
  naugroupState groupState;
};

typedef struct crnextendState crnextendState;

crnextendState states[MAXNN];




/******************************/
/** Group loading operations **/
/******************************/
void saveGroupState(int n) {
  // naugroupState *state = malloc (sizeof(naugroupState));

  states[n].groupState.group       = mygroup;
  states[n].groupState.group_depth = mygroup_depth;
  states[n].groupState.coset       = mycoset;
  states[n].groupState.coset_sz    = mycoset_sz;
  states[n].groupState.gens        = mygens;
  states[n].groupState.freelist    = myfreelist;
  states[n].groupState.freelist_n  = myfreelist_n;
  states[n].groupState.allp        = myallp;
  states[n].groupState.allp_sz     = myallp_sz;
  states[n].groupState.id          = myid;
  states[n].groupState.id_sz       = myid_sz;

  // return state;
}

void loadGroupState(int n) {
  mygroup       = states[n].groupState.group;
  mygroup_depth = states[n].groupState.group_depth;
  mycoset       = states[n].groupState.coset;
  mycoset_sz    = states[n].groupState.coset_sz;
  mygens        = states[n].groupState.gens;
  myfreelist    = states[n].groupState.freelist;
  myfreelist_n  = states[n].groupState.freelist_n;
  myallp        = states[n].groupState.allp;
  myallp_sz     = states[n].groupState.allp_sz;
  myid          = states[n].groupState.id;
  myid_sz       = states[n].groupState.id_sz;

  states[n].groupState.group       = NULL;
  states[n].groupState.group_depth = 0;
  states[n].groupState.coset       = NULL;
  states[n].groupState.coset_sz    = 0;
  states[n].groupState.gens        = NULL;
  states[n].groupState.freelist    = NULL;
  states[n].groupState.freelist_n  = 0;
  states[n].groupState.allp        = NULL;
  states[n].groupState.allp_sz     = 0;
  states[n].groupState.id          = NULL;
  states[n].groupState.id_sz       = 0;
}

void resetGroupState() {
  mygroup       = NULL;
  mygroup_depth = 0;
  mycoset       = 0;
  mycoset_sz    = 0;
  mygens        = NULL;
  myfreelist    = NULL;
  myfreelist_n  = 0;
  myallp        = NULL;
  myallp_sz     = 0;
  myid          = NULL;
  myid_sz       = 0;
}


/*****************************/
/** Adapted from naugroup.c **/
/*****************************/

permrec
*mynewpermrec(int n)
/* Get a permrec of order n.  This procedure and the next one are
designed to be efficient if lots of group ops are done with the
same value of n. */
{
  permrec *p;

  if (myfreelist_n != n)
  {
    while (myfreelist != NULL)
    {
      p = myfreelist;
      myfreelist = myfreelist->ptr;
      free(p);
    }
    myfreelist_n = n;
  }

  if (myfreelist != NULL)
  {
    p = myfreelist;
    myfreelist = myfreelist->ptr;
    return p;
  }

  p = (permrec*)malloc(sizeof(permrec) + (myfreelist_n - 2) * sizeof(int));

  if (p == NULL)
  {
    fprintf(ERRFILE, ">E malloc failed in newpermrec()\n");
    exit(1);
  }

  return p;
}

void
myfreepermrec(permrec *p, int n)
/* Free a permrec of given size. */
{
  permrec *q;

  if (p == NULL) return;

  if (myfreelist_n != n)
  {
    while (myfreelist)
    {
      q = myfreelist;
      myfreelist = myfreelist->ptr;
      free(q);
    }
    myfreelist_n = n;
  }

  p->ptr = myfreelist;
  myfreelist = p;
}

static void
mygroupelts3(levelrec *lr, int n, int level,
  void(*action)(int*, int, int*, void*), int *before,
  int *after, int *id, int *abort, void *userptr)
  /* Recursive routine used by allgroup3. */
{
  int i, j, orbsize;
  int *p, *cr;
  cosetrec *coset;

  coset = lr[level].replist;
  orbsize = lr[level].orbitsize;

  for (j = 0; j < orbsize; ++j)
  {
    cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
    if (before == NULL)
      p = cr;
    else if (cr == NULL)
      p = before;
    else
    {
      p = after;
      for (i = 0; i < n; ++i) p[i] = cr[before[i]];
    }

    if (level == 0)
      (*action)((p == NULL ? id : p), n, abort, userptr);
    else
      mygroupelts3(lr, n, level - 1, action, p, after + n, id, abort, userptr);
    if (*abort) return;
  }
}

DYNALLSTAT(int, queue, queue_sz);
DYNALLSTAT(int, lab, lab_sz);

void
mymakecosetreps(grouprec *grp)
/* Make all coset representatives for this group */
{
  int i, j, k, n, depth;
  int l, index;
  int *p, *q;
  permrec *gen, *g;
  cosetrec *cr;
  int head, tail;

  n = grp->n;
  depth = grp->depth;

  DYNALLOC1(int, queue, queue_sz, n, "malloc");
  DYNALLOC1(int, lab, lab_sz, n, "malloc");

  j = 0;
  for (i = 0; i < depth; ++i)
    j += grp->levelinfo[i].orbitsize;

  if (j > 0) DYNALLOC1(cosetrec, mycoset, mycoset_sz, j, "malloc");

  cr = mycoset;
  for (i = 0; i < depth; ++i)
  {
    grp->levelinfo[i].replist = cr;
    cr += grp->levelinfo[i].orbitsize;
  }

  for (i = 0; i < depth; ++i)
  {
    cr = grp->levelinfo[i].replist;
    gen = grp->levelinfo[i].gens;
    for (j = 0; j < n; ++j) lab[j] = -1;
    queue[0] = grp->levelinfo[i].fixedpt;
    lab[queue[0]] = 0;
    cr[0].image = queue[0];
    cr[0].rep = NULL;
    head = 0;
    tail = 1;
    index = 0;
    while (head < tail)
    {
      j = queue[head++];
      p = (cr[lab[j]].rep ? cr[lab[j]].rep->p : NULL);
      for (g = gen; g != NULL; g = g->ptr)
      {
        k = g->p[j];
        if (lab[k] < 0)
        {
          ++index;
          lab[k] = index;
          queue[tail++] = k;
          cr[index].image = k;
          cr[index].rep = mynewpermrec(n);
          q = cr[index].rep->p;
          if (p == NULL)
            for (l = 0; l < n; ++l) q[l] = g->p[l];
          else
            for (l = 0; l < n; ++l) q[l] = g->p[p[l]];
        }
      }
    }
  }
}


/**************************************************************************/

int
myallgroup3(grouprec *grp, void(*action)(int*, int, int*, void*), void *userptr)
/* Call action(p,n,&abort,userptr) for every element of the group,
   including the identity.  The identity is always the first call.
   If action() stores a non-zero value in abort, group generation is
   aborted and the abort value is returned by this procedure.  If no
   non-zero value is ever returned in abort by action(), this
   procedure returns 0. The pointer userptr is not interpretted and
   is passed to action() to use as it likes. */
{
  int i, depth, n, abort;

  depth = grp->depth;
  n = grp->n;

  DYNALLOC1(int, myid, myid_sz, n, "malloc");
  for (i = 0; i < n; ++i) myid[i] = i;

  abort = 0;
  if (depth == 0)
  {
    (*action)(myid, n, &abort, userptr);
    return abort;
  }

  DYNALLOC1(int, myallp, myallp_sz, n*depth, "malloc");

  mygroupelts3(grp->levelinfo, n, depth - 1, action, NULL, myallp, myid, &abort, userptr);

  return abort;
}

void
myfreegroup(grouprec *grp)
/* Free (or pretend to free) group structure. */
{
  int i, j;
  cosetrec *p;
  permrec *q, *qq;

  for (i = 0; i < grp->depth; ++i)
  {
    p = grp->levelinfo[i].replist;
    if (p != NULL)
      for (j = grp->levelinfo[i].orbitsize; --j >= 0; )
      {
        myfreepermrec(p[j].rep, grp->n);
        p[j].rep = NULL;
      }
  }

  if (grp->depth > 0)
  {
    p = grp->levelinfo[0].replist;
    if (p != NULL && p != mycoset)
    {
      free(p);
      grp->levelinfo[0].replist = NULL;
    }

    q = grp->levelinfo[0].gens;
    while (q != NULL)
    {
      qq = q;
      q = q->ptr;
      myfreepermrec(qq, grp->n);
    }
    grp->levelinfo[0].gens = NULL;
  }
}



void
mygrouplevelproc(int *lab, int *ptn, int level, int *orbits, statsblk *stats,
  int tv, int index, int tcellsize, int numcells, int cc, int n)
{
  int depth;
  size_t sz = 0;

  if (numcells == n)   /* first call */
  {
    depth = level - 1;

    // if (mygroup) myfreegroup(mygroup);

    if (depth > mygroup_depth)
    {
      if (depth <= 1) sz = sizeof(grouprec);
      else            sz = sizeof(grouprec) + (depth - 1) * sizeof(levelrec);
      if (mygroup) mygroup = (grouprec*)realloc((void*)mygroup, sz);
      else         mygroup = (grouprec*)malloc(sz);
      if (mygroup == NULL)
      {
        fprintf(ERRFILE, ">E malloc failed in mygrouplevelproc\n");
        exit(1);
      }
      mygroup_depth = depth;
    }
    if (depth <= 1) sz = sizeof(grouprec);
    else            sz = sizeof(grouprec) + (depth - 1) * sizeof(levelrec);
    if (mygroup == NULL) mygroup = (grouprec*)malloc(sz);
    mygroup->n = n;
    mygroup->depth = depth;
    mygens = NULL;
    return;
  }

  if (mygroup == NULL) mygroup = (grouprec*)malloc(sz);
  mygroup->levelinfo[level - 1].fixedpt = tv;
  mygroup->levelinfo[level - 1].orbitsize = index;
  mygroup->levelinfo[level - 1].gens = mygens;
  mygroup->levelinfo[level - 1].replist = NULL;

  if (level == 1) mygroup->numorbits = stats->numorbits;
}




/************************************/
/** Automorphism group invocations **/
/************************************/
static void
myuserautomproc(int count, int *p, int *orbits,
  int numorbits, int stabvertex, int n)
{
  // copied from naugroup.c groupautomproc
  // copy the automorphism into the group
  permrec *perm;
  int o;

  perm = mynewpermrec(n);
  for (o = 0; o < n; ++o) perm->p[o] = p[o];
  perm->ptr = mygens;
  mygens = perm;
}


/**************************************************************************/

static void
recomputeGroup(graph *g, int n, booleann isClassEnumeration)
{
  int lab[MAXNN], ptn[MAXNN], orbits[MAXNN];
  statsblk stats;
  static DEFAULTOPTIONS_GRAPH(options);
  setword workspace[50];

  // copy the initial blank CRN coloring; put assignments in the same cell
  int currClassIndex = graphSize + speciesCount - 1;
  int currClass      = 0;
  for (int z = 0; z < n; z++) {
    if (z < graphSize) {
      lab[z] = colorLab[z];
      ptn[z] = colorPtn[z];
    }
    else if (z < graphSize + speciesCount) { // species partitioning (all in the same cell)
      lab[z] = z;
      ptn[z] = 1;
    }
    else { // class partitioning (one cell per class)
      lab[z] = z;

      if (z == currClassIndex + speciesClasses[currClass]) {
        ptn[z]         = 0;
        currClassIndex = currClassIndex + speciesClasses[currClass];
        currClass++;
      }
      else ptn[z] = 1;
    }
  }
  if (n >= graphSize + speciesCount) // terminate the species assignment nodes partition
    ptn[graphSize + speciesCount - 1] = 0;
  lab[n] = n;
  ptn[n] = 0; // cell terminator

  options.digraph       = TRUE;
  options.getcanon      = FALSE;
  options.defaultptn    = FALSE;
  options.userautomproc = myuserautomproc;
  options.userlevelproc = mygrouplevelproc;
  options.invarproc     = adjacencies;
  options.maxinvarlevel = n;

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50, 1, n, NULL);

  // set rigid flag (i.e. all orbits have size 1)
  states[n].isRigid = TRUE;
  int graphMax = (!isClassEnumeration) ? graphSize : graphSize + speciesCount;
  for (int i = 0; i < graphMax; i++)
    states[n].isRigid = states[n].isRigid && (orbits[i] == i);

  // set coset representatives
  mymakecosetreps(mygroup);

}

/**************************************************************************/
// from https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
int binomialCoeff(int n, int k)
{
  int r = 1;

  // Since C(n, k) = C(n, n-k)
  if (k > n - k)
    k = n - k;

  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    r *= (n - i);
    r /= (i + 1);
  }

  return r;
}

/* count of combinations of k-tuples from s species */
int complexCount(int s, int k) {
  return binomialCoeff(s + k - 1, k);
}

/* assign total number of complexes by kind, encoded in the color.
   color 0: naught, 1 in total
   color 1: homomers, s in total
   color 2: homodimers, s in total
   color 3: heterodimers, n-choose-2 in total */
int bimolecularComplexesCount(int s, int color) {
  switch (color) {
  case 0: return 1; break;
  case 1: return s; break;
  case 2: return s; break;
  case 3: return (s*(s - 1)) / 2; break;
  default: exit(-2); /* exiting program */
  }
}

/* the total number of complexes from 0 to m-molecular comprising s species.
   This is the sum: \Sigma_{i=0}^m \binom{n+i-1}{i} = \binom{s+m}{m}

   (see https://math.stackexchange.com/questions/870219/sum-of-k-combination-with-repetitions) */
int complexesTotal(int s, int k) {
  return binomialCoeff(s + k, k);
}

/**************************************************************************/
// reverse the bits in an integer
int reverse(int x)
{
  x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
  x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
  x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
  x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
  x = ((x >> 16) & 0xffffu) | ((x & 0xffffu) << 16);
  return x;
}

/* compute next bit permutation according to lexicographic order
  from: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        https://stackoverflow.com/questions/8281951/bit-hack-to-generate-all-integers-with-a-given-number-of-1s
  */
static int next_perm(unsigned int v) {
  unsigned int w; // next permutation of bits
  unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change,
  // set to 0 the least significant ones, and add the necessary 1 bits.
  unsigned int i;
#ifdef _WIN32
  _BitScanForward(&i, v);
  w = (t + 1) | (((~t & -~t) - 1) >> (i + 1));
#else
  w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
#endif
  return w;
}

booleann skipEdge(xword x, xword assigned, graph* g, int n) {
  // check that x is not coloring an already assigned monomer or monodimer
  if (x & assigned & colorMask[1]
    || x & assigned & colorMask[2]
    || x & assigned & colorMask[3])
    return 1;

  // check that x doesn't form duplicate heterodimers
  unsigned int newHeterodimerAssignment = reverse(x & colorMask[3]); // reverse due to g[graphSize+c] bits being in reverse order
  if (newHeterodimerAssignment) {
    unsigned int reversedColorMask = reverse(colorMask[3]);
    for (int c = 0; c < n - graphSize; c++) {
      unsigned int oldHeterodimerAssignemnt = g[graphSize + c] & reversedColorMask;
      if (!oldHeterodimerAssignemnt) continue;
      unsigned int newHeterodimers = newHeterodimerAssignment & oldHeterodimerAssignemnt;
      if (XPOPCOUNT(newHeterodimers) > 1)
        return 1;
    }
  }
  return 0;
}

/*********************************************************************/

// prints a CRN's graph encoding into LBS format 
static void printCRN(graph* g, int n, int spCount) {
  char** nodeNames = (char**)malloc(sizeof(char*) * graphSize);
  for (int z = 0; z < graphSize; z++) {
    int nodeIdx = MAXNN - z - 1;
    if (colorMask[0] & bit[nodeIdx] || colorMask[1] & bit[nodeIdx]) {
      nodeNames[z] = malloc(2);
    }
    else if (colorMask[2] & bit[nodeIdx]) {
      nodeNames[z] = malloc(5);
    }
    if (!(colorMask[3] & bit[nodeIdx])) {
      if (colorMask[0] & bit[nodeIdx]) {
        sprintf(nodeNames[z], "");
      }
      else {
        int idx = -1;
        for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
          xword speciesWord = g[graphSize + speciesIdx];
          if (bit[z] & speciesWord)
          {
            idx = speciesIdx; break;
          }
        }
        if (colorMask[2] & bit[nodeIdx]) sprintf(nodeNames[z], "2%c", speciesNames[idx]);
        else sprintf(nodeNames[z], "%c", speciesNames[idx]);
      }
    }
    else {
      int idx1 = -1, idx2 = -1;
      // find heterodimer's indices
      for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
        xword speciesWord = g[graphSize + speciesIdx];
        if (bit[z] & speciesWord)
          if (idx1 == -1)
            idx1 = speciesIdx;
          else { idx2 = speciesIdx; break; }
      }

      if (idx1 > idx2) {
        int tmp = idx1;
        idx1 = idx2;
        idx2 = tmp;
      }

      // create name and color
      nodeNames[z] = malloc(7);
      sprintf(nodeNames[z], "%c+%c", speciesNames[idx1], speciesNames[idx2]);
    }
  }

  for (int i = 0; i < graphSize; i++) {
    xword w = g[i];
    while (w) {
      int j = XNEXTBIT(w);
      w ^= xbit[j];
      int node = MAXNN - j - 1;
      if (node < graphSize)  // skip species nodes
        printf("|%s->%s\n", nodeNames[i], nodeNames[node]);
    }
  }

  printf("---------\n");

  for (int i = 0; i < graphSize; i++) {
    free(nodeNames[i]);
  }
  free(nodeNames);
}

static void printDotGraph(graph* g, int n, int spCount, int* colorMask) {

  printf("digraph {\n{\nnode [style=filled];\n");
  char* speciesNames = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  // 26 colors supported by GraphViz
  char* cls[26] = { "aliceblue", "crimson", "antiquewhite",  "aquamarine", "azure", "beige", "bisque", "blanchedalmond", "blue", "brown", "burlywood", "chartreuse", "coral", "cornflowerblue", "cornsilk", "cyan", "darkgoldenrod", "ghostwhite", "gold", "floralwhite", "forestgreen", "gainsboro", "wheat4", "white", "whitesmoke", "yellow"
  };

  char** nodeNames = (char**) malloc(sizeof(char*) * n);
  for (int z = 0; z < graphSize; z++) {
    int nodeIdx = MAXNN - z - 1;
    if (colorMask[0] & bit[nodeIdx] || colorMask[1] & bit[nodeIdx]) {
      nodeNames[z] = malloc(2);
    }
    else if (colorMask[2] & bit[nodeIdx]) {
      nodeNames[z] = malloc(5);
    }
    if (!(colorMask[3] & bit[nodeIdx])) {
      if (colorMask[0] & bit[nodeIdx]) nodeNames[z] = "0";
      else {
        int idx = -1;
        for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
          xword speciesWord = g[graphSize + speciesIdx];
          if (bit[z] & speciesWord)
          {idx = speciesIdx; break;}
        }
        if (colorMask[2] & bit[nodeIdx]) sprintf(nodeNames[z], "\"2%c\"", speciesNames[idx]);
        else sprintf(nodeNames[z], "%c", speciesNames[idx]);
      }
    }
    else {
      int idx1 = -1, idx2 = -1;
      // find heterodimer's indices
      for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
        xword speciesWord = g[graphSize + speciesIdx];
        if (bit[z] & speciesWord)
          if (idx1 == -1)
            idx1 = speciesIdx;
          else { idx2 = speciesIdx; break; }
      }

      if (idx1 > idx2) {
        int tmp = idx1;
        idx1 = idx2;
        idx2 = tmp;
      }

      // create name and color
      nodeNames[z] = malloc(7);
      sprintf(nodeNames[z], "\"%c  %c\"", speciesNames[idx1], speciesNames[idx2]);
    }
  }

  for (int z = 0; z < graphSize; z++) {
    {
      int nodeIdx = MAXNN - z - 1;
      int idx = MAXNN - z - 1;

      char* shape = "star";
      if      (colorMask[1] & bit[idx]) { shape = "square";    }
      else if (colorMask[2] & bit[idx]) { shape = "circle";    }
      else if (colorMask[3] & bit[idx]) { shape = "rectangle"; }

      if ( !(colorMask[3] & bit[idx]) ) {
        if (colorMask[0] & bit[idx])
          printf("0 [shape=point, width=0.2, fillcolor=black];\n");
        else printf("%s [shape=%s, fillcolor=%s];\n", nodeNames[z], shape, cls[z]);
      }
      else {
        int idx1 = -1, idx2 = -1;
        // find heterodimer's indices
        for (int colorIdx = 0; colorIdx < n - graphSize; colorIdx++) {
          xword colorWord = g[graphSize + colorIdx];
          if (bit[z] & colorWord)
            if (idx1 == -1)
              idx1 = colorIdx;
            else { idx2 = colorIdx; break; }
        }

        if (idx1 > idx2) {
          int tmp = idx1;
          idx1 = idx2;
          idx2 = tmp;
        }
        // create name and color
        printf("%s [shape=%s, style=striped, fillcolor=\"%s:%s\"];\n", nodeNames[z], shape, cls[idx1], cls[idx2]);
      }
    }
  }
  printf("}\n\n");

  booleann b = FALSE;
  int mask = (int) pow(2, graphSize) - 1;
  for (int i = 0; i < graphSize; i++) {
    xword w = g[i];
    if (w) {
      while (w) {
        int j = XNEXTBIT(w);
        w ^= xbit[j];
        int node = MAXNN - j - 1;
        if (node < graphSize)  // skip species nodes
        {
          if (!b) {
            printf("%s -> ", nodeNames[i]);
            b = TRUE;
            printf("%s", nodeNames[node]);
          }
          else {
            printf(", %s", nodeNames[node]);
          }
        }
      }

      if (b) printf("\n");
    }
    b = FALSE;
  }

  printf("---------------------\n");
}

/* Hashing function:
\text{Let } N \text{ be the total number of species, and } m \text{ be the maximum number of complexes.}
\\
\\ \textup{hash}(C_1 \rightarrow C_2) = \textup{hash}(C_1)*m+\textup{hash}(C_2)
\\
\\ \textup{hash}(C) = \begin{cases}
0                     & \text{if } C = \emptyset\\
1+i                     & \text{if } C = A_i\\
N+1+i                   & \text{if } C = 2A_i\\
2N+1+\textup{hash}(i,j) & \text{if } C = A_i+ A_j, i < j\\
\end{cases}
\\
\\ \textup{hash}(i, j) = \Big(N i - \frac{i (i+1)}{2}\Big) + (j - i -1)
*/
/* compute the hash number of a complex */
// CRN reactions for N <= 4 can be encoded in a single byte
byte hashHeterodimerByte(int i, int j) {
  byte offset = speciesCount*i - (i * (i + 1) / 2);
  byte index = j - i - 1; // -1 to make the hashing 0-indexed
  return offset + index;
}

byte hashComplexByte(graph* g, int complexIndex) {
  int assignmentsCount = 0;
  boolean firstSet = FALSE;
  byte A_i, A_j;
  for (byte i = graphSize; i < graphSize + speciesCount; i++) { // for each species assignment
    xword speciesWord = g[i];
    if (bit[complexIndex] & speciesWord) {
      assignmentsCount++;
      if (!firstSet) {
        A_i = i - graphSize;
        firstSet = TRUE;
      }
      else A_j = i - graphSize;
    }
  }

  int nodeIdx = MAXNN - complexIndex - 1;
  if (assignmentsCount == 0)           return 0;
  else if (assignmentsCount == 1){
    if (colorMask[1] & bit[nodeIdx])   return (1 + A_i);
    else                               return (speciesCount + 1 + A_i);
  }
  else    /* heterodimer */            return 2 * speciesCount + 1 + hashHeterodimerByte(A_i, A_j);
}

byte hashReactionByte(graph* g, int i, int j) {
  byte C_1 = hashComplexByte(g, i);
  byte C_2 = hashComplexByte(g, j);
  return ((byte) maxComplexes) * C_1 + C_2;
}

// CRN reactions for N > 4 can be encoded in two bytes
unsigned short int  hashHeterodimerShort(int i, int j) {
  unsigned short int  offset = speciesCount - (i * (i + 1) / 2);
  unsigned short int  index = j - i - 1;
  return offset + index;
}

unsigned short int hashComplexShort(graph* g, int complexIndex) {
  int assignmentsCount = 0;
  boolean firstSet = FALSE;
  unsigned short int A_i, A_j;
  int nodeIdx = MAXNN - complexIndex - 1;
  for (int i = graphSize; i < graphSize + speciesCount; i++) { // for each species assignment
    xword speciesWord = g[i];
    if (bit[complexIndex] & speciesWord) {
      assignmentsCount++;
      if (!firstSet) {
        A_i = i - graphSize;
        firstSet = TRUE;
      }
      else A_j = i - graphSize;
    }
  }

  if (assignmentsCount == 0)         return 0;
  else if (assignmentsCount == 1){
    if (colorMask[1] & bit[nodeIdx]) return (1 + A_i);
    else                             return (speciesCount + 1 + A_i);
  }
  else    /* heterodimer */          return 2 * speciesCount + 1 + hashHeterodimerShort(A_i, A_j);
}

unsigned short int hashReactionShort(graph* g, int i, int j) {
  unsigned short int C_1 = hashComplexShort(g, i);
  unsigned short int C_2 = hashComplexShort(g, j);
  return ((unsigned short int) maxComplexes) * C_1 + C_2;
}

static void printByteEncoding(graph* g){
  if (speciesCount <= 4){ // 1 byte encoding
    for (int i = 0; i < graphSize; i++) {
      xword w = g[i];
      while (w) {
        int j = XNEXTBIT(w);
        w ^= xbit[j];
        int node = MAXNN - j - 1;
        if (node < graphSize) { // skip species nodes
          byte x = hashReactionByte(g, i, node);
          printf("%c", x);
        }
      }
    }
  }
  else{ // 2 bytes encoding
    for (int i = 0; i < graphSize; i++) {
      xword w = g[i];
      while (w) {
        int j = XNEXTBIT(w);
        w ^= xbit[j];
        int node = MAXNN - j - 1;
        if (node < graphSize) { // skip species nodes
          short unsigned int x = hashReactionShort(g, i, node);
          char c1 = x & 0xff;
          char c2 = x >> 8 & 0xff ;
          printf("%c%c", c2, c1);
        }
      }
    }
  }
}

static void printERODE(graph* g, int n, int spCount, int ne) {
  char** nodeNames = (char**)malloc(sizeof(char*) * graphSize);
  for (int z = 0; z < graphSize; z++) {
    int nodeIdx = MAXNN - z - 1;
    if (colorMask[0] & bit[nodeIdx] || colorMask[1] & bit[nodeIdx]) {
      nodeNames[z] = malloc(2);
    }
    else if (colorMask[2] & bit[nodeIdx]) {
      nodeNames[z] = malloc(5);
    }
    if (!(colorMask[3] & bit[nodeIdx])) {
      if (colorMask[0] & bit[nodeIdx]) {
        sprintf(nodeNames[z], "");
      }
      else {
        int idx = -1;
        for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
          xword speciesWord = g[graphSize + speciesIdx];
          if (bit[z] & speciesWord)
          {
            idx = speciesIdx; break;
          }
        }
        if (colorMask[2] & bit[nodeIdx])
          sprintf(nodeNames[z], "%c+%c", speciesNames[idx], speciesNames[idx]);
        else sprintf(nodeNames[z], "%c", speciesNames[idx]);
      }
    }
    else {
      int idx1 = -1, idx2 = -1;
      // find heterodimer's indices
      for (int speciesIdx = 0; speciesIdx < spCount; speciesIdx++) {
        xword speciesWord = g[graphSize + speciesIdx];
        if (bit[z] & speciesWord)
          if (idx1 == -1)
            idx1 = speciesIdx;
          else { idx2 = speciesIdx; break; }
      }

      if (idx1 > idx2) {
        int tmp = idx1;
        idx1 = idx2;
        idx2 = tmp;
      }

      // create name and color
      nodeNames[z] = malloc(7);
      sprintf(nodeNames[z], "%c+%c", speciesNames[idx1], speciesNames[idx2]);
    }
  }

  for (int i = 0; i < ne; i++) {
    printf("rate K%i = 1.0;\r\n", i);
  }

  for (int i = 0; i < spCount; i++) {
    printf("init %c 1.0 |\r\n", speciesNames[i]);
  }
  for (int i = 0; i < ne; i++) {
    printf("init P%i 1.0 |\r\n", i);
  }

  int mask = (int) pow(2, graphSize) - 1;
  int r = 0;
  for (int i = 0; i < graphSize; i++) {
    xword w = g[i];
    while (w) {
      int j = XNEXTBIT(w);
      w ^= xbit[j];
      int node = MAXNN - j - 1;
      if (node < graphSize) { // skip species nodes
        char* reactant = nodeNames[i];
        char* product  = nodeNames[node];
        // if the reactant is naught
        if (reactant[0] == '\0')     printf("P%i->{K%i}%s+P%i|\n", r, r, product, r);
        // the product is naught
        else if (product[0] == '\0') printf("%s+P%i->{K%i}P%i|\n", reactant, r, r, r);
        // all other reactions
        else                         printf("%s+P%i->{K%i}%s+P%i|\n", reactant, r, r, product, r);
        r++;
      }
    }
  }

  // printf("---------\n");

  for (int i = 0; i < graphSize; i++) {
    free(nodeNames[i]);
  }
  free(nodeNames);
}
/**************************************************************************/
/*
   TODO: better explanation
   TODO: adjust explanation to new rigidStreak algorithm

   Check that no duplicate assignment are generated due to isomorphism.
   The "streak" variable counts how many species assignments have the same degree xc so far. Recall that assignments are added in increasing order of degree.
   Our definition of canonical form requires that assignements of the same degree are sorted in the following way:
    Let G be the graph where e1, ..., eN are valid species assignements with the same degree. G is canonical if:
       * in G1 = G\\[e1, ..., eN], e1 is the least canonical representative of the orbits in G' w.r.t e2, ..., eN.
       * in G2 = G\\[e2,..., eN], e2 is the least canonical representative of the orbits in G'' w.r.t e3, ..., eN.
       * and so on

   Therefore, when adding a new assignemtn to an intermediate graph G_i, we need to make sure that e_{i+1} is not actually a lesser orbit representative
   then the existing assignments e1, ..., e_{i-1}.

   Orbit information is stored in data, therefore we can quickly check if this is the case without calling Nauty again.
*/
struct canonTestArgs {
  graph* g;
  int currentSpeciesCount;
  int n;
  xword x;
};

int cmpfunc(const void *a, const void *b) {
  return (*(int*)a - *(int*)b);
}

/***** canonical form test *****/
static void myFun(int *p, int k, int *abort, void* gvoid) {
  /* Called for each automorphism in the group */

  if (first2)
  {
    /* only the identity */
    first2 = FALSE;
    return;
  }

  struct canonTestArgs* args = (struct canonTestArgs*) gvoid;
  graph* g = args->g;
  int n = args->n;
  int streak = args->currentSpeciesCount;
  if (streak <= 0) return; // TODO: abort?

  for (int z = 0; z < streak; z++) {
    // apply p node-wise on the xset
    pxsetStreak[z] = 0;
    xword w = g[n-streak+z];
    while (w)
    {
      int j1 = XNEXTBIT(w);
      w ^= xbit[j1];
      pxsetStreak[z] |= xbit[p[MAXNN - j1 - 1]];
    }
  }


  // sort the assignments, to ensure they are in canonical form
  qsort(pxsetStreak, streak, sizeof(xword), cmpfunc);

  for (int z = 0; z < streak; z++) {
    int i  = reverse(g[n-streak+z]);
    int pi = pxsetStreak[z];
    if (i < pi) {
      *abort = FALSE;
      return;
    }
    if (i > pi) {
      *abort = TRUE;
      return;
    }
  }
}

void extendGraph(graph* g, xword x, int n) {
  // create complete CRN in gx
  int i;
  g[n] = 0;

  int xw = x;
  while (xw)
  {
    i = XNEXTBIT(xw);
    xw ^= xbit[i];
    g[n] |= bit[i];
  }
}

// Depth-first search of connected components (stops as soon as it has been confirmed that all species are connected)
int connectedSpecies(graph* g, int current, int n, booleann* visited, int missingSpecies) {
  visited[current] = TRUE;
  int currMissing = missingSpecies;
  if (current >= graphSize) currMissing--;
  if (currMissing == 0)     return currMissing;

  // check targets from outgoing edges
  xword targets = g[current];
  while (targets) {
    int target = 0;
    target  = XNEXTBIT(targets);
    targets ^= xbit[target];
    target  = MAXNN - target - 1;
    if (visited[target] == FALSE) {
      currMissing = connectedSpecies(g, target, n, visited, currMissing);
      if (currMissing == 0) return currMissing;
    }
  }

  // check incoming edges
  unsigned int currentBit = bit[current];
  for (int i = 0; i < n; i++) {
    if ((g[i] & currentBit) != 0 && visited[i] == FALSE) {
      currMissing = connectedSpecies(g, i, n, visited, currMissing);
      if (currMissing == 0) return currMissing;
    }
  }

  return currMissing;
}

booleann isConnected(graph *g, int n) {
  // set visited nodes array
  booleann* visited = malloc(n * sizeof(booleann));
  for (int i = 0; i < n; i++)
    // set the naught node as visited, as it doesn't count for connectedness
    if ((colorMask[0] & xbit[i]) != 0) visited[i] = TRUE;
    else visited[i] = FALSE;

  // explored graph with DPS, returns number of unconnected species
  int missingSpecies = connectedSpecies(g, n-1, n, visited, speciesCount);
  free(visited);

  return (missingSpecies == 0);
}

// count edges
int countEdges(graph* g, int n) {
  int ne = 0;
  for (size_t i = 0; i < graphSize; i++) {
    ne += XPOPCOUNT(g[i]);
  }

  return ne;
}

double** makeIncidenceMatrix(graph* g, int n, int ne) {
  // init species-reactions matrix
  double** matrix = (double**)malloc(speciesCount * sizeof(double*));
  for (int i = 0; i < speciesCount; i++) {
    matrix[i] = malloc(ne * sizeof(double));
  }

  for (int i = 0; i < speciesCount; i++) {
    for (int j = 0; j < ne; j++) {
      matrix[i][j] = 0;
    }
  }

  // populate matrix
  int j = 0; // reaction index in the matrix
  for (int source = 0; source < graphSize; source++) {
    xword targets = g[source];
    while (targets) {
      int target = 0;
      target = XNEXTBIT(targets);
      targets ^= xbit[target];

      for (int i = 0; i < speciesCount; i++) {
        if ((g[graphSize + i] & bit[source]) != 0) {
          double k = 1;
          if      ((xbit[source] & colorMask[0]) != 0) k = 0;
          else if ((xbit[source] & colorMask[2]) != 0) k = 2;
          matrix[i][j] = matrix[i][j] - k;
        }

        if ((g[graphSize + i] & xbit[target]) != 0) {
          double k = 1;
          int targe2 = bit[target];
          if      ((targe2 & colorMask[0]) != 0) k = 0;
          else if ((targe2 & colorMask[2]) != 0) k = 2;
          matrix[i][j] = matrix[i][j] + k;
        }
      }

      j++;
    }
  }

  return matrix;
}

void freeIncidenceMatrix(double** matrix, int n) {
  for (int i = 0; i < n; i++) {
    free(matrix[i]);
  }
  free(matrix);
}

// check that all species can be produced and consumed by some species
booleann allProducedAndConsumed(int** matrix, int ne) {
  booleann allProduced = TRUE;
  booleann allConsumed = TRUE;
  for (int i = 0; i < speciesCount; i++) {
    int canProduce = FALSE;
    int canConsume = FALSE;
    for (int j = 0; j < ne; j++) {
      if (matrix[i][j] < 0) canConsume = TRUE;
      if (matrix[i][j] > 0) canProduce = TRUE;
    }

    if (!canProduce || !canConsume) {
      allProduced = canProduce;
      allConsumed = canConsume;
      break;
    }
  }

  return (allProduced && allConsumed);
}

void flipRows(double** matrix, int r1, int r2, int m) {
  double tmp;
  for (int i = 0; i < m; i++) {
    tmp = matrix[r1][i];
    matrix[r1][i] = matrix[r2][i];
    matrix[r2][i] = tmp;
  }
}

// Gauss-Jordan elimination, following the algorithm described here: http://people.math.carleton.ca/~kcheung/math/notes/MATH1107/wk04/04_gaussian_elimination.html
void makeReducedRowEchelon(double** matrix, int m, int n) {
  int p = 0;
  for (int k = 0; k < n; k++) {
    for (int i = p; i < m; i++) {
      // find non-zero value
      if (matrix[i][k] != 0) {
        if (i != p) {
          flipRows(matrix, i, p, n);
        }
        // Scale
        double a = matrix[p][k];
        for (int j = 0; j < n; j++) {
          matrix[p][j] = matrix[p][j] / a;
        }
        // Subtract
        for (int q = 0; q < m; q++) {
          if (q != p) {
            double alpha = matrix[q][k];
            for (int r = 0; r < n; r++) {
              matrix[q][r] = matrix[q][r] - matrix[p][r] * alpha;
            }
          }
        }

        p++;

        if (p > m) return;
      }
    }
  }
}

double gcd(double a, double b) {
  double temp;
  while (b != 0) {
    temp = fmod(a, b);

    a = b;
    b = temp;
  }
  return a;
}

double** makeFarkasArray(double** matrix, int m, int* matrixRowsDim, int* allocatedRowsCount) {
  int n = speciesCount;

  // extend the matrix by bufferSize rows at a time
  int bufferSize = 10;

  // extend the incidence matrix with an n x n identity matrix
  for (int i = 0; i < *matrixRowsDim; i++) {
    double* tmpPointer = realloc(matrix[i], (m + n) * sizeof(double));
    if (tmpPointer == NULL)
      fprintf(ERRFILE, ">E realloc failed in Farkas algorithm\n");
    else
      matrix[i] = tmpPointer;

    // set identity matrix
    for (int j = m; j < m + n; j++) {
      if (j - m == i)
        matrix[i][j] = 1.0;
      else matrix[i][j] = 0.0;
    }
  }

  int matrixBufferSize = n; // instead of allocating a single row at a time, we allocate a block of them at a time

  for (int i = 0; i < m; i++) {
    int r = *matrixRowsDim;

    for (int j = 0; j < r; j++) {
      for (int k = j + 1; k < r; k++) {
        double* d1 = matrix[j];
        double* d2 = matrix[k];

        if ((d1[i] > 0 && d2[i] < 0) || (d1[i] < 0L && d2[i] > 0L)) {
          // extend the matrix 
          if (*matrixRowsDim == matrixBufferSize) { // allocate extra bufferSize rows in the matrix if required
            matrixBufferSize += bufferSize;
            double** tmpPointer = (double**) realloc(matrix, (matrixBufferSize) * sizeof(double*));
            if (tmpPointer == NULL)
              fprintf(ERRFILE, ">E realloc failed in Farkas algorithm\n");
            else {
              matrix = tmpPointer;
              for (int z = 0; z < bufferSize; z++) {
                matrix[matrixBufferSize - bufferSize + z] = malloc((n + m) * sizeof(double));
              }
            }
          }

          // extend matrix size
          (*matrixRowsDim)++;

          double d1a = fabs(d1[i]);
          double d2a = fabs(d2[i]);
          double g = -1;
          for (int l = 0; l < m + n; l++) {
            double dd = d2a * d1[l] + d1a * d2[l];
            matrix[*matrixRowsDim - 1][l] = dd;
            if (g == -1) g = dd;
            else g = gcd(g, dd);
          }
          for (int l = 0; l < m + n; l++) {
            matrix[*matrixRowsDim - 1][l] /= g;
          }
        }
      }
    }

    // Delete all rows whose i-th component is different from 0.
    int idx = 0;
    int tmpMatrixRowsSize = *matrixRowsDim;
    for (int j = 0; j < *matrixRowsDim; j++) {
      if (matrix[j][i] == 0) { // keep this row
        if (idx != j) {
          // overwrite to idx
          for (int k = 0; k < n + m; k++) {
            matrix[idx][k] = matrix[j][k];
          }
        }

        idx++;
      }
      else { // discard this row
        tmpMatrixRowsSize--;
      }
    }

    *matrixRowsDim = tmpMatrixRowsSize;
  }

  *allocatedRowsCount = matrixBufferSize;
  return matrix;
}

// Quick Non-trivial test: Check if there is a positive vector in the stoichiometry matrix
booleann isNonTrivialQuick(double** matrix, int numSpecies, int numReactions) {
  booleann res = TRUE;
  for (int i = 0; i < numSpecies; i++) {
    booleann allNegative = TRUE;
    booleann allPositive = TRUE;
    booleann allZero = TRUE;
    for (int j = 0; j < numReactions; j++) {
      if (matrix[i][j] < 0) {
        allPositive = FALSE;
        allZero = FALSE;
      }
      else if (matrix[i][j] > 0) {
        allNegative = FALSE;
        allZero = FALSE;
      }
    }
    if ((allPositive == TRUE || allNegative == TRUE) && allZero != TRUE) {
      res = FALSE;
      break;
    }
  }

  return res;
}

void printMatrix(double** matrix, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      printf("%1.12f ", matrix[i][j]);
    printf("\n");
  }
  printf("\n");
}

/* Application of the Fourier-Motzkin Elimination algorithm to the case (\Gamma^{T}x >= 0) with x \neq 0, where \Gamma is the stoichimetry matrix of the CRN.
     (note: the stoichimoetry matrix \Gamma is not transposed in the code, we just reverse apply FME on it with indices reversed) */
booleann isNonTrivialFME(double** matrix, int numSpecies, int numReactions) {

  double tol = 1e-12;

  boolean solved = FALSE;
  int matrixSize = numReactions;

  int* pos = (int*)malloc(matrixSize * sizeof(int));
  int* neg = (int*)malloc(matrixSize * sizeof(int));
  boolean moreInequalitiesAdded = FALSE;

  for (int j = 0; j < numSpecies; j++) {
    int psize = 0; int nsize = 0; // size of P and N

    //printMatrix(matrix, numSpecies, matrixSize);

    // dynamically augment the set of positive and negative coefficient indices if more inequalities were created last time
    if (moreInequalitiesAdded == TRUE) {
      int* tmppos = (int*) realloc(pos, matrixSize * sizeof(int));
      int* tmpneg = (int*) realloc(neg, matrixSize * sizeof(int));

      if (tmppos == NULL || (tmpneg == NULL))
        fprintf(ERRFILE, ">E realloc failed in Fourier-Motzkin Elimination algorithm\n");
      else {
        pos = tmppos;
        neg = tmpneg;
      }

      moreInequalitiesAdded = FALSE;
    }

    boolean allCoefficientsAreZero = TRUE;
    boolean allCoefficientsArePositive = TRUE;
    boolean allCoefficientsAreNegative = TRUE;
    // optimization: the following two booleans are used to detect if there are inequalities of the form 
    // x_j >= 0 and -x_j >= 0    (i.e. single variable inequalities)
    boolean someMonoCoefficientIsPositive = FALSE;
    boolean someMonoCoefficientIsNegative = FALSE;
    for (int i = 0; i < matrixSize; i++) {
      if (fabs(matrix[j][i]) < tol) continue;
      //if (matrix[j][i] == 0.0) continue;
      else allCoefficientsAreZero = FALSE;

      // check if all other coefficients are zero in the row
      boolean onlyNonZeroCoefficient = TRUE;
      for (int k = j + 1; k < numSpecies; k++) {
        if (fabs(matrix[k][i]) > tol) { onlyNonZeroCoefficient = FALSE; break; }
        //if (matrix[k][i] != 0.0) { onlyNonZeroCoefficient = FALSE; break; }
      }

      // mark down all positive and negative coefficients
      if (onlyNonZeroCoefficient == TRUE) {
        if (matrix[j][i] > tol) {
          allCoefficientsAreNegative = FALSE;
          someMonoCoefficientIsPositive = TRUE;
          pos[psize] = i; psize++;
        }
        else {
          allCoefficientsArePositive = FALSE;
          someMonoCoefficientIsNegative = TRUE;
          neg[nsize] = i; nsize++;
        }
      }
      else {
        if (matrix[j][i] > tol) { pos[psize] = i; psize++; allCoefficientsAreNegative = FALSE; }
        else                    { neg[nsize] = i; nsize++; allCoefficientsArePositive = FALSE; }
      }
    }

    // optimization: if x_j >= 0 and x_j <= 0, x_j must be 0
    if (allCoefficientsAreZero == TRUE
      || (someMonoCoefficientIsNegative == TRUE
        && someMonoCoefficientIsPositive == TRUE)) continue; // assign 0 to this element of the vector

    // found a positive vector in x_j, by putting x_j to 1 or -1 and all other coefficients x_i with i \neq j to 0
    if ((allCoefficientsAreNegative || allCoefficientsArePositive)) {
      solved = TRUE; break;
    }

    if ((psize != 0 || nsize != 0)) {
      if (psize == 0) {       continue; }
      else if (nsize == 0) {  continue; }
      else {
        // Simplify x_j and create new inequalities, as per FME
        moreInequalitiesAdded = TRUE;

        // allocate more space for the new equations
        int newEquationsSize = matrixSize + (psize * nsize);
        for (int k = 0; k < numSpecies; k++) {
          double* tmpPointer = realloc(matrix[k], sizeof(double) * (newEquationsSize));
          if (tmpPointer == NULL)
            fprintf(ERRFILE, ">E realloc failed in Fourier-Motzkin Elimination algorithm\n");
          else
            matrix[k] = tmpPointer;
        }

        // create new inequalities
        for (int pi = 0; pi < psize; pi++) {
          for (int ni = 0; ni < nsize; ni++) {
            int newColIndex = matrixSize + (pi*nsize) + ni;
            for (int row = j + 1; row < numSpecies; row++) {
              double pelem = matrix[row][pos[pi]] / fabs(matrix[j][pos[pi]]);
              double nelem = matrix[row][neg[ni]] / fabs(matrix[j][neg[ni]]);
              double diff = pelem + nelem;
              matrix[row][newColIndex] = diff;
            }
          }
        }

        for (int k = j + 1; k < numSpecies; k++) {
          for (int pi = 0; pi < psize; pi++) matrix[k][pos[pi]] = 0;
          for (int ni = 0; ni < nsize; ni++) matrix[k][neg[ni]] = 0;
        }
        matrixSize = newEquationsSize;
      }
    }
  }

  free(pos);
  free(neg);

  return (solved == FALSE);
}

booleann isNonTrivial(graph* g, int ne) {
  double** matrix = makeIncidenceMatrix(g, speciesCount, ne);

  //makeReducedRowEchelon(matrix, speciesCount, ne); // disabled due to numerical errors in the procedure

  booleann res = isNonTrivialQuick(matrix, speciesCount, ne);
  if (res != FALSE) {
    res = isNonTrivialFME(matrix, speciesCount, ne);
  }
  freeIncidenceMatrix(matrix, speciesCount);

  return res;
}

booleann isCanon(graph* g, int n, booleann isClassEnumeration) {
  booleann isCanon;

  if (!states[n].isRigid) {
    first2 = TRUE;
    struct canonTestArgs args = { g, states[n].streak, n, isClassEnumeration };
    isCanon = !myallgroup3(mygroup, myFun, &args);
  }
  else isCanon = TRUE;

  return isCanon;
}


void saveState(int n) {
  saveGroupState(n);
  if (states[n].currentCardinality != 0)
  {
    resetGroupState();
  }
}

void resumeState(int n) {
  myfreegroup(mygroup);
  free(mygroup);
  if (mycoset) {
    if (mycoset->rep) myfreepermrec(mycoset->rep, -1);
    free(mycoset);
  }

  if (myallp) free(myallp);
  if (myid) free(myid);
  if (myfreelist) myfreepermrec(myfreelist, -1);
  loadGroupState(n);
}

void makeNextState(graph* g, int n, xword x, int xc, int i, booleann isClassEnumeration) {
  // add the new assignment
  extendGraph(g, x, n);

  // check if the group has to be recomputed
  booleann hasCardinalityIncreased = (xc != states[n].currentCardinality);

  // update which blank nodes are now fully or partially assigned
  if (!isClassEnumeration) {
    states[n + 1].assigned = states[n].assigned | (x & (colorMask[1] | colorMask[2] | states[n].assignedHeterodimers));
    states[n + 1].assignedHeterodimers = states[n].assignedHeterodimers | x & ~states[n].assignedHeterodimers & colorMask[3];
    states[n + 1].streak = !hasCardinalityIncreased ? states[n].streak + 1 : 1;
  }
  else {
    states[n + 1].assigned = states[n].assigned | x;
    states[n + 1].streak = 1;
  }
  states[n + 1].isRigid = states[n].isRigid;
  states[n + 1].currentCardinality = xc;
  states[n + 1].currentXSetIndex = i;

  if (hasCardinalityIncreased)
    recomputeGroup(g, n, isClassEnumeration);
}



void printAcceptedCRN(graph* g, int n, booleann isClassEnumeration) {
  graph* g1;

  if (!isClassEnumeration)
    g1 = g;
  else {
    g1 = malloc(n * sizeof(xword));

    // copy CRN reactions
    for (int i = 0; i < graphSize; i++) {
      g1[i] = g[i];
    }

    int currClass = 0;
    for (int i = graphSize; i < graphSize + speciesCount; ) {
      xword w = g[graphSize + speciesCount + currClass];

      while (w) {
        int j = XNEXTBIT(w);
        w ^= xbit[j];
        int node = MAXNN - j - 1;
        g1[i] = g[node];
        i++;
      }

      currClass++;
    }
  }

  if (printCrnSwitch)
    if (byteEncodingSwitch) printByteEncoding(g1);
    else if (!isClassEnumeration)
      printCRN(g1, n, speciesCount);
    else printCRN(g1, graphSize + speciesCount, speciesCount);
}

void accept2(graph* g, int n, booleann isClassEnumeration) {
  if (!isClassEnumeration) {
    booleann con = TRUE;
    if (connectedSwitch == TRUE) con = isConnected(g, n);

    if (con)
    {
      int ne = countEdges(g, n);
      booleann nonTrivial = TRUE;

      if (nonTrivialSwitch == TRUE) {
        nonTrivial = isNonTrivial(g, ne);
      }
      if (nonTrivial == TRUE)
      {
        booleann isConserving     = TRUE;
        booleann isNotConserving  = TRUE;
        booleann isMassConserving = TRUE;
        if (conservationLawSwitch == TRUE || nonConservationLawSwitch == TRUE || massConservingSwitch == TRUE) {
          double** matrix = makeIncidenceMatrix(g, n, ne);

          // compute Farkas array
          int matrixRowsDim = speciesCount;
          int mSize         = -1;
          matrix = makeFarkasArray(matrix, ne, &matrixRowsDim, &mSize);
          if (conservationLawSwitch == TRUE && matrixRowsDim == 0)    isConserving = FALSE;
          if (nonConservationLawSwitch == TRUE && matrixRowsDim != 0) isNotConserving = FALSE;
          if (massConservingSwitch == TRUE) {
            if (matrixRowsDim == 0)
              isMassConserving = FALSE;
            else {
              for (int i = 0; i < speciesCount; i++) {
                int speciesNotConserved = FALSE;
                for (int j = 0; j < matrixRowsDim; j++) {
                  if (matrix[j][ne + i] != 0.0) {
                    speciesNotConserved = TRUE;
                    break;
                  }
                }

                if (speciesNotConserved == FALSE) {
                  isMassConserving = FALSE;
                  break;
                }
              }
            }
          }
          freeIncidenceMatrix(matrix, mSize);
        }

        if ((!conservationLawSwitch && !nonConservationLawSwitch && !massConservingSwitch)
          || (conservationLawSwitch && isConserving)
          || (nonConservationLawSwitch && isNotConserving)
          || (massConservingSwitch && isMassConserving))
        {

          if (!speciesClassesSwitch) {
            printAcceptedCRN(g, n, isClassEnumeration);
            counter++;
          }
          else {
            saveState(n);
            // makeNextState(g, n, 0, 0, 0, TRUE);
            recomputeGroup(g, n, TRUE);
            states[n].assigned = 0;
            states[n].assignedHeterodimers = 0;
            states[n].streak = 0;
            states[n].currentCardinality = 1;
            states[n].currentXSetIndex = 1;

            crnextend(g, n, 0, TRUE);
            resumeState(n);
          }
        }
      }
    }
  }
  else {
    printAcceptedCRN(g, n, isClassEnumeration);
    counter++;
  }
}



void
crnextend(graph* g, int n, int ne, booleann isClassEnumeration)
/* adapted from genc.c;  */
{
  xword x;
  int nx, xc, i;

  nx = n + 1; // new graph size OR next species class

  // get last assignment's cardinality and xset index
  int currXc, currIdx;
  booleann isLastAssignment;
  if (!isClassEnumeration) {
    currXc           = states[n].currentCardinality;
    currIdx          = states[n].currentXSetIndex;
    isLastAssignment = (nx == maxn);
  }
  else {
    currXc           = speciesClasses[0];
    currIdx          = pow(2, currXc) - 1;
    isLastAssignment = (nx - graphSize - speciesCount == classCount);
  }

  if (isLastAssignment)
  {
    /* infer the last assignment.
       After choosing the previous n-1 species, there is only one possible assignment left to choose.
       Therefore we infer it from the graph, and check that it is valid and canonical */
    if (!isClassEnumeration) {
      x  = leastNBitsOn ^ states[n].assigned;
      i  = data.xinv[x];
      xc = data.xcard[i];

      // safety checks 
      if (xc < currXc || i < currIdx) return;            // assignments are only added in increasing order of cardinality and xset encoding.
      if (ne + xc != maxColorEdges) return;              // check that all complexes have species assigned.
      if (skipEdge(x, states[n].assigned, g, n)) return; // skip invalid species assignments.
    }
    else {
      x  = leastNBitsOnClass ^ states[n].assigned;
      i  = classData.xinv[x >> graphSize];
      xc = classData.xcard[i];
      
      // safety checks 
      if (xc + ne != classCardTotal) return;            // check that all species have been partitioned into classes
      if (x & states[n].assigned) return;               // check that each species is assigned to one class only
    }

    // recompute the group if the cardinality has increased
    booleann hasCardinalityIncreased = (xc != currXc);
    if (hasCardinalityIncreased && !isClassEnumeration)
      saveState(n);
    makeNextState(g, n, x, xc, i, isClassEnumeration);

    // accept the state if canonical 
    if (isCanon(g, nx, isClassEnumeration))
      accept2(g, nx, isClassEnumeration);

    // backtrack to previous state
    if (hasCardinalityIncreased && !isClassEnumeration)
      resumeState(n);
  }
  else {
    int max = !isClassEnumeration ? imax : classMax;
    for (i = currIdx; i < max; ++i) {
      // take the next valid assignment
      if (!isClassEnumeration) {
        x = data.xset[i];
        xc = data.xcard[i];

        // skip invalid assignments
        if (skipEdge(x, states[n].assigned, g, n)) continue;  
      }
      else {
        x  = classData.xset[i];    
        xc = classData.xcard[i]; 

        // safety checks
        if (xc != speciesClasses[n - graphSize - speciesCount]) continue; // exit if all assignments for the current species class have been tried out
        if (x & states[n].assigned) continue;                             // check that each species is assigned to one class only
      }

    // recompute the group if the cardinality has increased
      booleann hasCardinalityIncreased = (xc != currXc);
      if (hasCardinalityIncreased && !isClassEnumeration)
        saveState(n);
      makeNextState(g, n, x, xc, i, isClassEnumeration);

      // explore the new state if canonical
      if (isCanon(g, nx, isClassEnumeration)) {
        if (isClassEnumeration) {
          saveState(nx);
          recomputeGroup(g, nx, isClassEnumeration);
        }

        crnextend(g, nx, ne + xc, isClassEnumeration);

        if (isClassEnumeration) resumeState(nx);
      }

      // backtrack to previous state
      if (hasCardinalityIncreased && !isClassEnumeration) {
        resumeState(n);
      }
    }
  }
}

/// checks whether a species assignment is valid (e.g. the same species cannot occur in two separate homomers)
static booleann validAssignment(xword i, xword* colorMask, int maxHeterodimers) {
  if (i & colorMask[0]
    || XPOPCOUNT(i & colorMask[1]) > 1
    || XPOPCOUNT(i & colorMask[2]) > 1
    || XPOPCOUNT(i & colorMask[3]) > maxHeterodimers)
    return 0;
  else return 1;
}

static void
makeleveldata(booleann restricted, xword* colorMask, int speciesCount, booleann isClassEnumeration)
/* make the level data for each level */
{
  long h;
  int n, nn, j, currMaxDeg;
  long ncj;
  leveldata* d;
  xword* xcard;
  xword* xset, tttn, nxsets;

  if (!isClassEnumeration) {
    n = graphSize;
    nn = maxdeg <= n ? maxdeg : n;
    currMaxDeg = maxdeg;
  }
  else {
    // calculate the largest class' size
    int maxClassCardinality = -1;
    for (int i = 0; i < classCount; i++)
      if (speciesClasses[i] > maxClassCardinality) maxClassCardinality = speciesClasses[i];

    n = speciesCount;
    currMaxDeg = maxClassCardinality;
    nn = currMaxDeg <= n ? currMaxDeg : n; // TODO
  }

  ncj = nxsets = 1;
  /* compute the total number of xsets (i.e. all possible sets of edges, or the powerset over the set of edges in the graph up to n nodes) */
  for (j = 1; j <= nn; ++j)
  {
    ncj = (ncj * (n - j + 1)) / j;
    nxsets += ncj;
  }
  tttn = 1L << n;

  if (!isClassEnumeration) d = &data;
  else                     d = &classData;

  d->ne = d->dmax = d->xlb = d->xub = -1;

  xset = (xword*)calloc(nxsets, sizeof(xword));
  xcard = (xword*)calloc(nxsets, sizeof(xword));

  if (xset == NULL || xcard == NULL)
  {
    fprintf(stderr, ">E geng: calloc failed in makeleveldata()\n");
    exit(2);
  }

  j = 0;

  int maxHeterodimers = bimolecularComplexesCount(speciesCount, 3);
  for (int bits = 0; bits <= n; bits++) {
    int total = binomialCoeff(n, bits);
    int perm = (int)pow(2, bits) - 1;
    for (int z = 0; z < total; z++) {
      if ((h = XPOPCOUNT(perm)) <= currMaxDeg
        && (isClassEnumeration || validAssignment(perm, colorMask, maxHeterodimers)))
      {
        if (!isClassEnumeration) xset[j] = perm;
        else                     xset[j] = perm << graphSize; // shift xset to target species assignment nodes
        xcard[j] = h;
        ++j;
      }
      perm = next_perm(perm);
    }
  }
  d->xset = (xword*)calloc(j, sizeof(xword));
  d->xcard = (xword*)calloc(j, sizeof(xword));
  d->xinv = (xword*)calloc(tttn, sizeof(xword));
  for (int z = 0; z < tttn; z++) {
    d->xinv[z] = -1;
  }
  for (int z = 0; z < j; z++) {
    d->xset[z] = xset[z];
    d->xcard[z] = xcard[z];
    if (!isClassEnumeration)
      d->xinv[xset[z]] = z;
    else d->xinv[xset[z] >> graphSize] = z;
  }

  if (d->xset == NULL || d->xcard == NULL || d->xinv == NULL)
  {
    fprintf(stderr, ">E geng: calloc failed in makeleveldata()\n");
    exit(2);
  }

  if (!isClassEnumeration) imax = j;
  else                     classMax = j;

  free(xset);
  free(xcard);
}

/**************************************************************************/

static void
writeautom(int* p, int n)
/* Called by allgroup. */
{
  int i;

  for (i = 0; i < n; ++i) printf(" %2d", p[i]);
  printf("\n");
}


/**************************************************************************/
/************ adapted from vcolg.c from here on ***************************/
/**************************************************************************/

static int
ismax(int *p, int n) {
  /* test if col^p <= col */
  int i, k;
  int fail;

  fail = 0;
  for (i = 0; i < n; ++i)
  {
    k = p[i];
    if (k > fail)
      fail = k;
    if (col[k] > col[i])
    {
      fail_level = fail;
      return FALSE;
    }
    else
      if (col[k] < col[i])
        return TRUE;
  }

  ++newgroupsize;
  return TRUE;
}


static void
testmax(int *p, int n, int *abort)
/* Called by allgroup2. */
{
  int i;

  if (first)
  {                       /* only the identity */
    first = FALSE;
    return;
  }

  if (!ismax(p, n))
  {
    *abort = 1;
    for (i = 0; i < n; ++i)
      lastreject[i] = p[i];
    lastrejok = TRUE;
  }
}


static int
trythisone(grouprec *group, graph *g, booleann digraph, int m, int n
  , int reactionsCount
  , int* lab, int* ptn, int* orbits, optionblk options, statsblk stats, int* workspace)
  /* Try one solution, accept if maximal. */
  /* Return value is level to return to. */
{
  booleann accept;
  newgroupsize = 1;

  if (!group || groupsize == 1)
    accept = TRUE;
  else if (lastrejok && !ismax(lastreject, n))
    accept = FALSE;
  else if (lastrejok && groupsize == 2)
    accept = TRUE;
  else
  {
    newgroupsize = 1;
    first = TRUE;

    if (allgroup2(group, testmax) == 0)
      accept = TRUE;
    else
      accept = FALSE;
  }

  if (accept)
  {
    {
      blanksCounter++;
      // create a bit mask for color assignments. Color 0 is for naught, 1 for homomers, 2 for dimers, 3 for heterodimers
      for (int i = 0; i < 4; i++)
        colorMask[i] = 0;

      for (int idx = 0; idx < n; idx++)
        colorMask[col[idx]] |= 1 << idx;


      int assigned = 0; // bitmask that stores which complexes have been fully assigned to some species
      int assignedHeterodimers = 0; // bitmask to store partially filled heterodimers

      // naught cannot be assigned any species; set it as assigned
      for (int idx = 0; idx < n; idx++)
        if (col[idx] == 0) {
          assigned = 1 << idx;
          break;
        }

      graphSize = n;
      makeleveldata(FALSE, colorMask, speciesCount, FALSE);
      if (speciesClassesSwitch)
        makeleveldata(FALSE, colorMask, speciesCount, TRUE);

      // bit mask to get the complexes or species from a graph
      leastNBitsOn      = n == 32 ? 0xffffffff : (1 << n) - 1;
      leastNBitsOnClass = n == 32 ? 0xffffffff : (1 << speciesCount) - 1 << graphSize;

      // set the initials cell partition to reflect the complex coloring.
      // Also, set the total number of edges that must be drawn from the species nodes to the complexes node to have a full CRN
      int idx = 0;
      maxColorEdges = 0;
      for (int c = 0; c < 4; c++) { // TODO: remove hard-coded number of colors (4)
        int color = colorMask[c];
        while (color) {
          int j = XNEXTBIT(color);
          color ^= xbit[j];
          lab[idx] = j;
          if (color) ptn[idx] = 1;
          else ptn[idx] = 0;
          idx++;

          if (c == 1 || c == 2) ++maxColorEdges;
          else if (c == 3) maxColorEdges += 2;
        }
      }

      static DEFAULTOPTIONS_GRAPH(options2);

      options2.digraph       = TRUE;
      options2.defaultptn    = FALSE;
      options2.userautomproc = myuserautomproc;
      options2.userlevelproc = mygrouplevelproc;
      options2.invarproc     = adjacencies;
      options2.maxinvarlevel = n;

      set active[MAXMM];
      active[0] = 0;

      // call nauty to populate xorb
      nauty(g, lab, ptn, active, orbits, &options2, &stats, workspace, 50, 1, n, NULL);
      mymakecosetreps(mygroup);
      int rigid = stats.numorbits == n;

      mindeg = 1;
      colorLab = lab;
      colorPtn = ptn;


      states[n].assigned             = assigned;
      states[n].assignedHeterodimers = assignedHeterodimers;
      states[n].streak               = 0;
      states[n].isRigid              = rigid;
      states[n].currentCardinality   = 1;
      states[n].currentXSetIndex     = 1;

      crnextend(g, n, 0, FALSE);

      free(data.xset);
      free(data.xcard);
      free(data.xinv);
    }
    return n - 1;
  }
  else
    return fail_level - 1;
}

/**************************************************************************/

static int
scan(int level, graph *g, booleann digraph, int *prev, long minedges, long maxedges,
  long sofar, long numcols, grouprec *group, int m, int n,
  int reactionsCount, int *currMol, int *maxMol,
  int* lab, int* ptn, int* orbits, optionblk options, statsblk stats, int* workspace)
  /* Recursive scan for default case */
  /* Returned value is level to return to. */
{
  int left;
  long min, max, k, ret;

  if (level == n)
    return trythisone(group, g, digraph, m, n, reactionsCount, lab, ptn, orbits, options, stats, workspace);

  left = n - level - 1;
  min = minedges - sofar - numcols * left;
  if (min < 0)
    min = 0;
  if (noZeroNodeSwitch) min = 1;
  max = maxedges - sofar;
  if (max >= numcols)
    max = numcols - 1;
  if (prev[level] >= 0 && col[prev[level]] < max)
    max = col[prev[level]];

  for (k = min; k <= max; ++k)
  {
    col[level] = k;

    // keep track of the number of k-complexes 
    currMol[k]++;
    if (currMol[k] > maxMol[k])
    {
      // quit if there are more k-complexes than the max
      currMol[k]--;
      continue;
    }
    ret = scan(level + 1, g, digraph, prev, minedges, maxedges, sofar + k, numcols, group, m, n, reactionsCount, currMol, maxMol,
      lab, ptn, orbits, options, stats, workspace);
    currMol[k]--;
    if (ret < level)
      return ret;
  }

  return level - 1;
}



static void
colourdigraph(graph *g, int nfixed, long minedges, long maxedges,
  long numcols, int m, int n, int speciesCount, int molecularity)
{
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  setword workspace[MAXNV];
  grouprec *group;
  int i, j, k, nloops;
  size_t ii;
  set *gi, *gj, *gci, *gcj;
  int lab[MAXNV], ptn[MAXNV], orbits[MAXNV];
  booleann loop[MAXNV];
  int prev[MAXNV]; /* If >= 0, earlier point that must have greater colour */
  int weight[MAXNV];
  int region, start, stop;
  DYNALLSTAT(graph, gconv, gconv_sz);

  if (n > MAXNV) gt_abort(">E gencrn: MAXNV exceeded\n");
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  DYNALLOC2(graph, gconv, gconv_sz, n, m, "colourdigraph");

  nloops = 0;
  for (i = 0, gi = g; i < n; ++i, gi += m)
    if (ISELEMENT(gi, i))
    {
      DELELEMENT(gi, i);
      loop[i] = TRUE;
      ++nloops;
    }
    else
      loop[i] = FALSE;

  for (ii = 0; ii < m*(size_t)n; ++ii) gconv[ii] = g[ii];
  converse(gconv, m, n);

  for (region = 0; region < 2; ++region)
  {
    if (region == 0)
    {
      if (nfixed == 0)
        continue;
      start = 0;
      stop = nfixed;
      if (stop > n)
        stop = n;
    }
    else
    {
      if (nfixed >= n)
        continue;
      start = nfixed;
      stop = n;
    }

    for (i = start, gi = g + m * (size_t)start, gci = gconv + m * (size_t)start
      ; i < stop
      ; ++i, gi += m, gci += m)
    {
      /* Find most recent equivalent j. */
      for (j = i - 1, gj = gi - m, gcj = gci - m
        ; j >= start
        ; --j, gj -= m, gcj -= m)
      {
        if (loop[j] != loop[i]
          || ISELEMENT(gi, j) != ISELEMENT(gj, i)) continue;
        for (k = 0; k < m; ++k)
          if (gi[k] != gj[k] || gci[k] != gcj[k])
            break;
        if (k < m)
        {
          FLIPELEMENT(gi, i); FLIPELEMENT(gj, j);
          FLIPELEMENT(gci, i); FLIPELEMENT(gcj, j);
          for (k = 0; k < m; ++k)
            if (gi[k] != gj[k] || gci[k] != gcj[k])
              break;
          FLIPELEMENT(gci, i); FLIPELEMENT(gcj, j);
          FLIPELEMENT(gi, i); FLIPELEMENT(gj, j);
        }
        if (k == m)
          break;
      }

      if (j >= start)
      {
        prev[i] = j;
        weight[i] = weight[j] + 1;
      }
      else
      {
        prev[i] = -1;
        weight[i] = 0;
      }
    }
  }

  for (i = nfixed; i < n; ++i)
    weight[i] += nfixed;

  if (maxedges == NOLIMIT || maxedges > n*numcols)
    maxedges = n * numcols;
  if (n*numcols < minedges)
    return;
  options.userautomproc = groupautomproc;
  options.userlevelproc = grouplevelproc;
  options.defaultptn = FALSE;
  options.digraph = TRUE;
  options.invarproc = adjacencies;
  options.maxinvarlevel = n;

  setlabptn(weight, lab, ptn, n);

  if (nloops > 0)
    for (i = 0, gi = g; i < n; ++i, gi += m)
      if (loop[i]) ADDELEMENT(gi, i);

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, MAXNV, m, n, NULL);

  if (stats.grpsize2 == 0)
    groupsize = stats.grpsize1 + 0.1;
  else
    groupsize = 0;

  group = groupptr(FALSE);
  makecosetreps(group);

  if (stats.numorbits < n)
  {
    j = n;
    for (i = 0; i < n; ++i)
      if (orbits[i] < i && orbits[i] < j) j = orbits[i];

    for (i = j + 1; i < n; ++i)
      if (orbits[i] == j) prev[i] = j;
  }

  lastrejok = FALSE;
  for (i = 0; i < n; ++i)
    col[i] = 0;

  int molStart = 0;
  if (noZeroNodeSwitch) {
    currMol[0] = 0;
    maxMol[0] = 0;
    molStart = 1;
  }
  for (int i = molStart; i < molecularity; i++)
  {
    currMol[i] = 0;
    maxMol[i] = bimolecularComplexesCount(speciesCount, i);
  }

  scan(0, g, TRUE, prev, minedges, maxedges, 0, molecularity, group, m, n,
    n, currMol, maxMol, lab, ptn, orbits, options, stats, workspace);
  freegroup(group);
}

/**************************************************************************/
int
mainMethod(int argc, char *argv[])
{
#ifdef _WIN32
  LARGE_INTEGER frequency;
  LARGE_INTEGER start;
  LARGE_INTEGER end;
#else
  double frequency;
  double start;
  double end;
  double t;
#endif
  double interval;
  graph *g;
  int m, n, codetype, /*reactionsCount,*/ molecularity = 2;
  int argnum, j, nfixed;
  char *arg, sw;
  booleann badargs, digraph;
  booleann fswitch, uswitch, eswitch, qswitch, mswitch;
  booleann hasSpeciesCount;
  long minedges, maxedges, numcols;
  char *infilename, *outfilename;
  FILE *infile;
  char msg[201];
  size_t msglen;

  HELP; PUTVERSION;

  nauty_check(WORDSIZE, 1, 1, NAUTYVERSIONID);

  fswitch = Tswitch = FALSE;
  uswitch = eswitch = mswitch = qswitch = noZeroNodeSwitch = FALSE;
  infilename = outfilename = NULL;

  argnum = 0;
  badargs = FALSE;
  for (j = 1; !badargs && j < argc; ++j)
  {
    arg = argv[j];
    if (arg[0] == '-' && arg[1] != '\0')
    {
      ++arg;
      while (*arg != '\0')
      {
        sw = *arg++;
        SWBOOLEAN('q', qswitch)
    else SWBOOLEAN('u', uswitch)
    else SWBOOLEAN('T', Tswitch)
    else SWBOOLEAN('c', connectedSwitch)          // connected
    else SWBOOLEAN('l', conservationLawSwitch)    // conservation law
    else SWBOOLEAN('x', nonConservationLawSwitch) // non-conservation law
    else SWBOOLEAN('z', noZeroNodeSwitch)         // do not include zero
    else SWBOOLEAN('t', nonTrivialSwitch)         // filter out trivial CRNs
    else SWBOOLEAN('m', massConservingSwitch)     // filter out non-mass conserving CRNs
    else SWBOOLEAN('b', byteEncodingSwitch)       // print CRNs using byte encoding
    else SWBOOLEAN('T', Tswitch)                  // non-conservation law
    else SWINT('n', hasSpeciesCount, speciesCount, "gencrn -n")
    else SWSEQUENCE('s', ';', speciesClassesSwitch, speciesClasses2, 100, speciesClasses, "gencrn -s")
    // else SWINT('m', hasMolecularity, molecularity, "gencrn -m")
    //else SWINT('f',fswitch,nfixed,"vcolg -f")
    //else SWRANGE('e',":-",eswitch,minedges,maxedges,"vcolg -e")
    else badargs = TRUE;
      }
    }
    else
    {
      ++argnum;
      if (argnum == 1) infilename = arg;
      else if (argnum == 2) outfilename = arg;
      else                  badargs = TRUE;
    }
  }

  printCrnSwitch = !qswitch;

  if (badargs || argnum > 2)
  {
    fprintf(stderr, ">E Usage: %s\n", USAGE);
    GETHELP;
    exit(1);
  }

  if (!mswitch) numcols = 2;
  if (!eswitch)
  {
    minedges = 0;
    maxedges = NOLIMIT;
  }
  if (!fswitch) nfixed = 0;

  if (!qswitch)
  {
    msg[0] = '\0';
    CATMSG0(">A gencrn");
    if (eswitch || mswitch || uswitch || (fswitch && nfixed > 0)
      || Tswitch)
      CATMSG0(" -");
    if (mswitch)                  CATMSG1("m%ld", numcols);
    if (uswitch)                  CATMSG0("u");
    if (Tswitch)                  CATMSG0("T");
    if (noZeroNodeSwitch)         CATMSG0("z");
    if (connectedSwitch)          CATMSG0("c");
    if (conservationLawSwitch)    CATMSG0("l");
    if (nonTrivialSwitch)         CATMSG0("t");
    if (nonConservationLawSwitch) CATMSG0("x");
    if (massConservingSwitch)     CATMSG0("m");
    if (byteEncodingSwitch)       CATMSG0("b");
    if (speciesClassesSwitch)     CATMSG0("s");
    if (fswitch) CATMSG1("f%d", nfixed);
    if (eswitch) CATMSG2("e%ld:%ld", minedges, maxedges);
    msglen = strlen(msg);
    if (argnum > 0) msglen += strlen(infilename);
    if (argnum > 1) msglen += strlen(outfilename);
    if (msglen >= 196)
    {
      fputs(msg, stderr);
      if (argnum > 0) fprintf(stderr, " %s", infilename);
      if (argnum > 1) fprintf(stderr, " %s", outfilename);
      fprintf(stderr, "\n");
    }
    else
    {
      if (argnum > 0) CATMSG1(" %s", infilename);
      if (argnum > 1) CATMSG1(" %s", outfilename);
      CATMSG0("\n");
      fputs(msg, stderr);
    }
    fflush(stderr);
  }

  if (infilename && infilename[0] == '-') infilename = NULL;
  infile = opengraphfile(infilename, &codetype, FALSE, 1);
  if (!infile) exit(1);
  if (!infilename) infilename = "stdin";

  if (uswitch)
    outfile = NULL;
  else
  {
    if (!outfilename || outfilename[0] == '-')
    {
      outfilename = "stdout";
      outfile = stdout;
    }
    else if ((outfile = fopen(outfilename, "w")) == NULL)
    {
      fprintf(stderr, "Can't open output file %s\n", outfilename);
      gt_abort(NULL);
    }
  }

  vc_nin = vc_nout = 0;

#ifndef _WIN32
  t = CPUTIME;
#else
  QueryPerformanceFrequency(&frequency);
  QueryPerformanceCounter(&start);
#endif

  // set species classes data if specified by the user
  if (speciesClassesSwitch) {
    classCount = speciesClasses[0];
    for (int i = 0; i < classCount; i++) {
      speciesClasses[i] = (int)speciesClasses2[i];
      classCardTotal += speciesClasses[i];
    }
		if (classCardTotal != speciesCount)
		{
			fprintf(stderr, "Species colouring sums to %d, which differs from the species count (-n%d). Exiting...\n", classCardTotal, speciesCount);
			gt_abort(NULL);
		}
  }

  boolean byteEncodingPrinted = FALSE;
  while (TRUE)
  {
    if ((g = readggcrn(infile, NULL, 0, &m, &n, &digraph, speciesCount, classCount)) == NULL) break;
    ++vc_nin;

    numcols = molecularity = 4; // molecularity;
    maxComplexes = complexesTotal(speciesCount, 2);
    if (byteEncodingSwitch && byteEncodingPrinted == FALSE) {
      int reactionsCount = 0;
      for (int i = 0; i < n; i++)
        reactionsCount += XPOPCOUNT(g[i]);
      byteEncodingPrinted = TRUE;
#ifdef _WIN32
      // set stdout in binary mode to escape ASCII characters
      int setmode_r = _setmode(_fileno(stdout), _O_BINARY); // avoids e.g. that \n automatically printed as \r\n and disrupt the byte output
      if (setmode_r == -1) {
        fprintf(stderr, "Error: cannot set stdout to binary mode.\n");
        exit(-1);
      }
      printf("%i\r\n%i\r\n", speciesCount, reactionsCount);
#else
      printf("%i\n%i\n", speciesCount, reactionsCount);
#endif

    }
    if (noZeroNodeSwitch == TRUE) maxComplexes--;
    if (n <= maxComplexes) {
      maxn = n + (speciesCount);
      mindeg = 0;
      // a species can occur at most in a homomer, a homodimer and (n-1) heterodimers in a CRN
      maxdeg = 2 + (speciesCount - 1);
      colourdigraph(g, nfixed, minedges, maxedges, numcols, m, n, speciesCount, molecularity);
    }

    if (!uswitch && ferror(outfile)) gt_abort(">E gencrn output error\n");
    FREES(g);
  }
#ifndef _WIN32
  t = CPUTIME - t;
  interval = t;
#else
  QueryPerformanceCounter(&end);
  interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
#endif

#ifdef SUMMARY
  SUMMARY();
#endif
  if(!byteEncodingPrinted){
    printf("CRNs generated: %llu\n", counter);
    printf("Time elapsed: %f seconds\n", interval);
  }

  exit(0);
}

// Non-trivial dynamics tests
void testTrivial1() {
  // A->B+C | ->B+C
  int numSpecies = 3;
  int numReactions = 2;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = -1.0;
  matrix[0][1] = 0.0;
  matrix[1][0] = 1.0;
  matrix[1][1] = 1.0;
  matrix[2][0] = 1.0;
  matrix[2][1] = 1.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == TRUE)
    printf("Test error: Trivial1 returned non-trivial in FME\n");
}

void testTrivial2() {
  // A + B -> | 2B -> A + C | 2A -> B + C | C -> B
  int numSpecies = 3;
  int numReactions = 4;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = -1.0;
  matrix[0][1] = 1.0;
  matrix[0][2] = -2.0;
  matrix[0][3] = 0.0;
  matrix[1][0] = -1.0;
  matrix[1][1] = -2.0;
  matrix[1][2] = 1.0;
  matrix[1][3] = 1.0;
  matrix[2][0] = 0.0;
  matrix[2][1] = 1.0;
  matrix[2][2] = 1.0;
  matrix[2][3] = -1.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == TRUE)
    printf("Test error: Trivial2 returned non-trivial in FME\n");
}

void testNontrivial1() {
  // A+B -> A+C | 2C -> B+C
  int numSpecies = 3;
  int numReactions = 2;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = 0.0;
  matrix[0][1] = 0.0;
  matrix[1][0] = -1.0;
  matrix[1][1] = 1.0;
  matrix[2][0] = 1.0;
  matrix[2][1] = -1.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == FALSE)
    printf("Test error: NonTrivial1 returned trivial in FME\n");
}

void testNontrivial2() {
  // A -> B+C | -> B+C | B+C -> A | B+C -> 
  int numSpecies = 3;
  int numReactions = 4;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = -1.0;
  matrix[0][1] = 0.0;
  matrix[0][2] = 1.0;
  matrix[0][3] = 0.0;
  matrix[1][0] = 1.0;
  matrix[1][1] = 1.0;
  matrix[1][2] = -1.0;
  matrix[1][3] = -1.0;
  matrix[2][0] = 1.0;
  matrix[2][1] = 1.0;
  matrix[2][2] = -1.0;
  matrix[2][3] = -1.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == FALSE)
    printf("Test error: NonTrivial2 returned trivial in FME\n");
}

void testNontrivial3() {
  // 2A->B+C | A+C->B+D | B+F->D+F | D+E->A+E
  int numSpecies = 6;
  int numReactions = 4;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = -2.0;
  matrix[0][1] = -1.0;
  matrix[0][2] = 0.0;
  matrix[0][3] = 1.0;

  matrix[1][0] = 1.0;
  matrix[1][1] = 1.0;
  matrix[1][2] = -1.0;
  matrix[1][3] = 0.0;

  matrix[2][0] = 1.0;
  matrix[2][1] = -1.0;
  matrix[2][2] = 0.0;
  matrix[2][3] = 0.0;

  matrix[3][0] = 0.0;
  matrix[3][1] = 1.0;
  matrix[3][2] = 1.0;
  matrix[3][3] = -1.0;

  matrix[4][0] = 0.0;
  matrix[4][1] = 0.0;
  matrix[4][2] = 0.0;
  matrix[4][3] = 0.0;

  matrix[5][0] = 0.0;
  matrix[5][1] = 0.0;
  matrix[5][2] = 0.0;
  matrix[5][3] = 0.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == FALSE)
    printf("Test error: NonTrivial3 returned trivial in FME\n");
}

void testNontrivial3_Hand() {
  // D+E->A+E | B+F->D+F | A+C->B+D | 2A->B+C
  int numSpecies = 6;
  int numReactions = 4;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  // E
  matrix[0][0] = 0.0;
  matrix[0][1] = 0.0;
  matrix[0][2] = 0.0;
  matrix[0][3] = 0.0;

  // F
  matrix[1][0] = 0.0;
  matrix[1][1] = 0.0;
  matrix[1][2] = 0.0;
  matrix[1][3] = 0.0;

  // C
  matrix[2][0] = 0.0;
  matrix[2][1] = 0.0;
  matrix[2][2] = -1.0;
  matrix[2][3] = 1.0;

  // A
  matrix[3][0] = 1.0;
  matrix[3][1] = 0.0;
  matrix[3][2] = -1.0;
  matrix[3][3] = -2.0;

  // D
  matrix[4][0] = -1.0;
  matrix[4][1] = 1.0;
  matrix[4][2] = 1.0;
  matrix[4][3] = 0.0;

  // B
  matrix[5][0] = 0.0;
  matrix[5][1] = -1.0;
  matrix[5][2] = 1.0;
  matrix[5][3] = 1.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == FALSE)
    printf("Test error: NonTrivial3_Hand returned trivial in FME\n");
}

void testNontrivial4() {
  // C+E->2A | A+E->B+C | B+F->E+F | C+D->B+D
  int numSpecies = 6;
  int numReactions = 4;
  double** matrix = (double**)malloc(numSpecies * sizeof(double*));
  for (int i = 0; i < numSpecies; i++) {
    matrix[i] = malloc(numReactions * sizeof(double));
  }

  matrix[0][0] = -1.0;
  matrix[0][1] = 1.0;
  matrix[0][2] = 0.0;
  matrix[0][3] = -1.0;

  matrix[1][0] = -1.0;
  matrix[1][1] = -1.0;
  matrix[1][2] = 1.0;
  matrix[1][3] = 0.0;

  matrix[2][0] = 2.0;
  matrix[2][1] = -1.0;
  matrix[2][2] = 0.0;
  matrix[2][3] = 0.0;

  matrix[3][0] = 0.0;
  matrix[3][1] = 1.0;
  matrix[3][2] = -1.0;
  matrix[3][3] = 1.0;

  matrix[4][0] = 0.0;
  matrix[4][1] = 0.0;
  matrix[4][2] = 0.0;
  matrix[4][3] = 0.0;

  matrix[5][0] = 0.0;
  matrix[5][1] = 0.0;
  matrix[5][2] = 0.0;
  matrix[5][3] = 0.0;

  if (isNonTrivialFME(matrix, numSpecies, numReactions) == FALSE)
    printf("Test error: NonTrivial4 returned trivial in FME\n");
}


/**************************************************************************/
int main(int argc, char *argv[]) {

  /*testTrivial1();
  testTrivial2();
  testNontrivial1();
  testNontrivial2();
  testNontrivial3();
  testNontrivial4();
  testNontrivial3_Hand();*/

  return mainMethod(argc, argv);
}