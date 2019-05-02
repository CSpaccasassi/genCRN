#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#ifndef _WIN32
#include <limits.h>
#endif
// only bimolecular reactions for now
int totalNumberOfNodes(int numberOfSpecies){
  return 1                                                  // naught node
         + numberOfSpecies                                  // homomers
         + numberOfSpecies                                  // homodimers
         + ((numberOfSpecies*(numberOfSpecies - 1)) / 2);   // heterodimers
}


int minNumberOfNodes(int s, int r){
  int n = 0;
  for (int i = 0; i < r; i += n) {
    n++;
  }

  return n;
}

int main(int argc, char *argv[]){
  if (argc < 3){
    printf("Usage: %s [number of species] [number of reactions]", argv[0]);
    exit(EXIT_FAILURE);
  }

  int s, rmax;
  char *arg, sw;
  int rswitch = 0;

  // Flags
  for (int j = 1; j < argc-2; ++j)
  {
	  arg = argv[j];
	  if (arg[0] == '-' && arg[1] != '\0')
	  {
		  ++arg;
		  while (*arg != '\0')
		  {
			  sw = *arg++;
			  if (sw == 'r')
				  rswitch = 1;
		  }
	  }
  }

  // arguments parsing and validation
  char *end = NULL;
  errno = 0;
  long temp = strtol(argv[argc - 2], &end, 10);
  if (end != argv[argc - 2] && errno != ERANGE && temp >= INT_MIN && temp <= INT_MAX)
  {
    s = (int)temp;
  }
  else {
    printf("Usage: %s [number of species] [number of reactions]", argv[0]);
    exit(EXIT_FAILURE);
  }

  temp = strtol(argv[argc - 1], &end, 10);
  if (end != argv[argc - 1] && errno != ERANGE && temp >= INT_MIN && temp <= INT_MAX)
  {
    rmax = (int)temp;
  }
  else {
    printf("Usage: %s [number of species] [number of reactions]", argv[0]);
    exit(EXIT_FAILURE);
  }

  // generation of geng and directg commands
  int nmin = minNumberOfNodes(s,rmax);
  int rmin = rmax%2 ? rmax/2 : (rmax+1)/2;
  int tot = totalNumberOfNodes(s);
  int nmax = 0;
  if (2*rmax < tot) // in case there are not enough reactions to connect all nodes
         nmax = 2*rmax;
    else nmax = tot;

  if (nmin > nmax) nmin = nmax;

  for (int n = nmin; n <= nmax; n++) {
    for (int r = rmin; r <= rmax; r++){
	  if (rswitch==1)
		printf("./geng.exe %i %i -d1 >> rcrn_%i_%i.txt\n", n, r, s, rmax);
	  else
		printf("./geng.exe %i %i -d1 | ./directg.exe -e%i >> crn_%i_%i.txt\n", n, r, rmax, s, rmax);
    }
  }

  exit(EXIT_SUCCESS);
}
