// splitCRN.cpp : Split an input byte-encoded CRN file into N chunks.

#include <iostream>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

int numberOfDigits(int n) {
  if (n == 0) return 1; else return floor(log10(abs(n))) + 1;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: %s [number of chunks] [input file name]", argv[0]);
    exit(EXIT_FAILURE);
  }

  int s, rmax;
  char* arg, sw;

  // arguments parsing and validation
  char* end = NULL;
  errno = 0;
  
  int chunks;
  long temp = strtol(argv[argc - 2], &end, 10);
  if (end != argv[argc - 2] && errno != ERANGE && temp >= INT_MIN && temp <= INT_MAX)
  {
    chunks = (int)temp;
  }
  else {
    printf("Usage: %s [number of chunks] [input file name]", argv[0]);
    exit(EXIT_FAILURE);
  }


  // generation of geng and directg commands
  FILE* input;
  int asd = fopen_s(&input, argv[2], "rb");

  if (input == NULL) {
    printf("Cannot read input file");
    exit(EXIT_FAILURE);
  }
  int speciesCount = 0;
  int reactionsCount = 0;
  long int czasdur = ftell(input);
  int res = fscanf_s(input, "%d\n", &speciesCount);
  int rekls = fscanf_s(input, "%d\n", &reactionsCount);

  int reactionSize = speciesCount < 4 ? 1 : 2;

  // calculate number of CRNs
  long int cur = ftell(input);
  fseek(input, 0, SEEK_END);
  long int endp = ftell(input);
  fseek(input, cur-1, SEEK_SET);

  // from https://stackoverflow.com/questions/2422712/rounding-integer-division-instead-of-truncating
  long int totalCRNs = reactionSize == 2 ? (endp - cur + 1) / (2 * reactionsCount) : (endp - cur + 1) / reactionsCount;
  long int remainder = totalCRNs % chunks;
  long int crnsPerFile = (totalCRNs - remainder) / chunks;

  int fileLen = strlen(argv[2]);
  
  int bufferLen = sizeof(char) * reactionSize * reactionsCount;
  char* buffer = (char*)malloc(bufferLen);
  
  for (int i = 0; i < chunks; i++) {
    // create file name
    int newPathLen = (sizeof(char) * fileLen + 1L) + (sizeof(int) * numberOfDigits(i));
    char* path = (char*)malloc(newPathLen);
    sprintf_s(path, newPathLen, "%s_%d", argv[2], i);
    
    // create file
    FILE* output;
    fopen_s(&output, path, "w");
    if (output == NULL) { printf("Cannot create output file %s", path);  exit(-1);} //error

    // write number of species and reactions
    fprintf(output, "%d\n%d\n", speciesCount, reactionsCount);
    fclose(output);
    
    // reopen the file in binary mode, to escape \n
    fopen_s(&output, path, "ab");
    if (output == NULL) { printf("Cannot create output file %s", path);  exit(-1); } //error
    fseek(output, 0, SEEK_END);

    // copy CRNs
    int crns;
    if (remainder > 0) {
      crns = crnsPerFile + 1;
      remainder--;
    }
    else {
      crns = crnsPerFile;
    }
    for (int j = 0; j < crns; j++){
      int read = fread(buffer, 1, bufferLen, input);
      if (read == EOF || (i == chunks -1 && read == 0)) break;
      if (read < bufferLen) { printf("Internal error: read misaligned CRN");  exit(-1); } //error
      int written = fwrite(buffer, 1, bufferLen, output);
      int asd = 0;
    }
    
    fclose(output);
    
    free(path);
  }

  free(buffer);
  fclose(input);

  exit(EXIT_SUCCESS);
}