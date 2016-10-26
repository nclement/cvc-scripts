#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#define BUFF_SIZE 1024

char* getCol(char* line, int colNo) {
  char* tok = strtok(line, " ");
  for (int i = 0; i < colNo; ++i) {
    if (!tok) {
      return NULL;
    }
    tok = strtok(NULL, " ");
  }
  return tok;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    fprintf(stderr, "usage: summaryStatCol <col#> <infile> [outlier SD]\n");
    return -1;
  }

  int col = atoi(argv[1]);
  int outlierSD = -1;
  if (argc == 4) {
    outlierSD = atoi(argv[3]);
    //printf("Updating outlierSD to %d\n", outlierSD);
  }

  FILE* inf = fopen(argv[2], "r");

  if(!inf) {
    perror("Could not open input file");
    return -2;
  }

  char line[BUFF_SIZE + 1];

  std::vector<double> dat;
  while(fgets(line, BUFF_SIZE, inf) != NULL) {
    char* tok = getCol(line, col);
    if (tok != NULL) {
      dat.push_back(atof(tok));
    }
  }

  if (dat.size() == 0) {
    printf("%s 0 nan nan 0\n", argv[2]);
    return 0;
  }
  // SD is 0
  if (dat.size() == 1) {
    printf("%s %u %f 0.0 0\n", argv[2], dat.size(), dat[1]);
    return 0;
  }
  
  // Calculate the SD, removing outliers
  double sd = -1;
  double avg = 0;
  // Initailize so we'll keep going
  int numOutside = dat.size();
  int numPrevOutside = 0;
  int numRounds = 0;

  do {
    numRounds++;
    numPrevOutside = numOutside;
    numOutside = 0;
    double count = 0, sum = 0;
    for (double f : dat) {
      if (sd == -1 || abs(f - avg) < outlierSD * sd) {
        count++;
        sum += f;
      }
    }
    avg = sum / count;

    double sqtotal = 0;
    for (double f : dat) {
      if (sd == -1 || abs(f - avg) < outlierSD * sd) {
        sqtotal += pow(avg - f, 2);
      } else {
        numOutside++;
      }
    }

    sd = pow(sqtotal / (count - 1), 0.5);
    //printf(" - %f %f %d/%d\n", avg, sd, numPrevOutside, numOutside);
  // Keep going as long as we're removing things -- or if we've not explicitely stated to do this
  } while (numPrevOutside != numOutside && outlierSD != -1);

  printf("%s %u %f %f %d %d\n", argv[2], dat.size(), avg, sd, numOutside, numRounds);
}
