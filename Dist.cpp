#include <iostream>
#include <cmath>
#include <cfloat>
#include <time.h>
#include "nr3.h"
#include "ran.h"


using namespace std;
Ran ran(time(NULL));

double generateGaussianNoise(const double& mean, const double &stdDev);


int main ( int argc, char *argv[] ){
  int N = atoi(argv[1]);  /* Sample numbers */
  double mean = 0.5;    /* Mean of Gaussian distribution*/
  double sv = 0.04;    /* Standart variation of Gaussian distribution*/
  FILE* outfile = fopen("gaussian.dat", "w");
  FILE* outfile2 = fopen("uniform.dat", "w");
  for (int i = 0; i < N; i++){
    fprintf(outfile, "%.5lf\n", generateGaussianNoise(mean, sv));
    fprintf(outfile2, "%.5lf\n", ran.doub());
  }
  fclose(outfile);
  fclose(outfile2);
  return 0;
}

double generateGaussianNoise(const double& mean, const double &stdDev) { /* Marsaglia polar method */
  static bool hasSpare = false;
  static double spare;

  if(hasSpare) {
  	hasSpare = false;
  	return mean + stdDev * spare;
  }
  hasSpare = true;
  static double u, v, s;
  do {
   	u = ran.doub() * 2.0 - 1.0;
  	v = ran.doub() * 2.0 - 1.0;
    s = u * u + v * v;
  } while( (s >= 1.0) || (s == 0.0) );
  s = sqrt(-2.0 * log(s) / s);
  spare = v * s;
  return mean + stdDev * u * s;
}
