#ifndef _CONSTANTES_H_
#define _CONSTANTES_H_

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define EPSILON             DBL_EPSILON
#define EARTH_RADIUS        6378136.6

using namespace std;

#define NUM_HEBRASX_ARI 8
#define NUM_HEBRASY_ARI 8
#define NUM_HEBRASX_EST 8
#define NUM_HEBRASY_EST 8

#define iDivUp(a,b)  (((a)%(b) != 0) ? ((a)/(b) + 1) : ((a)/(b)))

typedef struct {
	double2 *areaYCosPhi;
	double *anchoVolumenes;
	double *altoVolumenes;
	double *longitud;
	double *latitud;
} tipoDatosSubmalla;

#endif

