#ifndef CONSTANTES_HXX
#define CONSTANTES_HXX

#include <iostream>
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>

constexpr double EPSILON       = std::numeric_limits<double>::epsilon();
constexpr double EARTH_RADIUS  = 6378136.6;
constexpr int    MAX_FAULTS    = 5000;   // Maximum number of Okada faults
constexpr int    SPONGE_SIZE   = 4;      // Size of the sponge layer (>= 0)
constexpr double SEA_LEVEL     = 0.0;    // Sea level in meters. Used in sponge layer
constexpr int    DEFLATE_LEVEL = 5;      // Level of compression of the NetCDF files (0-9)

constexpr int SEA_SURFACE_FROM_FILE = 0;
constexpr int OKADA_STANDARD        = 1;
constexpr int NO_CROP               = 0;
constexpr int CROP_RELATIVE         = 1;
constexpr int CROP_ABSOLUTE         = 2;

constexpr int NUM_HEBRASX_ARI   = 8;
constexpr int NUM_HEBRASY_ARI   = 8;
constexpr int NUM_HEBRASX_EST   = 8;
constexpr int NUM_HEBRASY_EST   = 8;
constexpr int NUM_HEBRAS_PUNTOS = 128;

constexpr int iDivUp(int a, int b) { return 1 + (a-1)/b; }

using tipoDatosSubmalla = struct {
	double2 *areaYCosPhi;
	double *anchoVolumenes;
	double *altoVolumenes;
	double *longitud;
	double *latitud;
};

#endif

