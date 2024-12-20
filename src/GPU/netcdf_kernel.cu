#ifndef _NETCDF_KERNELS_H_
#define _NETCDF_KERNELS_H_

#include "Constantes.hxx"

__global__ void obtenerEtaNetCDF(double2 *d_datosVolumenesNivel_1, float *d_vec, int num_volx,
			int num_voly, double epsilon_h, double Hmin, double H)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double2 h;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		h = d_datosVolumenesNivel_1[pos];
		d_vec[pos] = (float) ( (h.x < epsilon_h) ? -9999.0f : (h.x - h.y - Hmin)*H );
	}
}

__global__ void obtenerUxNetCDF(double2 *d_datosVolumenesNivel_1, double2 *d_datosVolumenesNivel_2, float *d_vec,
			int num_volx, int num_voly, double epsilon_h, double Hmin, double H, double Q)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double maxheps, factor;
	double h;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		h = d_datosVolumenesNivel_1[pos].x;
		if (h < epsilon_h) {
			maxheps = epsilon_h;
			factor = M_SQRT2*h / sqrt(h*h*h*h + maxheps*maxheps*maxheps*maxheps)*(Q/H);
		}
		else {
			factor = Q/(h*H);
		}
		d_vec[pos] = (float) (d_datosVolumenesNivel_2[pos].x*factor);
	}
}

__global__ void obtenerUyNetCDF(double2 *d_datosVolumenesNivel_1, double2 *d_datosVolumenesNivel_2, float *d_vec,
			int num_volx, int num_voly, double epsilon_h, double Hmin, double H, double Q)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double maxheps, factor;
	double h;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		h = d_datosVolumenesNivel_1[pos].x;
		if (h < epsilon_h) {
			maxheps = epsilon_h;
			factor = M_SQRT2*h / sqrt(h*h*h*h + maxheps*maxheps*maxheps*maxheps)*(Q/H);
		}
		else {
			factor = Q/(h*H);
		}
		d_vec[pos] = (float) (d_datosVolumenesNivel_2[pos].y*factor);
	}
}

__global__ void obtenerBatimetriaModificadaNetCDF(double2 *d_datosVolumenesNivel_1, float *d_vec, int num_volx,
			int num_voly, double Hmin, double H)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double prof;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		prof = d_datosVolumenesNivel_1[pos].y;
		d_vec[pos] = (float) ((prof+Hmin)*H);
	}
}

/***********/
/* Metrics */
/***********/

__global__ void obtenerEta1MaximaNetCDF(double *d_eta1_maxima, float *d_vec, int num_volx,
			int num_voly, double Hmin, double H)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double eta;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		eta = d_eta1_maxima[pos];
		d_vec[pos] = (float) ((eta < -1e20) ? -9999.0f : (eta-Hmin)*H);
	}
}

__global__ void obtenerTiemposLlegadaNetCDF(double *d_tiempos_llegada, float *d_vec, int num_volx,
			int num_voly, double T)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double t;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		t = d_tiempos_llegada[pos];
		d_vec[pos] = (float) ((t < 0.0) ? -1.0f : t*T);
	}
}

#endif
