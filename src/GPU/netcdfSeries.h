#ifndef NETCDF_SERIES_H
#define NETCDF_SERIES_H

#include <cstring>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>

namespace netcdfSeries
{
	extern bool ErrorEnNetCDF;
	extern int ncid_ts;
	extern int time_ts_id;
	extern int eta_ts_id;
	extern int ux_ts_id;
	extern int uy_ts_id;
	extern int eta_max_ts_id;
	extern int eta_min_ts_id;
	extern int lon_ts_id;
	extern int lat_ts_id;
	extern int grid_ts_id_okada;

	void check_err(int iret);
	void escribirDatosCroppingTimeSeriesNC(int ncid_ts, int crop_flag, double crop_value, double H);
	int initTimeSeriesNC(char *nombre_bati, char *nombre_fich, int num_points, double *lonPuntos, double *latPuntos,
		double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax, double hpos, double cvis, double dif_at,
		double borde_sup, double borde_inf, double borde_izq, double borde_der, bool es_periodica, int numFaults,
		double *defTime, double *LON_C, double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE,
		double *DIP, double *RAKE, double *SLIP, int okada_flag, int crop_flag, double crop_value, double H, char *version);
	void writeStateTimeSeriesNC(int num, float tiempo_act, int num_points, float *eta, float *ux, float *uy);
	void guardarBatimetriaModificadaTimeSeriesNC(float *vec);
	void guardarAmplitudesTimeSeriesNC(float *eta_min_puntos, float *eta_max_puntos);
	void closeTimeSeriesNC();
}

#endif
