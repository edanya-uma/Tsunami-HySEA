#ifndef NETCDFTH_H
#define NETCDFTH_H

#include "netcdfSeries.h"
#include <cstring>

using namespace netcdfSeries;

namespace netcdfTH
{
	extern int ncid;
	extern int time_id;
	extern int eta_id;
	extern int eta_max_id;
	extern int arr_times_id;
	extern int ux_id;
	extern int uy_id;
	extern int grid_id_okada;

	void abrirGRD(const char *nombre_fich, int *nvx, int *nvy);
	void leerLongitudGRD(double *lon);
	void leerLatitudGRD(double *lat);
	void leerBatimetriaGRD(float *bati);
	void leerEtaGRD(float *eta);
	void cerrarGRD();
	void escribirDatosCroppingNC(int ncid, int crop_flag, double crop_value, double H);
	void crearFicheroNC(double *lon_grid, double *lat_grid, double *lon, double *lat, char *nombre_bati, int okada_flag,
		char *nombre_fich, int *p_ncid, int *p_time_id, int *p_eta_id, int *p_ux_id, int *p_uy_id, int num_volx,
		int num_voly, double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax, double hpos,
		double cvis, double dif_at, double borde_sup, double borde_inf, double borde_izq, double borde_der,
		bool es_periodica, int numFaults, double *defTime, double *LON_C, double *LAT_C, double *DEPTH_C,
		double *FAULT_L, double *FAULT_W, double *STRIKE, double *DIP, double *RAKE, double *SLIP, int crop_flag,
		double crop_value, double H, float *bati, char *version);
	int crearFicherosNC(char *nombre_bati, int okada_flag, char *nombre_fich, int num_volx, int num_voly, double *lon_grid,
		double *lat_grid, double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax, double hpos, double cvis,
		double dif_at, double borde_sup, double borde_inf, double borde_izq, double borde_der, bool es_periodica, int numFaults,
		double *defTime, double *LON_C, double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE,
		double *DIP, double *RAKE, double *SLIP, int crop_flag, double crop_value, double H, float *bati, char *version);
	void escribirTiempoNC(int num, float tiempo_act);
	void escribirEtaNC(int num_volx, int num_voly, int num, float *eta);
	void escribirUxNC(int num_volx, int num_voly, int num, float *ux);
	void escribirUyNC(int num_volx, int num_voly, int num, float *uy);
	void guardarBatimetriaModificadaNC(float *vec);
	void guardarEta1MaximaNC(float *vec);
	void guardarTiemposLlegadaNC(float *vec);
	void cerrarFicheroNC();
}

#endif
