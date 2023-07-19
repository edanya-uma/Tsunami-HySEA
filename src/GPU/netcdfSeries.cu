#ifndef _NETCDF_SERIES_H_
#define _NETCDF_SERIES_H_

#include "Constantes.hxx"
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>

bool ErrorEnNetCDF;
int ncid_ts;
int time_ts_id;
int eta_ts_id;
int ux_ts_id;
int uy_ts_id;
int eta_max_ts_id;
int eta_min_ts_id;
int lon_ts_id;
int lat_ts_id;
int grid_ts_id_okada;

void check_err(int iret)
{
	if ((iret != NC_NOERR) && (! ErrorEnNetCDF)) {
		fprintf(stderr, "%s\n", nc_strerror(iret));
		ErrorEnNetCDF = true;
	}
}

/*********************/
/* Series de tiempos */
/*********************/

int initTimeSeriesNC(char *nombre_bati, char *prefijo, int num_points, double *lonPuntos, double *latPuntos,
		double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax, double hpos, double cvis,
		double borde_sup, double borde_inf, double borde_izq, double borde_der, double LON_C, double LAT_C,
		double DEPTH_C, double FAULT_L, double FAULT_W, double STRIKE, double DIP, double RAKE, double SLIP,
		int okada_flag)
{
	char nombre_fich[256];
	char cadena[256];
	int grid_npoints_dim;
	int grid_dims[1];
	int var_dims[2];
	int time_dim;
	int deflate_level = DEFLATE_LEVEL;
	double fill_double;
	float val_float, fill_float;
	struct timeval tv;
	char fecha_act[24];
	int iret;

	sprintf(nombre_fich, "%s_ts.nc", prefijo);
	iret = nc_create(nombre_fich, NC_CLOBBER|NC_NETCDF4, &ncid_ts);
	check_err(iret);

	iret = nc_def_dim(ncid_ts, "grid_npoints", num_points, &grid_npoints_dim);
	check_err(iret);
	iret = nc_def_dim(ncid_ts, "time", NC_UNLIMITED, &time_dim);
	check_err(iret);

	fill_double = -9999.0;
	fill_float = -9999.0f;
	grid_dims[0] = grid_npoints_dim;
	if (okada_flag == OKADA_STANDARD) {
		iret = nc_def_var(ncid_ts, "deformed_bathy", NC_FLOAT, 1, grid_dims, &grid_ts_id_okada);
		check_err(iret);
		iret = nc_def_var_deflate(ncid_ts, grid_ts_id_okada, NC_SHUFFLE, 1, deflate_level);
		check_err(iret);
	}
	else {
		iret = nc_def_var(ncid_ts, "bathymetry", NC_FLOAT, 1, grid_dims, &grid_ts_id_okada);
		check_err(iret);
		iret = nc_def_var_deflate(ncid_ts, grid_ts_id_okada, NC_SHUFFLE, 1, deflate_level);
		check_err(iret);
	}
	var_dims[0] = time_dim;
	var_dims[1] = grid_npoints_dim;
	grid_dims[0] = grid_npoints_dim;
	iret = nc_def_var(ncid_ts, "time", NC_FLOAT, 1, &time_dim, &time_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, time_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "longitude", NC_DOUBLE, 1, grid_dims, &lon_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, lon_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, lon_ts_id, "units", 12, "degrees_east");
	check_err(iret);
	iret = nc_put_att_double(ncid_ts, lon_ts_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid_ts, lon_ts_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "latitude", NC_DOUBLE, 1, grid_dims, &lat_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, lat_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, lat_ts_id, "units", 13, "degrees_north");
	check_err(iret);
	iret = nc_put_att_double(ncid_ts, lat_ts_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid_ts, lat_ts_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "min_height", NC_FLOAT, 1, grid_dims, &eta_min_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, eta_min_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_min_ts_id, "long_name", 22, "Minimum wave amplitude");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_min_ts_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_min_ts_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_min_ts_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "max_height", NC_FLOAT, 1, grid_dims, &eta_max_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, eta_max_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_max_ts_id, "long_name", 22, "Maximum wave amplitude");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_max_ts_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_max_ts_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_max_ts_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "eta", NC_FLOAT, 2, var_dims, &eta_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, eta_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_ts_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, eta_ts_id, "long_name", 14, "Wave amplitude");
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_ts_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, eta_ts_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "ux", NC_FLOAT, 2, var_dims, &ux_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, ux_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, ux_ts_id, "units", 13, "meters/second");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, ux_ts_id, "long_name", 33, "Velocity of water along longitude");
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, ux_ts_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, ux_ts_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_def_var(ncid_ts, "uy", NC_FLOAT, 2, var_dims, &uy_ts_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid_ts, uy_ts_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, uy_ts_id, "units", 13, "meters/second");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, uy_ts_id, "long_name", 32, "Velocity of water along latitude");
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, uy_ts_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid_ts, uy_ts_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	if (okada_flag == OKADA_STANDARD) {
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "long_name", 24, "Grid deformed bathymetry");
		check_err(iret);
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "standard_name", 14, "deformed depth");
		check_err(iret);
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "units", 6, "meters");
		check_err(iret);
		iret = nc_put_att_float(ncid_ts, grid_ts_id_okada, "missing_value", NC_FLOAT, 1, &fill_float);
		check_err(iret);
		iret = nc_put_att_float(ncid_ts, grid_ts_id_okada, "_FillValue", NC_FLOAT, 1, &fill_float);
		check_err(iret);
	}
	else {
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "long_name", 15, "Grid bathymetry");
		check_err(iret);
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "standard_name", 5, "depth");
		check_err(iret);
		iret = nc_put_att_text(ncid_ts, grid_ts_id_okada, "units", 6, "meters");
		check_err(iret);
		iret = nc_put_att_float(ncid_ts, grid_ts_id_okada, "missing_value", NC_FLOAT, 1, &fill_float);
		check_err(iret);
		iret = nc_put_att_float(ncid_ts, grid_ts_id_okada, "_FillValue", NC_FLOAT, 1, &fill_float);
		check_err(iret);
	}

	iret = nc_put_att_text(ncid_ts, time_ts_id, "long_name", 4, "Time");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, time_ts_id, "units", 24, "seconds since 1970-01-01");
	check_err(iret);

	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "Conventions", 6, "CF-1.0");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "title", 40, "Time series output of TsunamiHySEA model");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "Tsunami-HySEA_version", 16, "open source v1.1");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "creator_name", 12, "EDANYA Group");
	check_err(iret);
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "institution", 20, "University of Malaga");
	check_err(iret);
	sprintf(cadena, " ");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "comments", strlen(cadena), cadena);
	check_err(iret);
	sprintf(cadena, " ");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "references", strlen(cadena), cadena);
	check_err(iret);

	gettimeofday(&tv, NULL);
	strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
	check_err(iret);

	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "grid_name", strlen(nombre_bati), nombre_bati);
	check_err(iret);
	val_float = (float) tiempo_tot;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "simulation_time", NC_FLOAT, 1, &val_float);
	check_err(iret);

	sprintf(cadena, (fabs(borde_sup-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "upper_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_inf-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "lower_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_izq-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "left_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_der-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "right_border", 4, cadena);
	check_err(iret);

	val_float = CFL;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "CFL", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = epsilon_h;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "epsilon_h", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = mf0;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "water_bottom_friction", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = vmax;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "max_velocity_water", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = hpos;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "threshold_2swaf", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = cvis;
	iret = nc_put_att_float(ncid_ts, NC_GLOBAL, "stability_coefficient", NC_FLOAT, 1, &val_float);
	check_err(iret);

	if (okada_flag == SEA_SURFACE_FROM_FILE) {
		iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "initialization_mode", 21, "sea_surface_from_file");
		check_err(iret);
	}
	else if (okada_flag == OKADA_STANDARD) {
		iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "initialization_mode", 14, "okada_standard");
		check_err(iret);
		sprintf(cadena, "lon: %.4f, lat: %.4f, depth: %.2f, length: %.2f, width: %.2f, strike: %.2f, dip: %.2f, rake: %.2f, slip: %.2f",
			LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP);
		iret = nc_put_att_text(ncid_ts, NC_GLOBAL, "fault", strlen(cadena), cadena);
		check_err(iret);
	}

	iret = nc_enddef(ncid_ts);
	check_err(iret);

	iret = nc_put_var_double(ncid_ts, lon_ts_id, lonPuntos);
	check_err(iret);
	iret = nc_put_var_double(ncid_ts, lat_ts_id, latPuntos);
	check_err(iret);

	return 0;
}

void writeStateTimeSeriesNC(int num, float tiempo_act, int num_points, float *eta, float *ux, float *uy)
{
	int iret;
	float t_act = tiempo_act;
	const size_t num_cst = num;
	const size_t npoints_cst = num_points;
	const size_t uno = 1;
	const size_t start[] = {num_cst, 0};
	const size_t count[] = {1, npoints_cst};

	iret = nc_put_vara_float(ncid_ts, time_ts_id, &num_cst, &uno, &t_act);
	check_err(iret);

	iret = nc_put_vara_float(ncid_ts, eta_ts_id, start, count, eta);
	check_err(iret);

	iret = nc_put_vara_float(ncid_ts, ux_ts_id, start, count, ux);
	check_err(iret);

	iret = nc_put_vara_float(ncid_ts, uy_ts_id, start, count, uy);
	check_err(iret);

	iret = nc_sync(ncid_ts);
	check_err(iret);
}

void guardarBatimetriaModificadaTimeSeriesNC(float *vec)
{
	int iret;

	iret = nc_put_var_float(ncid_ts, grid_ts_id_okada, vec);
	check_err(iret);
}

void guardarAmplitudesTimeSeriesNC(float *eta_min_puntos, float *eta_max_puntos)
{
	int iret;

	iret = nc_put_var_float(ncid_ts, eta_min_ts_id, eta_min_puntos);
	check_err(iret);
	iret = nc_put_var_float(ncid_ts, eta_max_ts_id, eta_max_puntos);
	check_err(iret);
}

void closeTimeSeriesNC()
{
	int iret;

	iret = nc_close(ncid_ts);
	check_err(iret);
}

#endif
