#ifndef _NETCDF_H_
#define _NETCDF_H_

#include "Constantes.hxx"
#include "netcdfSeries.cu"

int ncid;
int time_id;
int eta_id, eta_max_id;
int ux_id, uy_id;
int grid_id_okada;

void abrirGRD(const char *nombre_fich, int *nvx, int *nvy)
{
	int iret;
	int lon_dim, lat_dim;
	size_t nx, ny;

	iret = nc_open(nombre_fich, NC_NOWRITE, &ncid);
	check_err(iret);
	iret = nc_inq_dimid(ncid, "lon", &lon_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(ncid, "x", &lon_dim);
		check_err(iret);
	}
	iret = nc_inq_dimid(ncid, "lat", &lat_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(ncid, "y", &lat_dim);
		check_err(iret);
	}
	iret = nc_inq_dimlen(ncid, lon_dim, &nx);
	check_err(iret);
	iret = nc_inq_dimlen(ncid, lat_dim, &ny);
	check_err(iret);
	*nvx = nx;
	*nvy = ny;
}

void leerLongitudGRD(double *lon)
{
	int iret, lon_id;

	iret = nc_inq_varid(ncid, "lon", &lon_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(ncid, "x", &lon_id);
		check_err(iret);
	}
	iret = nc_get_var_double(ncid, lon_id, lon);
	check_err(iret);
}

void leerLatitudGRD(double *lat)
{
	int iret, lat_id;

	iret = nc_inq_varid(ncid, "lat", &lat_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(ncid, "y", &lat_id);
		check_err(iret);
	}
	iret = nc_get_var_double(ncid, lat_id, lat);
	check_err(iret);
}

void leerBatimetriaGRD(float *bati)
{
	int iret, z_id;

	iret = nc_inq_varid(ncid, "z", &z_id);
	check_err(iret);
	iret = nc_get_var_float(ncid, z_id, bati);
	check_err(iret);
}

void leerEtaGRD(float *eta)
{
	int iret, var_id;

	iret = nc_inq_varid(ncid, "z", &var_id);
	check_err(iret);
	iret = nc_get_var_float(ncid, var_id, eta);
	check_err(iret);
}

void cerrarGRD()
{
	int iret;

	iret = nc_close(ncid);
	check_err(iret);
}

void crearFicheroNC(double *lon_grid, double *lat_grid, double *lon, double *lat, char *nombre_bati, int okada_flag,
			char *prefijo, int *p_ncid, int *p_time_id, int *p_eta_id, int *p_ux_id, int *p_uy_id, int num_volx,
			int num_voly, double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax, double hpos,
			double cvis, double borde_sup, double borde_inf, double borde_izq, double borde_der, double LON_C,
			double LAT_C, double DEPTH_C, double FAULT_L, double FAULT_W, double STRIKE, double DIP,
			double RAKE, double SLIP, float *bati)
{
	char nombre_fich[256];
	char cadena[256];
	int grid_lon_dim, grid_lat_dim;
	int grid_dims[2];
	int lon_dim, lat_dim;
	int var_dims[3];
	int time_dim;
	int deflate_level = DEFLATE_LEVEL;
	int ncid;
	int grid_id, grid_lon_id, grid_lat_id;
	int lon_id, lat_id;
	int val_int;
	float val_float, fill_float;
	struct timeval tv;
	char fecha_act[24];
	int iret;

	sprintf(nombre_fich, "%s.nc", prefijo);
	iret = nc_create(nombre_fich, NC_CLOBBER|NC_NETCDF4, p_ncid);
	check_err(iret);
	ncid = *p_ncid;

	iret = nc_def_dim(ncid, "lon", num_volx, &lon_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "lat", num_voly, &lat_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "grid_lon", num_volx, &grid_lon_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "grid_lat", num_voly, &grid_lat_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	check_err(iret);

	iret = nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dim, &lon_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, lon_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dim, &lat_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, lat_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_def_var(ncid, "grid_lon", NC_DOUBLE, 1, &grid_lon_dim, &grid_lon_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, grid_lon_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_def_var(ncid, "grid_lat", NC_DOUBLE, 1, &grid_lat_dim, &grid_lat_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, grid_lat_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	grid_dims[0] = grid_lat_dim;
	grid_dims[1] = grid_lon_dim;
	if (okada_flag == OKADA_STANDARD) {
		iret = nc_def_var(ncid, "original_bathy", NC_FLOAT, 2, grid_dims, &grid_id);
		check_err(iret);
		iret = nc_def_var_deflate(ncid, grid_id, NC_SHUFFLE, 1, deflate_level);
		check_err(iret);
		iret = nc_def_var(ncid, "deformed_bathy", NC_FLOAT, 2, grid_dims, &grid_id_okada);
		check_err(iret);
		iret = nc_def_var_deflate(ncid, grid_id_okada, NC_SHUFFLE, 1, deflate_level);
		check_err(iret);
	}
	else {
		iret = nc_def_var(ncid, "bathymetry", NC_FLOAT, 2, grid_dims, &grid_id);
		check_err(iret);
		iret = nc_def_var_deflate(ncid, grid_id, NC_SHUFFLE, 1, deflate_level);
		check_err(iret);
	}
	var_dims[0] = time_dim;
	var_dims[1] = lat_dim;
	var_dims[2] = lon_dim;
	grid_dims[0] = lat_dim;
	grid_dims[1] = lon_dim;
	iret = nc_def_var(ncid, "time", NC_FLOAT, 1, &time_dim, p_time_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, *p_time_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_def_var(ncid, "max_height", NC_FLOAT, 2, grid_dims, &eta_max_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, eta_max_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_def_var(ncid, "eta", NC_FLOAT, 3, var_dims, p_eta_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, *p_eta_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_eta_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_eta_id, "long_name", 14, "Wave amplitude");
	check_err(iret);
	iret = nc_def_var(ncid, "ux", NC_FLOAT, 3, var_dims, p_ux_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, *p_ux_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_ux_id, "units", 13, "meters/second");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_ux_id, "long_name", 33, "Velocity of water along longitude");
	check_err(iret);
	iret = nc_def_var(ncid, "uy", NC_FLOAT, 3, var_dims, p_uy_id);
	check_err(iret);
	iret = nc_def_var_deflate(ncid, *p_uy_id, NC_SHUFFLE, 1, deflate_level);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_uy_id, "units", 13, "meters/second");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_uy_id, "long_name", 32, "Velocity of water along latitude");
	check_err(iret);

	fill_float = -9999.0;
	if (okada_flag == OKADA_STANDARD) {
		iret = nc_put_att_text(ncid, grid_id, "long_name", 24, "Grid original bathymetry");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id, "standard_name", 14, "original depth");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id, "units", 6, "meters");
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id, "missing_value", NC_FLOAT, 1, &fill_float);
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id, "_FillValue", NC_FLOAT, 1, &fill_float);
		check_err(iret);

		iret = nc_put_att_text(ncid, grid_id_okada, "long_name", 24, "Grid deformed bathymetry");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id_okada, "standard_name", 14, "deformed depth");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id_okada, "units", 6, "meters");
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id_okada, "missing_value", NC_FLOAT, 1, &fill_float);
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id_okada, "_FillValue", NC_FLOAT, 1, &fill_float);
		check_err(iret);
	}
	else {
		iret = nc_put_att_text(ncid, grid_id, "long_name", 15, "Grid bathymetry");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id, "standard_name", 5, "depth");
		check_err(iret);
		iret = nc_put_att_text(ncid, grid_id, "units", 6, "meters");
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id, "missing_value", NC_FLOAT, 1, &fill_float);
		check_err(iret);
		iret = nc_put_att_float(ncid, grid_id, "_FillValue", NC_FLOAT, 1, &fill_float);
		check_err(iret);
	}

	iret = nc_put_att_text(ncid, grid_lon_id, "long_name", 14, "Grid longitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lon_id, "units", 12, "degrees_east");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lat_id, "long_name", 13, "Grid latitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lat_id, "units", 13, "degrees_north");
	check_err(iret);

	iret = nc_put_att_text(ncid, eta_max_id, "long_name", 22, "Maximum wave amplitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, eta_max_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_float(ncid, eta_max_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, eta_max_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
	check_err(iret);

	iret = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
	check_err(iret);

	iret = nc_put_att_text(ncid, *p_time_id, "long_name", 4, "Time");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_time_id, "units", 24, "seconds since 1970-01-01");
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_eta_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_eta_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_ux_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_ux_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_uy_id, "missing_value", NC_FLOAT, 1, &fill_float);
	check_err(iret);
	iret = nc_put_att_float(ncid, *p_uy_id, "_FillValue", NC_FLOAT, 1, &fill_float);
	check_err(iret);

	iret = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "title", 25, "TsunamiHySEA model output");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "Tsunami-HySEA_version", 16, "open source v1.1");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "creator_name", 12, "EDANYA Group");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "institution", 20, "University of Malaga");
	check_err(iret);
	sprintf(cadena, " ");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "comments", strlen(cadena), cadena);
	check_err(iret);
	sprintf(cadena, " ");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "references", strlen(cadena), cadena);
	check_err(iret);

	gettimeofday(&tv, NULL);
	strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
	iret = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
	check_err(iret);

	iret = nc_put_att_text(ncid, NC_GLOBAL, "grid_name", strlen(nombre_bati), nombre_bati);
	check_err(iret);
	val_float = tiempo_tot;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "simulation_time", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_int = 1;
	iret = nc_put_att_int(ncid, NC_GLOBAL, "output_grid_interval", NC_INT, 1, &val_int);
	check_err(iret);

	sprintf(cadena, (fabs(borde_sup-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "upper_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_inf-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "lower_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_izq-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "left_border", 4, cadena);
	check_err(iret);
	sprintf(cadena, (fabs(borde_der-1.0) < EPSILON) ? "open" : "wall");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "right_border", 4, cadena);
	check_err(iret);

	val_float = CFL;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "CFL", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = epsilon_h;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "epsilon_h", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = mf0;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "water_bottom_friction", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = vmax;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "max_velocity_water", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = hpos;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "threshold_2swaf", NC_FLOAT, 1, &val_float);
	check_err(iret);
	val_float = cvis;
	iret = nc_put_att_float(ncid, NC_GLOBAL, "stability_coefficient", NC_FLOAT, 1, &val_float);
	check_err(iret);

	if (okada_flag == SEA_SURFACE_FROM_FILE) {
		iret = nc_put_att_text(ncid, NC_GLOBAL, "initialization_mode", 21, "sea_surface_from_file");
		check_err(iret);
	}
	else if (okada_flag == OKADA_STANDARD) {
		iret = nc_put_att_text(ncid, NC_GLOBAL, "initialization_mode", 14, "okada_standard");
		check_err(iret);
		sprintf(cadena, "lon: %.4f, lat: %.4f, depth: %.4f, length: %.4f, width: %.4f, strike: %.4f, dip: %.4f, rake: %.4f, slip: %.4f",
			LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP);
		iret = nc_put_att_text(ncid, NC_GLOBAL, "fault", strlen(cadena), cadena);
		check_err(iret);
	}

	iret = nc_enddef(ncid);
	check_err(iret);

	iret = nc_put_var_double(ncid, lon_id, lon);
	check_err(iret);
	iret = nc_put_var_double(ncid, lat_id, lat);
	check_err(iret);

	iret = nc_put_var_double(ncid, grid_lon_id, lon_grid);
	check_err(iret);
	iret = nc_put_var_double(ncid, grid_lat_id, lat_grid);
	check_err(iret);
	iret = nc_put_var_float(ncid, grid_id, bati);
	check_err(iret);
}

int crearFicherosNC(char *nombre_bati, int okada_flag, char *prefijo, int num_volx, int num_voly, double *lon_grid,
			double *lat_grid, double tiempo_tot, double CFL, double epsilon_h, double mf0, double vmax,
			double hpos, double cvis, double borde_sup, double borde_inf, double borde_izq, double borde_der,
			double LON_C, double LAT_C, double DEPTH_C, double FAULT_L, double FAULT_W, double STRIKE,
			double DIP, double RAKE, double SLIP, float *bati)
{
	double *lon, *lat;
	int i;

	ErrorEnNetCDF = false;
	lon = (double *) malloc(num_volx*sizeof(double));
	lat = (double *) malloc(num_voly*sizeof(double));
	if (lat == NULL) {
		if (lon != NULL) free(lon);
		return 1;
	}

	for (i=0; i<num_volx; i++)
		lon[i] = lon_grid[i];
	for (i=0; i<num_voly; i++)
		lat[i] = lat_grid[i];

	crearFicheroNC(lon_grid, lat_grid, lon, lat, nombre_bati, okada_flag, prefijo, &ncid, &time_id, &eta_id,
		&ux_id, &uy_id, num_volx, num_voly, tiempo_tot, CFL, epsilon_h, mf0, vmax, hpos, cvis, borde_sup,
		borde_inf, borde_izq, borde_der, LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP, bati);

	free(lon);
	free(lat);

	return 0;
}

void escribirTiempoNC(int num, float tiempo_act)
{
	int iret;
	float t_act = tiempo_act;
	const size_t paso = num;
	const size_t uno = 1;

	iret = nc_put_vara_float(ncid, time_id, &paso, &uno, &t_act);
	check_err(iret);
}

void escribirEtaNC(int num_volx, int num_voly, int num, float *eta)
{
	int iret;
	const size_t var_start[3] = {(size_t) num, 0, 0};
	const size_t var_count[3] = {1, (size_t) num_voly, (size_t) num_volx};

	iret = nc_put_vara_float(ncid, eta_id, var_start, var_count, eta);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void escribirUxNC(int num_volx, int num_voly, int num, float *ux)
{
	int iret;
	const size_t var_start[3] = {(size_t) num, 0, 0};
	const size_t var_count[3] = {1, (size_t) num_voly, (size_t) num_volx};

	iret = nc_put_vara_float(ncid, ux_id, var_start, var_count, ux);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void escribirUyNC(int num_volx, int num_voly, int num, float *uy)
{
	int iret;
	const size_t var_start[3] = {(size_t) num, 0, 0};
	const size_t var_count[3] = {1, (size_t) num_voly, (size_t) num_volx};

	iret = nc_put_vara_float(ncid, uy_id, var_start, var_count, uy);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void guardarBatimetriaModificadaNC(float *vec)
{
	int iret;

	iret = nc_put_var_float(ncid, grid_id_okada, vec);
	check_err(iret);
}

void cerrarFicheroNC(float *eta_max)
{
	int iret;

	iret = nc_put_var_float(ncid, eta_max_id, eta_max);
	check_err(iret);
	iret = nc_close(ncid);
	check_err(iret);
}

#endif
