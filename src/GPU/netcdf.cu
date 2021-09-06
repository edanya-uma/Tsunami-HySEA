#include "Constantes.hxx"
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>

bool ErrorEnNetCDF;
int ncid;
int time_id;
int eta_id, eta_max_id;
int qx_id, qy_id;

void check_err(int iret)
{
	if ((iret != NC_NOERR) && (! ErrorEnNetCDF)) {
		fprintf(stderr, "%s\n", nc_strerror(iret));
		ErrorEnNetCDF = true;
	}
}

void crearFicheroNC(double *lon_grid, double *lat_grid, double *lon, double *lat, char *nombre_bati,
			char *prefijo, int *p_ncid, int *p_time_id, int *p_eta_id, int *p_qx_id, int *p_qy_id,
			int num_volx, int num_voly, double tiempo_tot, double CFL, double epsilon_h, double mf0,
			double vmax, double hpos, double cvis, double *bati)
{
	char nombre_fich[256];
	char cadena[256];
	int grid_lon_dim, grid_lat_dim;
	int grid_dims[2];
	int lon_dim, lat_dim;
	int var_dims[3];
	int time_dim;
	int ncid;
	int grid_id, grid_lon_id, grid_lat_id;
	int lon_id, lat_id;
	int val_int;
	double val_double, fill_double;
	struct timeval tv;
	char fecha_act[24];
	int iret;

	sprintf(nombre_fich, "%s.nc", prefijo);
	iret = nc_create(nombre_fich, NC_CLOBBER, p_ncid);
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
	iret = nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dim, &lat_id);
	check_err(iret);
	iret = nc_def_var(ncid, "grid_lon", NC_DOUBLE, 1, &grid_lon_dim, &grid_lon_id);
	check_err(iret);
	iret = nc_def_var(ncid, "grid_lat", NC_DOUBLE, 1, &grid_lat_dim, &grid_lat_id);
	check_err(iret);
	grid_dims[0] = grid_lat_dim;
	grid_dims[1] = grid_lon_dim;
	iret = nc_def_var(ncid, "bathymetry", NC_DOUBLE, 2, grid_dims, &grid_id);
	check_err(iret);
	var_dims[0] = time_dim;
	var_dims[1] = lat_dim;
	var_dims[2] = lon_dim;
	grid_dims[0] = lat_dim;
	grid_dims[1] = lon_dim;
	iret = nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, p_time_id);
	check_err(iret);
	iret = nc_def_var(ncid, "max_height", NC_DOUBLE, 2, grid_dims, &eta_max_id);
	check_err(iret);
	iret = nc_def_var(ncid, "eta", NC_DOUBLE, 3, var_dims, p_eta_id);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_eta_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_eta_id, "long_name", 14, "Wave amplitude");
	check_err(iret);
	iret = nc_def_var(ncid, "qx", NC_DOUBLE, 3, var_dims, p_qx_id);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_qx_id, "units", 15, "meters^2/second");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_qx_id, "long_name", 34, "Mass flow of water along longitude");
	check_err(iret);
	iret = nc_def_var(ncid, "qy", NC_DOUBLE, 3, var_dims, p_qy_id);
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_qy_id, "units", 15, "meters^2/second");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_qy_id, "long_name", 33, "Mass flow of water along latitude");
	check_err(iret);

	fill_double = -9999.0;
	iret = nc_put_att_text(ncid, grid_id, "long_name", 15, "Grid bathymetry");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_id, "standard_name", 5, "depth");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_double(ncid, grid_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, grid_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);

	iret = nc_put_att_text(ncid, grid_lon_id, "long_name", 14, "Grid longitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lon_id, "units", 7, "degrees");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lat_id, "long_name", 13, "Grid latitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, grid_lat_id, "units", 7, "degrees");
	check_err(iret);

	iret = nc_put_att_text(ncid, eta_max_id, "long_name", 22, "Maximum wave amplitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, eta_max_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_double(ncid, eta_max_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, eta_max_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);

	iret = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lon_id, "units", 7, "degrees");
	check_err(iret);

	iret = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
	check_err(iret);
	iret = nc_put_att_text(ncid, lat_id, "units", 7, "degrees");
	check_err(iret);

	iret = nc_put_att_text(ncid, *p_time_id, "long_name", 4, "Time");
	check_err(iret);
	iret = nc_put_att_text(ncid, *p_time_id, "units", 24, "seconds since 1970-01-01");
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_eta_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_eta_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_qx_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_qx_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_qy_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *p_qy_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);

	iret = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "title", 25, "TsunamiHySEA model output");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "creator_name", 12, "EDANYA Group");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "institution", 20, "University of Malaga");
	check_err(iret);
	sprintf(cadena, "No paper");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "comments", strlen(cadena), cadena);
	check_err(iret);
	sprintf(cadena, "http://path.to.paper/paper.pdf");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "references", strlen(cadena), cadena);
	check_err(iret);

	gettimeofday(&tv, NULL);
	strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
	iret = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
	check_err(iret);

	iret = nc_put_att_text(ncid, NC_GLOBAL, "grid_name", strlen(nombre_bati), nombre_bati);
	check_err(iret);
	val_double = tiempo_tot;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "simulation_time", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_int = 1;
	iret = nc_put_att_int(ncid, NC_GLOBAL, "output_grid_interval", NC_INT, 1, &val_int);
	check_err(iret);
	val_double = CFL;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "CFL", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = epsilon_h;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "epsilon_h", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = mf0;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "water_bottom_friction", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = vmax;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "max_velocity_water", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = hpos;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "threshold_2swaf", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = cvis;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "stability_coefficient", NC_DOUBLE, 1, &val_double);
	check_err(iret);

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
	iret = nc_put_var_double(ncid, grid_id, bati);
	check_err(iret);
}

int crearFicherosNC(char *nombre_bati, char *prefijo, int num_volx, int num_voly, double *lon_grid,
			double *lat_grid, double tiempo_tot, double CFL, double epsilon_h, double mf0,
			double vmax, double hpos, double cvis, double *bati)
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

	crearFicheroNC(lon_grid, lat_grid, lon, lat, nombre_bati, prefijo, &ncid, &time_id, &eta_id,
		&qx_id, &qy_id, num_volx, num_voly, tiempo_tot, CFL, epsilon_h, mf0, vmax, hpos, cvis, bati);

	free(lon);
	free(lat);

	return 0;
}

void escribirTiempoNC(int num, double tiempo_act)
{
	int iret;
	double t_act = tiempo_act;
	const size_t paso = num;
	const size_t uno = 1;

	iret = nc_put_vara_double(ncid, time_id, &paso, &uno, &t_act);
	check_err(iret);
}

void escribirEtaNC(int num_volx, int num_voly, int num, double *eta)
{
	int iret;
	const size_t var_start[3] = {num, 0, 0};
	const size_t var_count[3] = {1, num_voly, num_volx};

	iret = nc_put_vara_double(ncid, eta_id, var_start, var_count, eta);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void escribirQxNC(int num_volx, int num_voly, int num, double *qx)
{
	int iret;
	const size_t var_start[3] = {num, 0, 0};
	const size_t var_count[3] = {1, num_voly, num_volx};

	iret = nc_put_vara_double(ncid, qx_id, var_start, var_count, qx);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void escribirQyNC(int num_volx, int num_voly, int num, double *qy)
{
	int iret;
	const size_t var_start[3] = {num, 0, 0};
	const size_t var_count[3] = {1, num_voly, num_volx};

	iret = nc_put_vara_double(ncid, qy_id, var_start, var_count, qy);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}

void cerrarFicheroNC(double *eta_max)
{
	int iret;

	iret = nc_put_var_double(ncid, eta_max_id, eta_max);
	check_err(iret);
	iret = nc_close(ncid);
	check_err(iret);
}

