#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <fstream>
#include "Timestep.cu"
#include "helper_timer.h"
#include "Deformacion.cu"
#include "Metrics.cu"
#include "netcdf_kernel.cu"
#include "netcdf.cu"

/*****************/
/* NetCDF saving */
/*****************/

void guardarNetCDFNivel0(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2, float *d_vec, float *vec,
			int num, int num_volx, int num_voly, double Hmin, double tiempo_act, double epsilon_h,
			dim3 blockGridEstNivel, dim3 threadBlockEst, double H, double Q, double T)
{
	double tiempo = tiempo_act*T;
	int tam_datosVolFloat = num_volx*num_voly*sizeof(float);

	escribirTiempoNC(num, tiempo);
	obtenerEtaNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenes_1, d_vec, num_volx,
		num_voly, epsilon_h, Hmin, H);
	cudaMemcpy(vec, d_vec, tam_datosVolFloat, cudaMemcpyDeviceToHost);
	escribirEtaNC(num_volx, num_voly, num, vec);

	obtenerUxNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenes_1, d_datosVolumenes_2,
		d_vec, num_volx, num_voly, epsilon_h, Hmin, H, Q);
	cudaMemcpy(vec, d_vec, tam_datosVolFloat, cudaMemcpyDeviceToHost);
	escribirUxNC(num_volx, num_voly, num, vec);

	obtenerUyNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenes_1, d_datosVolumenes_2,
		d_vec, num_volx, num_voly, epsilon_h, Hmin, H, Q);
	cudaMemcpy(vec, d_vec, tam_datosVolFloat, cudaMemcpyDeviceToHost);
	escribirUyNC(num_volx, num_voly, num, vec);
}

/***************/
/* Time series */
/***************/

// Format of the NetCDF time series file:
// For each point:
//   <longitude>
//   <latitude>
//   <bathymetry (original if okada_flag is INITIAL_FROM_FILE; deformed if okada_flag is OKADA_STANDARD)
//   <minimum eta>
//   <maximum eta>
// For each time:
//   For each point:
//     <eta point 1> <u point 1> <v point 1> ... <eta point n> <u point n> <v point n>
void obtenerBatimetriaParaSerieTiempos(double2 *datosVolumenesNivel_1, int numPuntosGuardar,
		int *posicionesVolumenesGuardado, float *profPuntosGuardado, double Hmin, double H)
{
	int i;

	for (i=0; i<numPuntosGuardar; i++) {
		if (posicionesVolumenesGuardado[i] != -1)
			profPuntosGuardado[i] = (float) ((datosVolumenesNivel_1[i].y + Hmin)*H);
	}
}

void guardarSerieTiemposNivel0(double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2,
			int numPuntosGuardar, int *posicionesVolumenesGuardado, float *etaPuntosGuardado,
			float *uPuntosGuardado, float *vPuntosGuardado, float *etaMinPuntosGuardado,
			float *etaMaxPuntosGuardado, double Hmin, int num_ts, double tiempo_act,
			double H, double Q, double T)
{
	double tiempo = tiempo_act*T;
	double h;
	int i;

	for (i=0; i<numPuntosGuardar; i++) {
		if (posicionesVolumenesGuardado[i] != -1) {
			h = datosVolumenesNivel_1[i].x;
			etaPuntosGuardado[i] = (float) ((h - datosVolumenesNivel_1[i].y - Hmin)*H);
			uPuntosGuardado[i] = (float) (datosVolumenesNivel_2[i].x*Q/H);
			vPuntosGuardado[i] = (float) (datosVolumenesNivel_2[i].y*Q/H);
			etaMinPuntosGuardado[i] = min(etaPuntosGuardado[i], etaMinPuntosGuardado[i]);
			etaMaxPuntosGuardado[i] = max(etaPuntosGuardado[i], etaMaxPuntosGuardado[i]);
		}
	}
	writeStateTimeSeriesNC(num_ts, tiempo, numPuntosGuardar, etaPuntosGuardado, uPuntosGuardado, vPuntosGuardado);
}

/**************/
/* Block size */
/**************/

void obtenerTamBloquesKernel(int num_volx, int num_voly, dim3 *blockGridVer1, dim3 *blockGridVer2,
							dim3 *blockGridHor1, dim3 *blockGridHor2, dim3 *blockGridEst)
{
	int num_aristas_ver1, num_aristas_ver2;
	int num_aristas_hor1, num_aristas_hor2;

	num_aristas_ver1 = (num_volx/2 + 1)*num_voly;
	num_aristas_ver2 = ((num_volx&1) == 0) ? num_volx*num_voly/2 : num_aristas_ver1;
	num_aristas_hor1 = (num_voly/2 + 1)*num_volx;
	num_aristas_hor2 = ((num_voly&1) == 0) ? num_volx*num_voly/2 : num_aristas_hor1;

	blockGridVer1->x = iDivUp(num_aristas_ver1/num_voly, NUM_HEBRASX_ARI);
	blockGridVer1->y = iDivUp(num_voly, NUM_HEBRASY_ARI);
	blockGridVer2->x = iDivUp(num_aristas_ver2/num_voly, NUM_HEBRASX_ARI);
	blockGridVer2->y = iDivUp(num_voly, NUM_HEBRASY_ARI);

	blockGridHor1->x = iDivUp(num_volx, NUM_HEBRASX_ARI);
	blockGridHor1->y = iDivUp(num_aristas_hor1/num_volx, NUM_HEBRASY_ARI);
	blockGridHor2->x = iDivUp(num_volx, NUM_HEBRASX_ARI);
	blockGridHor2->y = iDivUp(num_aristas_hor2/num_volx, NUM_HEBRASY_ARI);

	blockGridEst->x = iDivUp(num_volx, NUM_HEBRASX_EST);
	blockGridEst->y = iDivUp(num_voly, NUM_HEBRASY_EST);
}

/*********************/
/* Okada deformation */
/*********************/

void comprobarYAplicarDeformaciones(double2 *d_datosVolumenesNivel_1, double *d_eta1InicialNivel, double tiempoActSubmalla,
		int numFaults, int *fallaOkada, double *defTime, int crop_flag, double crop_value, double *d_deltaTVolumenesNivel,
		float *d_vec, double *LON_C, double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE,
		double *DIP, double *RAKE, double *SLIP, int numVolxNivel0, int numVolyNivel0, double lon_ini, double incx,
		double lat_ini, double incy, dim3 blockGridEstNivel, dim3 threadBlockEst, double H, double T)
{
	bool encontrado;

	if (*fallaOkada < numFaults) {
		encontrado = true;
		while ((*fallaOkada < numFaults) && encontrado) {
			if (tiempoActSubmalla*T >= defTime[*fallaOkada]) {
				fprintf(stdout, "Applying Okada with fault %d\n", (*fallaOkada)+1);
				aplicarOkada(d_datosVolumenesNivel_1, d_eta1InicialNivel, crop_flag, crop_value, d_deltaTVolumenesNivel,
					d_vec, numVolxNivel0, numVolyNivel0, lon_ini, incx, lat_ini, incy, LON_C[*fallaOkada], LAT_C[*fallaOkada],
					DEPTH_C[*fallaOkada], FAULT_L[*fallaOkada], FAULT_W[*fallaOkada], STRIKE[*fallaOkada], DIP[*fallaOkada],
					RAKE[*fallaOkada], SLIP[*fallaOkada], blockGridEstNivel, threadBlockEst, H);
				(*fallaOkada)++;
			}
			else encontrado = false;
		}
	}
}

/*******************/
/* Free GPU memory */
/*******************/

void liberarMemoria(int numNiveles, double2 *d_datosVolumenesNivel_1, double2 *d_datosVolumenesNivel_2,
			tipoDatosSubmalla d_datosNivel, double *d_eta1MaximaNivel, double *d_tiemposLlegadaNivel,
			double *d_eta1InicialNivel, double *d_deltaTVolumenesNivel, double2 *d_acumuladorNivel1_1,
			double2 *d_acumuladorNivel1_2, int leer_fichero_puntos, int *d_posicionesVolumenesGuardado,
			double2 *d_datosVolumenesGuardado_1, double2 *d_datosVolumenesGuardado_2, float *d_vec)
{
	cudaFree(d_datosNivel.areaYCosPhi);
	cudaFree(d_datosNivel.anchoVolumenes);
	cudaFree(d_datosNivel.altoVolumenes);
	cudaFree(d_datosVolumenesNivel_1);
	cudaFree(d_datosVolumenesNivel_2);
	cudaFree(d_eta1MaximaNivel);
	cudaFree(d_tiemposLlegadaNivel);
	cudaFree(d_eta1InicialNivel);
	cudaFree(d_deltaTVolumenesNivel);
	cudaFree(d_acumuladorNivel1_1);
	cudaFree(d_acumuladorNivel1_2);
	cudaFree(d_vec);
	if (leer_fichero_puntos == 1) {
		cudaFree(d_posicionesVolumenesGuardado);
		cudaFree(d_datosVolumenesGuardado_1);
		cudaFree(d_datosVolumenesGuardado_2);
	}
}

/*****************/
/* Main function */
/*****************/

extern "C" int shallowWater(int numNiveles, int okada_flag, int numFaults, int crop_flag, double crop_value, double *LON_C,
			double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE, double *DIP, double *RAKE,
			double *SLIP, double *defTime, double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2, tipoDatosSubmalla datosNivel,
			int leer_fichero_puntos, int numPuntosGuardar, int *posicionesVolumenesGuardado, double *lonPuntos, double *latPuntos,
			int numVolxNivel0, int numVolyNivel0, int64_t numVolumenesNivel, double Hmin, char *nombre_bati, string prefijo,
			double borde_sup, double borde_inf, double borde_izq, double borde_der, int tam_spongeSup, int tam_spongeInf,
			int tam_spongeIzq, int tam_spongeDer, double tiempo_tot, double tiempoGuardarNetCDF, double tiempoGuardarSeries,
			double CFL, double mf0, double vmax, double epsilon_h, double hpos, double cvis, double dif_at, double L, double H,
			double Q, double T, char *version, double *tiempo)
{
	double2 *d_datosVolumenesNivel_1;
	double2 *d_datosVolumenesNivel_2;
	double2 *d_acumuladorNivel_1, *d_acumuladorNivel_2;
	tipoDatosSubmalla d_datosNivel;
	double *d_deltaTVolumenesNivel;
	double *d_eta1MaximaNivel;
	double *d_tiemposLlegadaNivel;
	double *d_eta1InicialNivel;
	float *vec, *d_vec;
	cudaError_t err;
	int num = 0;
	int fallaOkada = 0;
	bool saltar_deformacion = false;
	double sea_level = SEA_LEVEL/H;
	dim3 blockGridVer1Nivel;
	dim3 blockGridVer2Nivel;
	dim3 blockGridHor1Nivel;
	dim3 blockGridHor2Nivel;
	dim3 blockGridEstNivel;
	dim3 threadBlockAri(NUM_HEBRASX_ARI, NUM_HEBRASY_ARI);
	dim3 threadBlockEst(NUM_HEBRASX_EST, NUM_HEBRASY_EST);
	dim3 blockGridPuntos(iDivUp(numPuntosGuardar, NUM_HEBRAS_PUNTOS), 1);
	dim3 threadBlockPuntos(NUM_HEBRAS_PUNTOS, 1);
	double lon_ini = datosNivel.longitud[0];
	double lat_ini = datosNivel.latitud[0];
	double incx = (datosNivel.longitud[numVolxNivel0-1] - datosNivel.longitud[0])/(numVolxNivel0-1);
	double incy = (datosNivel.latitud[numVolyNivel0-1] - datosNivel.latitud[0])/(numVolyNivel0-1);

	int *d_posicionesVolumenesGuardado;
	double2 *d_datosVolumenesGuardado_1, *d_datosVolumenesGuardado_2;
	float *etaPuntosGuardado;
	float *uPuntosGuardado;
	float *vPuntosGuardado;
	float *etaMinPuntosGuardado, *etaMaxPuntosGuardado;
	int num_ts = 0;

	int64_t tam_datosVolDoubleNivel;
	int64_t tam_datosVolDouble2Nivel;
	int64_t tam_datosVolFloatNivel;
	int64_t tam_datosAcumDouble2Nivel;
	int64_t tam_datosVolGuardadoDouble2 = ((int64_t) numPuntosGuardar)*sizeof(double2);
	double deltaTNivel;
	double tiempoActSubmalla;
	double sigTiempoGuardarNetCDF = 0.0;
	double sigTiempoGuardarSeries = 0.0;
	StopWatchInterface *timer = NULL;
	int i, iter;

	cudaSetDevice(0);

	tiempoActSubmalla = 0.0;
	obtenerTamBloquesKernel(numVolxNivel0, numVolyNivel0, &blockGridVer1Nivel, &blockGridVer2Nivel,
		&blockGridHor1Nivel, &blockGridHor2Nivel, &blockGridEstNivel);

	cudaMalloc((void **)&(d_datosNivel.areaYCosPhi), numVolyNivel0*sizeof(double2));
	cudaMalloc((void **)&(d_datosNivel.anchoVolumenes), (numVolyNivel0+1)*sizeof(double));
	cudaMalloc((void **)&(d_datosNivel.altoVolumenes), numVolyNivel0*sizeof(double));
	tam_datosVolDoubleNivel = numVolumenesNivel*sizeof(double);
	tam_datosVolDouble2Nivel = numVolumenesNivel*sizeof(double2);
	tam_datosVolFloatNivel = numVolumenesNivel*sizeof(float);
	tam_datosAcumDouble2Nivel = tam_datosVolDouble2Nivel;
	cudaMalloc((void **)&d_datosVolumenesNivel_1, tam_datosVolDouble2Nivel);
	cudaMalloc((void **)&d_datosVolumenesNivel_2, tam_datosVolDouble2Nivel);
	cudaMalloc((void **)&d_eta1MaximaNivel, tam_datosVolDoubleNivel);
	cudaMalloc((void **)&d_tiemposLlegadaNivel, tam_datosVolDoubleNivel);
	cudaMalloc((void **)&d_eta1InicialNivel, tam_datosVolDoubleNivel);
	cudaMalloc((void **)&d_vec, tam_datosVolFloatNivel);
	cudaMalloc((void **)&d_deltaTVolumenesNivel, tam_datosVolDoubleNivel);
	cudaMalloc((void **)&d_acumuladorNivel_1, tam_datosAcumDouble2Nivel);
	err = cudaMalloc( (void **)&d_acumuladorNivel_2, tam_datosAcumDouble2Nivel);
	if (err == cudaErrorMemoryAllocation) {
		cudaFree(d_datosNivel.areaYCosPhi);
		cudaFree(d_datosNivel.anchoVolumenes);
		cudaFree(d_datosNivel.altoVolumenes);
		cudaFree(d_datosVolumenesNivel_1);
		cudaFree(d_datosVolumenesNivel_2);
		cudaFree(d_eta1MaximaNivel);
		cudaFree(d_tiemposLlegadaNivel);
		cudaFree(d_eta1InicialNivel);
		cudaFree(d_vec);
		cudaFree(d_deltaTVolumenesNivel);
		cudaFree(d_acumuladorNivel_1);
		return 1;
	}
	if (leer_fichero_puntos == 1) {
		cudaMalloc((void **)&d_posicionesVolumenesGuardado, ((int64_t) numPuntosGuardar)*sizeof(int));
		cudaMalloc((void **)&d_datosVolumenesGuardado_1, tam_datosVolGuardadoDouble2);
		err = cudaMalloc((void **)&d_datosVolumenesGuardado_2, tam_datosVolGuardadoDouble2);
		if (err == cudaErrorMemoryAllocation) {
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel, d_eta1MaximaNivel,
				d_tiemposLlegadaNivel, d_eta1InicialNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2,
				1, d_posicionesVolumenesGuardado, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2, d_vec);
			return 1;
		}
	}

	cudaFuncSetCacheConfig(procesarAristasVerDeltaTInicialNivel0GPU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(procesarAristasHorDeltaTInicialNivel0GPU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(procesarAristasVerNivel0Paso1GPU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(procesarAristasVerNivel0Paso2GPU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(procesarAristasHorNivel0Paso1GPU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(procesarAristasHorNivel0Paso2GPU, cudaFuncCachePreferL1);

	cudaMemcpy(d_datosNivel.areaYCosPhi, datosNivel.areaYCosPhi, numVolyNivel0*sizeof(double2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_datosNivel.anchoVolumenes, datosNivel.anchoVolumenes, (numVolyNivel0+1)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_datosNivel.altoVolumenes, datosNivel.altoVolumenes, numVolyNivel0*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_datosVolumenesNivel_1, datosVolumenesNivel_1, tam_datosVolDouble2Nivel, cudaMemcpyHostToDevice);
	cudaMemcpy(d_datosVolumenesNivel_2, datosVolumenesNivel_2, tam_datosVolDouble2Nivel, cudaMemcpyHostToDevice);
	if (leer_fichero_puntos == 1) {
		cudaMemcpy(d_posicionesVolumenesGuardado, posicionesVolumenesGuardado, ((int64_t) numPuntosGuardar)*sizeof(int), cudaMemcpyHostToDevice);
	}

	inicializarEta1MaximaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
		d_eta1MaximaNivel, numVolxNivel0, numVolyNivel0, epsilon_h);
	inicializarEta1MaximaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
		d_eta1InicialNivel, numVolxNivel0, numVolyNivel0, epsilon_h);
	inicializarTiemposLlegadaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_tiemposLlegadaNivel, numVolxNivel0, numVolyNivel0);

	// INICIO NETCDF
	double mf0_ini = sqrt(mf0*pow(H,4.0/3.0)/(9.81*L));
	if (tiempoGuardarNetCDF >= 0.0) {
		vec = (float *) malloc(tam_datosVolFloatNivel);
		if (vec == NULL) {
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel, d_eta1MaximaNivel,
				d_tiemposLlegadaNivel, d_eta1InicialNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2,
				leer_fichero_puntos, d_posicionesVolumenesGuardado, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2, d_vec);
			return 2;
		}
		for (i=0; i<numVolumenesNivel; i++)
			vec[i] = (float) ((datosVolumenesNivel_1[i].y + Hmin)*H);
		iter = crearFicherosNC(nombre_bati, okada_flag, (char *) prefijo.c_str(), numVolxNivel0, numVolyNivel0,
					datosNivel.longitud, datosNivel.latitud, tiempo_tot*T, CFL, epsilon_h*H, mf0_ini, vmax*Q/H,
					hpos*H, 1.0-cvis, dif_at*H, borde_sup, borde_inf, borde_izq, borde_der, numFaults, defTime,
					LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP, crop_flag, crop_value,
					H, vec, version);
		if (iter == 1) {
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel, d_eta1MaximaNivel,
				d_tiemposLlegadaNivel, d_eta1InicialNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2,
				leer_fichero_puntos, d_posicionesVolumenesGuardado, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2, d_vec);
			free(vec);
			return 2;
		}
	}
	if (leer_fichero_puntos == 1) {
		etaPuntosGuardado = (float *) malloc(numPuntosGuardar*sizeof(float));
		uPuntosGuardado = (float *) malloc(numPuntosGuardar*sizeof(float));
		vPuntosGuardado = (float *) malloc(numPuntosGuardar*sizeof(float));
		etaMinPuntosGuardado = (float *) malloc(numPuntosGuardar*sizeof(float));
		etaMaxPuntosGuardado = (float *) malloc(numPuntosGuardar*sizeof(float));
		if (etaMaxPuntosGuardado == NULL) {
			if (etaPuntosGuardado != NULL)			free(etaPuntosGuardado);
			if (uPuntosGuardado != NULL)			free(uPuntosGuardado);
			if (vPuntosGuardado != NULL)			free(vPuntosGuardado);
			if (etaMinPuntosGuardado != NULL)		free(etaMinPuntosGuardado);
			if (tiempoGuardarNetCDF >= 0.0)
				free(vec);
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel, d_eta1MaximaNivel,
				d_tiemposLlegadaNivel, d_eta1InicialNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2,
				leer_fichero_puntos, d_posicionesVolumenesGuardado, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2, d_vec);
			return 2;
		}
		initTimeSeriesNC(nombre_bati, (char *) prefijo.c_str(), numPuntosGuardar, lonPuntos, latPuntos, tiempo_tot*T,
			CFL, epsilon_h*H, mf0_ini, vmax*Q/H, hpos*H, 1.0-cvis, dif_at*H, borde_sup, borde_inf, borde_izq, borde_der,
			numFaults, defTime, LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP, okada_flag, crop_flag,
			crop_value, H, version);
		for (i=0; i<numPuntosGuardar; i++) {
			if (posicionesVolumenesGuardado[i] == -1) {
				etaPuntosGuardado[i] = -9999.0f;
				uPuntosGuardado[i] = -9999.0f;
				vPuntosGuardado[i] = -9999.0f;
				etaMinPuntosGuardado[i] = -9999.0f;
				etaMaxPuntosGuardado[i] = -9999.0f;
			}
			else {
				etaMinPuntosGuardado[i] = 1e30f;
				etaMaxPuntosGuardado[i] = -1e30f;
			}
		}
	}
	// FIN NETCDF

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel);
	cudaMemset(d_acumuladorNivel_2, 0, tam_datosVolDouble2Nivel);

	sdkCreateTimer(&timer);
	sdkStartTimer(&timer);
	deltaTNivel = obtenerDeltaTInicialNivel0(d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, &d_datosNivel,
					d_deltaTVolumenesNivel, d_acumuladorNivel_1, numVolxNivel0, numVolyNivel0, borde_izq,
					borde_der, borde_sup, borde_inf, CFL, epsilon_h, blockGridVer1Nivel, blockGridVer2Nivel,
					blockGridHor1Nivel, blockGridHor2Nivel, threadBlockAri, blockGridEstNivel, threadBlockEst);
	if (numFaults > 0) {
		if (defTime[0] < EPSILON) {
			truncarDeltaTNivel0ParaSigDeformacionSaltando(okada_flag, numFaults, fallaOkada, tiempoActSubmalla,
				defTime, T, &deltaTNivel, &saltar_deformacion);
		}
		else {
			truncarDeltaTNivel0ParaSigDeformacion(okada_flag, numFaults, fallaOkada, tiempoActSubmalla,
				defTime, T, &deltaTNivel, &saltar_deformacion);
		}
	}
	fprintf(stdout, "Initial deltaT = %e sec\n", deltaTNivel*T);

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel);

	iter = 1;
	while (tiempoActSubmalla < tiempo_tot) {
		comprobarYAplicarDeformaciones(d_datosVolumenesNivel_1, d_eta1InicialNivel, tiempoActSubmalla, numFaults,
			&fallaOkada, defTime, crop_flag, crop_value, d_deltaTVolumenesNivel, d_vec, LON_C, LAT_C, DEPTH_C,
			FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP, numVolxNivel0, numVolyNivel0, lon_ini, incx, lat_ini,
			incy, blockGridEstNivel, threadBlockEst, H, T);

		// INICIO NETCDF
		if ((tiempoGuardarNetCDF >= 0.0) && (tiempoActSubmalla >= sigTiempoGuardarNetCDF)) {
			sigTiempoGuardarNetCDF += tiempoGuardarNetCDF;
			guardarNetCDFNivel0(d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_vec, vec, num, numVolxNivel0,
				numVolyNivel0, Hmin, tiempoActSubmalla, epsilon_h, blockGridEstNivel, threadBlockEst, H, Q, T);
			num++;
		}
		if ((tiempoGuardarSeries >= 0.0) && (tiempoActSubmalla >= sigTiempoGuardarSeries)) {
			sigTiempoGuardarSeries += tiempoGuardarSeries;
			escribirVolumenesGuardadoGPU<<<blockGridPuntos, threadBlockPuntos>>>(d_datosVolumenesNivel_1,
				d_datosVolumenesNivel_2, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2,
				d_posicionesVolumenesGuardado, numPuntosGuardar, epsilon_h);
			cudaMemcpy(datosVolumenesNivel_1, d_datosVolumenesGuardado_1, tam_datosVolGuardadoDouble2, cudaMemcpyDeviceToHost);
			cudaMemcpy(datosVolumenesNivel_2, d_datosVolumenesGuardado_2, tam_datosVolGuardadoDouble2, cudaMemcpyDeviceToHost);
			guardarSerieTiemposNivel0(datosVolumenesNivel_1, datosVolumenesNivel_2, numPuntosGuardar, posicionesVolumenesGuardado,
				etaPuntosGuardado, uPuntosGuardado, vPuntosGuardado, etaMinPuntosGuardado, etaMaxPuntosGuardado,
				Hmin, num_ts, tiempoActSubmalla, H, Q, T);
			num_ts++;
		}
		// FIN NETCDF
		actualizarEta1MaximaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1, d_eta1MaximaNivel,
			numVolxNivel0, numVolyNivel0, tiempoActSubmalla, epsilon_h);
		actualizarTiemposLlegadaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1, d_tiemposLlegadaNivel,
			d_eta1InicialNivel, numVolxNivel0, numVolyNivel0, tiempoActSubmalla, dif_at);

		deltaTNivel = siguientePasoNivel0(d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, &d_datosNivel, d_acumuladorNivel_1,
						d_acumuladorNivel_2, numVolxNivel0, numVolyNivel0, d_deltaTVolumenesNivel, borde_sup, borde_inf,
						borde_izq, borde_der, Hmin, tam_spongeSup, tam_spongeInf, tam_spongeIzq, tam_spongeDer, sea_level,
						&tiempoActSubmalla, CFL, deltaTNivel, mf0, vmax, epsilon_h, hpos, cvis, L, H, tam_datosVolDouble2Nivel,
						blockGridVer1Nivel, blockGridVer2Nivel, blockGridHor1Nivel, blockGridHor2Nivel, threadBlockAri,
						blockGridEstNivel, threadBlockEst);

		if (! saltar_deformacion) {
			truncarDeltaTNivel0ParaSigDeformacion(okada_flag, numFaults, fallaOkada, tiempoActSubmalla, defTime,
				T, &deltaTNivel, &saltar_deformacion);
		}
		else {
			truncarDeltaTNivel0ParaSigDeformacionSaltando(okada_flag, numFaults, fallaOkada, tiempoActSubmalla,
				defTime, T, &deltaTNivel, &saltar_deformacion);
		}

		cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel);
		cudaMemset(d_acumuladorNivel_2, 0, tam_datosVolDouble2Nivel);

		fprintf(stdout, "Iteration %3d, deltaT = %e sec, ", iter, deltaTNivel*T);
		fprintf(stdout, "Time = %g sec\n", tiempoActSubmalla*T);
		iter++;
	}
	cudaDeviceSynchronize();
	sdkStopTimer(&timer);
	*tiempo = sdkGetTimerValue(&timer)*0.001;

	// INICIO NETCDF
	if (tiempoGuardarNetCDF >= 0.0) {
		if (okada_flag == OKADA_STANDARD) {
			obtenerBatimetriaModificadaNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1, d_vec,
				numVolxNivel0, numVolyNivel0, Hmin, H);
			cudaMemcpy(vec, d_vec, tam_datosVolFloatNivel, cudaMemcpyDeviceToHost);
			guardarBatimetriaModificadaNC(vec);
		}

		obtenerEta1MaximaNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_eta1MaximaNivel, d_vec,
			numVolxNivel0, numVolyNivel0, Hmin, H);
		cudaMemcpy(vec, d_vec, tam_datosVolFloatNivel, cudaMemcpyDeviceToHost);
		guardarEta1MaximaNC(vec);

		obtenerTiemposLlegadaNetCDF<<<blockGridEstNivel, threadBlockEst>>>(d_tiemposLlegadaNivel, d_vec,
			numVolxNivel0, numVolyNivel0, T);
		cudaMemcpy(vec, d_vec, tam_datosVolFloatNivel, cudaMemcpyDeviceToHost);
		guardarTiemposLlegadaNC(vec);

		cerrarFicheroNC();
		free(vec);
	}
	if (leer_fichero_puntos == 1) {
		escribirVolumenesGuardadoGPU<<<blockGridPuntos, threadBlockPuntos>>>(d_datosVolumenesNivel_1,
			d_datosVolumenesNivel_2, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2,
			d_posicionesVolumenesGuardado, numPuntosGuardar, epsilon_h);
		cudaMemcpy(datosVolumenesNivel_1, d_datosVolumenesGuardado_1, tam_datosVolGuardadoDouble2, cudaMemcpyDeviceToHost);
		obtenerBatimetriaParaSerieTiempos(datosVolumenesNivel_1, numPuntosGuardar, posicionesVolumenesGuardado,
			etaPuntosGuardado, Hmin, H);
		guardarBatimetriaModificadaTimeSeriesNC(etaPuntosGuardado);
		guardarAmplitudesTimeSeriesNC(etaMinPuntosGuardado, etaMaxPuntosGuardado);
		closeTimeSeriesNC();
		free(etaPuntosGuardado);
		free(uPuntosGuardado);
		free(vPuntosGuardado);
		free(etaMinPuntosGuardado);
		free(etaMaxPuntosGuardado);
	}
	// FIN NETCDF

	sdkDeleteTimer(&timer);
	liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel, d_eta1MaximaNivel,
		d_tiemposLlegadaNivel, d_eta1InicialNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2,
		leer_fichero_puntos, d_posicionesVolumenesGuardado, d_datosVolumenesGuardado_1, d_datosVolumenesGuardado_2, d_vec);

	return 0;
}

