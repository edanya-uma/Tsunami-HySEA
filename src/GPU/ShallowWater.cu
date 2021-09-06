#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <fstream>
#include "helper_timer.h"
#include "Reduccion_kernel.cu"
#include "Volumen_kernel.cu"
#include "Arista_kernel.cu"
#include "netcdf.cu"

void guardarNetCDFNivel0(double2 *datosVolumenes_1, double2 *datosVolumenes_2, double *vec, int num, int num_volx,
				int num_voly, double Hmin, double tiempo_act, double epsilon_h, double H, double Q, double T)
{
	double tiempo = tiempo_act*T;
	double2 datos;
	int i;
	int num_volumenes = num_volx*num_voly;

	escribirTiempoNC(num, tiempo);
	for (i=0; i<num_volumenes; i++) {
		datos = datosVolumenes_1[i];
		vec[i] = ( (datos.x < epsilon_h) ? -9999.0 : (datos.x - datos.y - Hmin)*H );
	}
	escribirEtaNC(num_volx, num_voly, num, vec);
	for (i=0; i<num_volumenes; i++) {
		vec[i] = datosVolumenes_2[i].x*Q;
	}
	escribirQxNC(num_volx, num_voly, num, vec);
	for (i=0; i<num_volumenes; i++) {
		vec[i] = datosVolumenes_2[i].y*Q;
	}
	escribirQyNC(num_volx, num_voly, num, vec);
}

void obtenerTamBloquesKernel(int num_volx, int num_voly, dim3 *blockGridVer1, dim3 *blockGridVer2,
							dim3 *blockGridHor1, dim3 *blockGridHor2, dim3 *blockGridEst)
{
	int num_aristas_ver1, num_aristas_ver2;
	int num_aristas_hor1, num_aristas_hor2;

	num_aristas_ver1 = (num_volx/2 + 1)*num_voly;
	num_aristas_ver2 = (num_volx&1 == 0) ? num_volx*num_voly/2 : num_aristas_ver1;
	num_aristas_hor1 = (num_voly/2 + 1)*num_volx;
	num_aristas_hor2 = (num_voly&1 == 0) ? num_volx*num_voly/2 : num_aristas_hor1;

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

double obtenerDeltaTInicialNivel0(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2, tipoDatosSubmalla *d_datosNivel0,
		double *d_deltaTVolumenesNivel0, double2 *d_acumulador_1, int numVolxNivel0, int numVolyNivel0, double borde_izq,
		double borde_der, double borde_sup, double borde_inf, double CFL, double epsilon_h, dim3 blockGridVer1, dim3 blockGridVer2,
		dim3 blockGridHor1, dim3 blockGridHor2, dim3 threadBlockAri, dim3 blockGridEst, dim3 threadBlockEst)
{
	double delta_T;

	procesarAristasVerDeltaTInicialNivel0GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, d_datosNivel0->altoVolumenes,
		d_acumulador_1, epsilon_h, 1);
	procesarAristasVerDeltaTInicialNivel0GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, d_datosNivel0->altoVolumenes,
		d_acumulador_1, epsilon_h, 2);
	procesarAristasHorDeltaTInicialNivel0GPU<<<blockGridHor1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		numVolxNivel0, numVolyNivel0, borde_sup, borde_inf, d_datosNivel0->anchoVolumenes, d_acumulador_1, epsilon_h, 1);
	procesarAristasHorDeltaTInicialNivel0GPU<<<blockGridHor2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		numVolxNivel0, numVolyNivel0, borde_sup, borde_inf, d_datosNivel0->anchoVolumenes, d_acumulador_1, epsilon_h, 2);
	obtenerDeltaTVolumenesGPU<<<blockGridEst, threadBlockEst>>>(d_datosNivel0->areaYCosPhi, d_acumulador_1,
			d_deltaTVolumenesNivel0, numVolxNivel0, numVolyNivel0, CFL);
	delta_T = obtenerMinimoReduccion<double>(d_deltaTVolumenesNivel0, numVolxNivel0*numVolyNivel0);

	return delta_T;
}

double siguientePasoNivel0(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2, tipoDatosSubmalla *d_datosNivel0,
		double2 *d_acumuladorNivel_1, double2 *d_acumuladorNivel_2, int numVolxNivel0, int numVolyNivel0,
		double *d_deltaTVolumenesNivel0, double borde_sup, double borde_inf, double borde_izq, double borde_der,
		double *tiempo_act, double CFL, double delta_T, double mf0, double vmax, double epsilon_h, double hpos,
		double cvis, double L, double H, int tam_datosVolDouble2Nivel0, dim3 blockGridVer1, dim3 blockGridVer2,
		dim3 blockGridHor1, dim3 blockGridHor2, dim3 threadBlockAri, dim3 blockGridEst, dim3 threadBlockEst)
{
	// PASO 1
	procesarAristasVerNivel0Paso1GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, delta_T,
		d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasVerNivel0Paso1GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, delta_T,
		d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 2);
	procesarAristasHorNivel0Paso1GPU<<<blockGridHor1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf, delta_T,
		d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasHorNivel0Paso1GPU<<<blockGridHor2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf, delta_T,
		d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 2);
	obtenerEstadosPaso1Nivel0GPU<<<blockGridEst, threadBlockEst>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_acumuladorNivel_1, numVolxNivel0, numVolyNivel0, delta_T);

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel0);

	// PASO 2
	procesarAristasVerNivel0Paso2GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasVerNivel0Paso2GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 2);
	procesarAristasHorNivel0Paso2GPU<<<blockGridHor1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasHorNivel0Paso2GPU<<<blockGridHor2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 2);
	obtenerEstadoYDeltaTVolumenesNivel0GPU<<<blockGridEst, threadBlockEst>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_acumuladorNivel_1, d_acumuladorNivel_2, d_deltaTVolumenesNivel0,
		numVolxNivel0, numVolyNivel0, CFL, delta_T, mf0, vmax, hpos, epsilon_h);

	*tiempo_act += delta_T;

	delta_T = obtenerMinimoReduccion<double>(d_deltaTVolumenesNivel0, numVolxNivel0*numVolyNivel0);

	return delta_T;
}

void liberarMemoria(int numNiveles, double2 *d_datosVolumenesNivel_1, double2 *d_datosVolumenesNivel_2,
			tipoDatosSubmalla d_datosNivel, double *d_eta1MaximaNivel, double *d_deltaTVolumenesNivel,
			double2 *d_acumuladorNivel1_1, double2 *d_acumuladorNivel1_2)
{
	cudaFree(d_datosNivel.areaYCosPhi);
	cudaFree(d_datosNivel.anchoVolumenes);
	cudaFree(d_datosNivel.altoVolumenes);
	cudaFree(d_datosVolumenesNivel_1);
	cudaFree(d_datosVolumenesNivel_2);
	cudaFree(d_eta1MaximaNivel);
	cudaFree(d_deltaTVolumenesNivel);
	cudaFree(d_acumuladorNivel1_1);
	cudaFree(d_acumuladorNivel1_2);
}

extern "C" int shallowWater(int numNiveles, double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2,
			tipoDatosSubmalla datosNivel, int leer_fichero_puntos, int numVolxNivel0, int numVolyNivel0,
			int numVolumenesNivel, double Hmin, char *nombre_bati, string prefijo, double borde_sup, double borde_inf,
			double borde_izq, double borde_der, double tiempo_tot, double tiempoGuardarNetCDF, double CFL, double mf0,
			double vmax, double epsilon_h, double hpos, double cvis, double L, double H, double Q, double T, double *tiempo)
{
	double2 *d_datosVolumenesNivel_1;
	double2 *d_datosVolumenesNivel_2;
	double2 *d_acumuladorNivel_1, *d_acumuladorNivel_2;
	tipoDatosSubmalla d_datosNivel;
	double *d_deltaTVolumenesNivel;
	double *d_eta1MaximaNivel;
	double *vec, *p_eta;
	cudaError_t err;
	int num = 0;
	dim3 blockGridVer1Nivel;
	dim3 blockGridVer2Nivel;
	dim3 blockGridHor1Nivel;
	dim3 blockGridHor2Nivel;
	dim3 blockGridEstNivel;
	dim3 blockGridFanNivel;
	dim3 threadBlockAri(NUM_HEBRASX_ARI, NUM_HEBRASY_ARI);
	dim3 threadBlockEst(NUM_HEBRASX_EST, NUM_HEBRASY_EST);

	int tam_datosVolDoubleNivel;
	int tam_datosVolDouble2Nivel;
	int tam_datosDeltaTDouble;
	int tam_datosAcumDouble2Nivel;
	double deltaTNivel;
	double tiempoActSubmalla;
	double sigTiempoGuardarNetCDF = 0.0;
	StopWatchInterface *timer = NULL;
	int i, iter;

	cudaSetDevice(0);

	tiempoActSubmalla = 0.0;
	obtenerTamBloquesKernel(numVolxNivel0, numVolyNivel0, &blockGridVer1Nivel, &blockGridVer2Nivel,
		&blockGridHor1Nivel, &blockGridHor2Nivel, &blockGridEstNivel);

	tam_datosDeltaTDouble = 0;
	tam_datosAcumDouble2Nivel = 0;
	cudaMalloc((void **)&(d_datosNivel.areaYCosPhi), numVolyNivel0*sizeof(double2));
	cudaMalloc((void **)&(d_datosNivel.anchoVolumenes), (numVolyNivel0+1)*sizeof(double));
	cudaMalloc((void **)&(d_datosNivel.altoVolumenes), numVolyNivel0*sizeof(double));
	tam_datosVolDoubleNivel = numVolumenesNivel*sizeof(double);
	tam_datosVolDouble2Nivel = numVolumenesNivel*sizeof(double2);
	tam_datosDeltaTDouble = max(tam_datosDeltaTDouble, tam_datosVolDoubleNivel);
	tam_datosAcumDouble2Nivel = max(tam_datosAcumDouble2Nivel, tam_datosVolDouble2Nivel);
	cudaMalloc((void **)&d_datosVolumenesNivel_1, tam_datosVolDouble2Nivel);
	cudaMalloc((void **)&d_datosVolumenesNivel_2, tam_datosVolDouble2Nivel);
	cudaMalloc((void **)&d_eta1MaximaNivel, tam_datosVolDoubleNivel);
	cudaMalloc((void **)&d_deltaTVolumenesNivel, tam_datosDeltaTDouble);
	cudaMalloc((void **)&d_acumuladorNivel_1, tam_datosAcumDouble2Nivel);
	err = cudaMalloc( (void **)&d_acumuladorNivel_2, tam_datosAcumDouble2Nivel);
	if (err == cudaErrorMemoryAllocation) {
		cudaFree(d_datosNivel.areaYCosPhi);
		cudaFree(d_datosNivel.anchoVolumenes);
		cudaFree(d_datosNivel.altoVolumenes);
		cudaFree(d_datosVolumenesNivel_1);
		cudaFree(d_datosVolumenesNivel_2);
		cudaFree(d_eta1MaximaNivel);
		cudaFree(d_deltaTVolumenesNivel);
		cudaFree(d_acumuladorNivel_1);
		return 1;
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
	cudaDeviceSynchronize();

	inicializarEta1MaximaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
		d_eta1MaximaNivel, numVolxNivel0, numVolyNivel0, epsilon_h);

	// INICIO NETCDF
	if (tiempoGuardarNetCDF >= 0.0) {
		vec = (double *) malloc(tam_datosDeltaTDouble);
		if (vec == NULL) {
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel,
				d_eta1MaximaNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2);
			return 2;
		}
		double mf0_ini = sqrt(mf0*pow(H,4.0/3.0)/(9.81*L));
		for (i=0; i<numVolumenesNivel; i++)
			vec[i] = (datosVolumenesNivel_1[i].y + Hmin)*H;
		iter = crearFicherosNC(nombre_bati, (char *) prefijo.c_str(), numVolxNivel0, numVolyNivel0, datosNivel.longitud,
					datosNivel.latitud, tiempo_tot*T, CFL, epsilon_h*H, mf0_ini, vmax*Q/H, hpos*H, 1.0-cvis, vec);
		if (iter == 1) {
			liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel,
				d_eta1MaximaNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2);
			free(vec);
			return 2;
		}
	}
	// FIN NETCDF

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel);
	cudaMemset(d_acumuladorNivel_2, 0, tam_datosVolDouble2Nivel);

	deltaTNivel = obtenerDeltaTInicialNivel0(d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, &d_datosNivel,
					d_deltaTVolumenesNivel, d_acumuladorNivel_1, numVolxNivel0, numVolyNivel0, borde_izq,
					borde_der, borde_sup, borde_inf, CFL, epsilon_h, blockGridVer1Nivel, blockGridVer2Nivel,
					blockGridHor1Nivel, blockGridHor2Nivel, threadBlockAri, blockGridEstNivel, threadBlockEst);
	fprintf(stdout, "Initial deltaT = %e sec\n", deltaTNivel*T);

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel);

	iter = 1;
	sdkCreateTimer(&timer);
	sdkStartTimer(&timer);
	while (tiempoActSubmalla < tiempo_tot) {
		// INICIO NETCDF
		if ((tiempoGuardarNetCDF >= 0.0) && (tiempoActSubmalla >= sigTiempoGuardarNetCDF)) {
			cudaMemcpy(datosVolumenesNivel_1, d_datosVolumenesNivel_1, tam_datosVolDouble2Nivel, cudaMemcpyDeviceToHost);
			cudaMemcpy(datosVolumenesNivel_2, d_datosVolumenesNivel_2, tam_datosVolDouble2Nivel, cudaMemcpyDeviceToHost);
			guardarNetCDFNivel0(datosVolumenesNivel_1, datosVolumenesNivel_2, vec, num, numVolxNivel0, numVolyNivel0,
				Hmin, tiempoActSubmalla, epsilon_h, H, Q, T);
			num++;
			sigTiempoGuardarNetCDF += tiempoGuardarNetCDF;
		}
		// FIN NETCDF
		actualizarEta1MaximaNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
			d_eta1MaximaNivel, numVolxNivel0, numVolyNivel0, tiempoActSubmalla, epsilon_h);

		deltaTNivel = siguientePasoNivel0(d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, &d_datosNivel,
						d_acumuladorNivel_1, d_acumuladorNivel_2, numVolxNivel0, numVolyNivel0, d_deltaTVolumenesNivel,
						borde_sup, borde_inf, borde_izq, borde_der, &tiempoActSubmalla, CFL, deltaTNivel,
						mf0, vmax, epsilon_h, hpos, cvis, L, H, tam_datosVolDouble2Nivel, blockGridVer1Nivel,
						blockGridVer2Nivel, blockGridHor1Nivel, blockGridHor2Nivel, threadBlockAri,
						blockGridEstNivel, threadBlockEst);

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
		cudaMemcpy(datosVolumenesNivel_1, d_eta1MaximaNivel, tam_datosVolDoubleNivel, cudaMemcpyDeviceToHost);
		p_eta = (double *) datosVolumenesNivel_1;
		for (i=0; i<numVolumenesNivel; i++) {
			if (p_eta[i] < -1e20)
				vec[i] = -9999.0;
			else
				vec[i] = max((p_eta[i] - Hmin)*H, 0.0);
		}
		cerrarFicheroNC(vec);
		free(vec);
	}
	// FIN NETCDF

	sdkDeleteTimer(&timer);
	liberarMemoria(numNiveles, d_datosVolumenesNivel_1, d_datosVolumenesNivel_2, d_datosNivel,
		d_eta1MaximaNivel, d_deltaTVolumenesNivel, d_acumuladorNivel_1, d_acumuladorNivel_2);

	return 0;
}

