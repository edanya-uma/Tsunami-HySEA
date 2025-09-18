#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include "Reduccion_kernel.cu"
#include "Volumen_kernel.cu"
#include "Arista_kernel.cu"

double obtenerDeltaTInicialNivel0(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2, tipoDatosSubmalla *d_datosNivel0,
		double *d_deltaTVolumenesNivel0, double2 *d_acumulador_1, int numVolxNivel0, int numVolyNivel0, double borde_izq,
		double borde_der, double borde_sup, double borde_inf, bool es_periodica, double CFL, double epsilon_h, dim3 blockGridVer1,
		dim3 blockGridVer2, dim3 blockGridHor1, dim3 blockGridHor2, dim3 threadBlockAri, dim3 blockGridEst, dim3 threadBlockEst)
{
	double delta_T;

	procesarAristasVerDeltaTInicialNivel0GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, es_periodica, d_datosNivel0->altoVolumenes,
		d_acumulador_1, epsilon_h, 1);
	procesarAristasVerDeltaTInicialNivel0GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, numVolxNivel0, numVolyNivel0, borde_izq, borde_der, es_periodica, d_datosNivel0->altoVolumenes,
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

void truncarDeltaTNivel0ParaSigDeformacion(int okada_flag, int numFaults, int fallaOkada, double tiempoActSubmalla,
		double *defTime, double T, double *delta_T, bool *saltar_deformacion)
{
	if (fallaOkada < numFaults) {
		if ((tiempoActSubmalla + (*delta_T))*T > defTime[fallaOkada]) {
			*delta_T = defTime[fallaOkada]/T - tiempoActSubmalla;
			*saltar_deformacion = true;
		}
	}
}

void truncarDeltaTNivel0ParaSigDeformacionSaltando(int okada_flag, int numFaults, int fallaOkada, double tiempoActSubmalla,
		double *defTime, double T, double *delta_T, bool *saltar_deformacion)
{
	int i;
	bool encontrado = false;

	*saltar_deformacion = false;
	if (fallaOkada < numFaults) {
		i = fallaOkada+1;
		while ((i < numFaults) && (! encontrado)) {
			if (fabs(defTime[i] - defTime[fallaOkada]) > EPSILON)
				encontrado = true;
			else
				i++;
		}
		if (encontrado) {
			if ((tiempoActSubmalla + (*delta_T))*T > defTime[i]) {
				*delta_T = defTime[i]/T - tiempoActSubmalla;
				*saltar_deformacion = true;
			}
		}
	}
}

double siguientePasoNivel0(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2, tipoDatosSubmalla *d_datosNivel0,
		double2 *d_acumuladorNivel_1, double2 *d_acumuladorNivel_2, int numVolxNivel0, int numVolyNivel0, double *d_deltaTVolumenesNivel0,
		double borde_sup, double borde_inf, double borde_izq, double borde_der, bool es_periodica, double Hmin, int tam_spongeSup,
		int tam_spongeInf, int tam_spongeIzq, int tam_spongeDer, double sea_level, double *tiempo_act, double CFL, double delta_T,
		double mf0, double vmax, double epsilon_h, double hpos, double cvis, double L, double H, int64_t tam_datosVolDouble2Nivel0,
		dim3 blockGridVer1, dim3 blockGridVer2, dim3 blockGridHor1, dim3 blockGridHor2, dim3 threadBlockAri, dim3 blockGridEst,
		dim3 threadBlockEst)
{
	// PASO 1
	procesarAristasVerNivel0Paso1GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		es_periodica, delta_T, d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasVerNivel0Paso1GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		es_periodica, delta_T, d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 2);
	procesarAristasHorNivel0Paso1GPU<<<blockGridHor1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasHorNivel0Paso1GPU<<<blockGridHor2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, CFL, epsilon_h, hpos, cvis, 2);
	obtenerEstadosPaso1Nivel0GPU<<<blockGridEst, threadBlockEst>>>(d_datosVolumenesNivel0_1, d_datosNivel0->areaYCosPhi,
		d_datosNivel0->anchoVolumenes, d_datosNivel0->altoVolumenes, d_acumuladorNivel_1, numVolxNivel0, numVolyNivel0,
		delta_T, Hmin, tam_spongeSup, tam_spongeInf, tam_spongeIzq, tam_spongeDer, sea_level);

	cudaMemset(d_acumuladorNivel_1, 0, tam_datosVolDouble2Nivel0);

	// PASO 2
	procesarAristasVerNivel0Paso2GPU<<<blockGridVer1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		es_periodica, delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasVerNivel0Paso2GPU<<<blockGridVer2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->altoVolumenes, numVolxNivel0, numVolyNivel0, borde_izq, borde_der,
		es_periodica, delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 2);
	procesarAristasHorNivel0Paso2GPU<<<blockGridHor1, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 1);
	procesarAristasHorNivel0Paso2GPU<<<blockGridHor2, threadBlockAri>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, numVolxNivel0, numVolyNivel0, borde_sup, borde_inf,
		delta_T, d_acumuladorNivel_1, d_acumuladorNivel_2, CFL, epsilon_h, hpos, cvis, 2);
	obtenerEstadoYDeltaTVolumenesNivel0GPU<<<blockGridEst, threadBlockEst>>>(d_datosVolumenesNivel0_1, d_datosVolumenesNivel0_2,
		d_datosNivel0->areaYCosPhi, d_datosNivel0->anchoVolumenes, d_datosNivel0->altoVolumenes, d_acumuladorNivel_1,
		d_acumuladorNivel_2, d_deltaTVolumenesNivel0, numVolxNivel0, numVolyNivel0, CFL, delta_T, mf0, vmax, hpos,
		epsilon_h, Hmin, tam_spongeSup, tam_spongeInf, tam_spongeIzq, tam_spongeDer, sea_level);

	*tiempo_act += delta_T;

	delta_T = obtenerMinimoReduccion<double>(d_deltaTVolumenesNivel0, numVolxNivel0*numVolyNivel0);

	return delta_T;
}

#endif
