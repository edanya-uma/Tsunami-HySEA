#ifndef _DEFORMACION_H_
#define _DEFORMACION_H_

#include <stdio.h>
#include "Constantes.hxx"
#include "Reduccion_kernel.cu"
#include "prtoxy.cu"
#include "DC3D.cu"

/******************/
/* Okada standard */
/******************/

__global__ void convertirAMenosValorAbsolutoGPU(double *d_in, float *d_out, int n)
{
	int pos = blockIdx.x*NUM_HEBRAS_PUNTOS + threadIdx.x;

	if (pos < n) {
		d_out[pos] = (float) (-fabs(d_in[pos]));
	}
}

__global__ void truncarDeformacionGPU(double *d_def, int num_volx, int num_voly, double crop_value)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double U_Z;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;
		U_Z = d_def[pos];
		if (fabs(U_Z) < crop_value)
			d_def[pos] = 0.0;
	}
}

__global__ void sumarDeformacionADatosGPU(double2 *d_datosVolumenes_1, double *d_def, double *d_eta1Inicial,
				int num_volx, int num_voly)
{
	int pos, pos_x_hebra, pos_y_hebra;
	double U_Z;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;
		U_Z = d_def[pos];
		d_datosVolumenes_1[pos].y -= U_Z;
		d_eta1Inicial[pos] += U_Z;
	}
}

__global__ void aplicarOkadaStandardGPU(double2 *d_datosVolumenes_1, double *d_def, int num_volx, int num_voly, double lon_ini,
				double incx, double lat_ini, double incy, double LON_C_ent, double LAT_C_ent, double DEPTH_C_ent, double FAULT_L,
				double FAULT_W, double STRIKE, double DIP_ent, double RAKE, double SLIP, double H)
{
	double LON_P, LAT_P;
	double LON_C, LAT_C;
	double DEPTH_C, DIP;
	double S_RAKE;
	double C_RAKE;
	double S_STRIKE;
	double C_STRIKE;
	double Z;
	double AL1, AL2, AW1, AW2;
	double X_OKA, Y_OKA;
	double XP, YP;
	double RAD = M_PI/180.0;
	double alfa = 2.0/3.0;
	int i0 = 0;
	int IRET;
	double U_X, U_Y, U_Z, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ;
	double DISL1, DISL2, DISL3;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;
		LON_C = LON_C_ent;
		LAT_C = LAT_C_ent;
		DEPTH_C = DEPTH_C_ent;
		DIP = DIP_ent;

		LON_P = lon_ini + pos_x_hebra*incx;
		LAT_P = lat_ini + pos_y_hebra*incy;

		S_RAKE = sin(RAD*RAKE);
		C_RAKE = cos(RAD*RAKE);

		S_STRIKE = sin(RAD*STRIKE);
		C_STRIKE = cos(RAD*STRIKE);

		DISL2 = SLIP*S_RAKE;
		DISL1 = SLIP*C_RAKE;
		DISL3 = 0.0;

		Z = 0.0;
		AL1 = -0.5*FAULT_L;
		AL2 = 0.5*FAULT_L;
		AW1 = -0.5*FAULT_W;
		AW2 = 0.5*FAULT_W;

		prtoxy_(&LAT_P, &LON_P, &LAT_C, &LON_C, &XP, &YP, &i0);
		X_OKA = XP*S_STRIKE + YP*C_STRIKE;
		Y_OKA = -XP*C_STRIKE + YP*S_STRIKE;
		dc3d_(&alfa, &X_OKA, &Y_OKA, &Z, &DEPTH_C, &DIP, &AL1, &AL2, &AW1, &AW2, &DISL1, &DISL2, &DISL3,
			&U_X, &U_Y, &U_Z, &UXX, &UYX, &UZX, &UXY, &UYY, &UZY, &UXZ, &UYZ, &UZZ, &IRET);

		U_Z /= H;
		d_def[pos] = U_Z;
	}
}

void aplicarOkada(double2 *d_datosVolumenes_1, double *d_eta1Inicial, int crop_flag, double crop_value,
		double *d_deltaTVolumenes, float *d_vec, int num_volx, int num_voly, double lon_ini, double incx,
		double lat_ini, double incy, double LON_C, double LAT_C, double DEPTH_C, double FAULT_L, double FAULT_W,
		double STRIKE, double DIP, double RAKE, double SLIP, dim3 blockGridEstNivel, dim3 threadBlockEst, double H)
{
	int num_volumenes = num_volx*num_voly;
	dim3 blockGridVec(iDivUp(num_volumenes, NUM_HEBRAS_PUNTOS), 1);
	dim3 threadBlockVec(NUM_HEBRAS_PUNTOS, 1);
	float def_max;
	double crop_value_final;

	aplicarOkadaStandardGPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenes_1, d_deltaTVolumenes,
		num_volx, num_voly, lon_ini, incx, lat_ini, incy, LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W,
		STRIKE, DIP, RAKE, SLIP, H);
	if (crop_flag == CROP_RELATIVE) {
		convertirAMenosValorAbsolutoGPU<<<blockGridVec, threadBlockVec>>>(d_deltaTVolumenes, d_vec, num_volumenes);
		def_max = -obtenerMinimoReduccion<float>(d_vec, num_volumenes);
		crop_value_final = crop_value*((double) def_max);
		truncarDeformacionGPU<<<blockGridEstNivel, threadBlockEst>>>(d_deltaTVolumenes, num_volx, num_voly, crop_value_final);
	}
	else if (crop_flag == CROP_ABSOLUTE) {
		truncarDeformacionGPU<<<blockGridEstNivel, threadBlockEst>>>(d_deltaTVolumenes, num_volx, num_voly, crop_value);
	}
	sumarDeformacionADatosGPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenes_1, d_deltaTVolumenes,
		d_eta1Inicial, num_volx, num_voly);
}

#endif
