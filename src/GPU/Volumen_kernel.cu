#ifndef _VOLUMEN_KERNEL_H_
#define _VOLUMEN_KERNEL_H_

#include <stdio.h>
#include "Constantes.hxx"

__device__ void disImplicita(double2 Want1, double2 Want2, double2 *acum1, double2 *acum2,
							double delta_T, double mf0, double epsilon_h)
{
	double h1, h1m;
	double u1, u1x, u1y;
	double factor, maxheps;

	h1 = acum1->x;
	maxheps = max(h1,epsilon_h);
	h1m = sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
	factor = M_SQRT2*h1/h1m;
	u1x = factor*acum2->x;
	u1y = factor*acum2->y;
	maxheps = max(Want1.x,epsilon_h);
	u1 = M_SQRT2*sqrt(Want2.x*Want2.x + Want2.y*Want2.y)*Want1.x/sqrt(Want1.x*Want1.x*Want1.x*Want1.x + maxheps*maxheps*maxheps*maxheps);

	if (h1 > 0.0)  {
		double c1 = delta_T*mf0*fabs(u1) / (pow(h1,4.0/3.0) + EPSILON);
		factor = h1/(1.0+c1);
		acum2->x = factor*u1x;
		acum2->y = factor*u1y;
	}
}

__device__ void filtroEstado(double2 *acum1, double2 *acum2, double vmax, double delta_T, double epsilon_h)
{
	double aux, aux0;
	double factor, maxheps;

	factor = acum1->x/epsilon_h;
	aux0 = 1.0/(factor*factor*factor*factor + EPSILON);
	aux = exp(-delta_T*aux0);
	acum2->x *= aux;
	acum2->y *= aux;
	if (vmax > 0.0) {
		double h, hm, ux, uy, u;

		h = acum1->x;
		maxheps = max(h,epsilon_h);
		hm = sqrt(h*h*h*h + maxheps*maxheps*maxheps*maxheps);
		factor = M_SQRT2*h/hm;
		ux = factor*acum2->x;
		uy = factor*acum2->y;
		u = sqrt(ux*ux + uy*uy);
		if (u > vmax) {
			ux *= vmax/u;  
			uy *= vmax/u;
		}
		acum2->x = ux*h;
		acum2->y = uy*h;
	}
}

__global__ void inicializarEta1MaximaNivel0GPU(double2 *d_datosVolumenes_1, double *d_eta1_maxima,
					int num_volx, int num_voly, double epsilon_h)
{
	double2 W;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		W = d_datosVolumenes_1[pos];
		d_eta1_maxima[pos] = ((W.x < epsilon_h) ? -1e30 : W.x - W.y);
	}
}

__global__ void actualizarEta1MaximaNivel0GPU(double2 *d_datosVolumenes_1, double *d_eta1_maxima, int num_volx,
					int num_voly, double tiempo_act, double epsilon_h)
{
	double2 W;
	double val, eta1;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;

		eta1 = d_eta1_maxima[pos];
		W = d_datosVolumenes_1[pos];
		val = W.x - W.y;
		if ((val > eta1) && (W.x > epsilon_h))
			d_eta1_maxima[pos] = val;
	}
}

__global__ void obtenerDeltaTVolumenesGPU(double2 *d_datosVolumenesNivel0_3, double2 *d_acumulador_1,
				double *d_deltaTVolumenes, int numVolxNivel0, int numVolyNivel0, float CFL)
{
	double deltaT, area, paso;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < numVolxNivel0) && (pos_y_hebra < numVolyNivel0)) {
		pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
		deltaT = d_acumulador_1[pos].y;
		area = d_datosVolumenesNivel0_3[pos_y_hebra].x;
		paso = ((deltaT < EPSILON) ? 1e30 : (2.0*CFL*area)/deltaT);
		d_deltaTVolumenes[pos] = paso;
	}
}

__global__ void obtenerEstadosPaso1Nivel0GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double2 *d_acumulador_1, int num_volx, int num_voly, double delta_T)
{
	double2 W1, acum1;
	double val, area;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;
		area = d_datosVolumenes_3[pos_y_hebra].x;
		val = delta_T / area;
		acum1 = d_acumulador_1[pos];

		W1 = d_datosVolumenes_1[pos];
		acum1.x = W1.x + val*acum1.x;
		acum1.y = W1.y;
		acum1.x *= (acum1.x > 0.0);
		d_datosVolumenes_1[pos] = acum1;
	}
}

__global__ void obtenerEstadoYDeltaTVolumenesNivel0GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double2 *d_acumulador_1, double2 *d_acumulador_2, double *d_deltaTVolumenes,
				int num_volx, int num_voly, double CFL, double delta_T, double mf0, double vmax, double hpos, double epsilon_h)
{
	double2 W1, W2, acum1, acum2;
	double val, area, paso;
	int pos, pos_x_hebra, pos_y_hebra;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra < num_voly)) {
		pos = pos_y_hebra*num_volx + pos_x_hebra;
		area = d_datosVolumenes_3[pos_y_hebra].x;
		val = delta_T / area;

		acum1 = d_acumulador_1[pos];
		paso = ((acum1.y < EPSILON) ? 1e30 : (2.0*CFL*area)/acum1.y);
		d_deltaTVolumenes[pos] = paso;

		acum2 = d_acumulador_2[pos];
		W1 = d_datosVolumenes_1[pos];
		W2 = d_datosVolumenes_2[pos];
		acum2.x = W2.x + val*acum2.x;
		acum2.y = W2.y + val*acum2.y;

		if (W1.x <= hpos) {
			filtroEstado(&W1, &acum2, vmax, delta_T, epsilon_h);
			disImplicita(W1, W2, &W1, &acum2, delta_T, mf0, epsilon_h);
		}

		d_datosVolumenes_2[pos] = acum2;
	}
}

#endif
