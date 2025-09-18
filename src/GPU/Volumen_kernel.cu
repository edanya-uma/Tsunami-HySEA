#ifndef VOLUMEN_KERNEL_H
#define VOLUMEN_KERNEL_H

#include "Constantes.hxx"

__device__ void spongeLayerPaso1(double *h, double dist, double ls, double h0)
{
	double co, val;

	val = dist/ls;
	co = 0.9*val*val;
	*h -= co*((*h)-h0);
	*h *= (*h >= 0);
}

__device__ void spongeLayerPaso2(double *q, double dist, double ls, double h0, double v0)
{
	double co, val;

	val = dist/ls;
	co = 0.9*val*val;
	*q -= co*((*q)-h0*v0);
}

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

__global__ void obtenerEstadosPaso1Nivel0GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_3,
				double *d_anchoVolumenes, double *d_altoVolumenes, double2 *d_acumulador_1, int num_volx,
				int num_voly, double delta_T, double Hmin, int tam_spongeSup, int tam_spongeInf,
				int tam_spongeIzq, int tam_spongeDer, double sea_level)
{
	double2 W1, acum1;
	double val, area;
	double dist, dx, dy;
	double dist_izq, dist_der, dist_sup, dist_inf;
	double h0_sponge, ls;
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
		// Start sponge layer
		h0_sponge = max(W1.y+Hmin+sea_level, 0.0);
		if ((tam_spongeIzq > 0) || (tam_spongeDer > 0) || (tam_spongeSup > 0) || (tam_spongeInf > 0)) {
			if ((pos_x_hebra < tam_spongeIzq) || (pos_x_hebra >= num_volx-tam_spongeDer) || (pos_y_hebra < tam_spongeSup) || (pos_y_hebra >= num_voly-tam_spongeInf)) {
				dx = d_anchoVolumenes[pos_y_hebra];
				dy = d_altoVolumenes[pos_y_hebra];
				dist_izq = (tam_spongeIzq-pos_x_hebra-0.5)*dx;
				if (dist_izq <= 0.0)  dist_izq = 1e30;
				dist_der = (pos_x_hebra-(num_volx-1-tam_spongeDer)-0.5)*dx;
				if (dist_der <= 0.0)  dist_der = 1e30;
				dist_sup = (tam_spongeSup-pos_y_hebra-0.5)*dy;
				if (dist_sup <= 0.0)  dist_sup = 1e30;
				dist_inf = (pos_y_hebra-(num_voly-1-tam_spongeInf)-0.5)*dy;
				if (dist_inf <= 0.0)  dist_inf = 1e30;

				if ((dist_izq <= dist_der) && (dist_izq <= dist_sup) && (dist_izq <= dist_inf)) {
					dist = dist_izq;
					ls = dx*tam_spongeIzq;
				}
				else if ((dist_der <= dist_izq) && (dist_der <= dist_sup) && (dist_der <= dist_inf)) {
					dist = dist_der;
					ls = dx*tam_spongeDer;
				}
				else if ((dist_sup <= dist_izq) && (dist_sup <= dist_der) && (dist_sup <= dist_inf)) {
					dist = dist_sup;
					ls = dy*tam_spongeSup;
				}
				else {
					dist = dist_inf;
					ls = dy*tam_spongeInf;
				}
				spongeLayerPaso1(&(acum1.x), dist, ls, h0_sponge);
			}
		}
		// End sponge layer
		d_datosVolumenes_1[pos] = acum1;
	}
}

__global__ void obtenerEstadoYDeltaTVolumenesNivel0GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double *d_anchoVolumenes, double *d_altoVolumenes, double2 *d_acumulador_1,
				double2 *d_acumulador_2, double *d_deltaTVolumenes, int num_volx, int num_voly, double CFL,
				double delta_T, double mf0, double vmax, double hpos, double epsilon_h, double Hmin,
				int tam_spongeSup, int tam_spongeInf, int tam_spongeIzq, int tam_spongeDer, double sea_level)
{
	double2 W1, W2, acum1, acum2;
	double val, area, paso;
	double dist, dx, dy;
	double h0_sponge, v0_sponge;
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
		// Start sponge layer
		h0_sponge = max(W1.y+Hmin+sea_level, 0.0);
		v0_sponge = 0.0;
		if ((tam_spongeIzq > 0) || (tam_spongeDer > 0) || (tam_spongeSup > 0) || (tam_spongeInf > 0)) {
			dx = d_anchoVolumenes[pos_y_hebra];
			dy = d_altoVolumenes[pos_y_hebra];
			if (pos_x_hebra < tam_spongeIzq) {
				dist = (tam_spongeIzq-pos_x_hebra-0.5)*dx;
				spongeLayerPaso2(&(acum2.x), dist, dx*tam_spongeIzq, h0_sponge, v0_sponge);
			}
			if (pos_x_hebra >= num_volx-tam_spongeDer) {
				dist = (pos_x_hebra-(num_volx-1-tam_spongeDer)-0.5)*dx;
				spongeLayerPaso2(&(acum2.x), dist, dx*tam_spongeDer, h0_sponge, v0_sponge);
			}
			if (pos_y_hebra < tam_spongeSup) {
				dist = (tam_spongeSup-pos_y_hebra-0.5)*dy;
				spongeLayerPaso2(&(acum2.y), dist, dy*tam_spongeSup, h0_sponge, v0_sponge);
			}
			if (pos_y_hebra >= num_voly-tam_spongeInf) {
				dist = (pos_y_hebra-(num_voly-1-tam_spongeInf)-0.5)*dy;
				spongeLayerPaso2(&(acum2.y), dist, dy*tam_spongeInf, h0_sponge, v0_sponge);
			}
		}
		// End sponge layer

		if (W1.x <= hpos) {
			filtroEstado(&W1, &acum2, vmax, delta_T, epsilon_h);
			disImplicita(W1, W2, &W1, &acum2, delta_T, mf0, epsilon_h);
		}

		d_datosVolumenes_2[pos] = acum2;
	}
}

__global__ void escribirVolumenesGuardadoGPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
					double2 *d_datosVolumenesGuardado_1, double2 *d_datosVolumenesGuardado_2,
					int *d_posicionesVolumenesGuardado, int num_puntos_guardar, double epsilon_h)
{
	int pos_hebra;
	int pos;
	double2 W1, W2;
	double h;
	double maxheps, factor;

	pos_hebra = blockIdx.x*NUM_HEBRAS_PUNTOS + threadIdx.x;

	if (pos_hebra < num_puntos_guardar) {
		pos = d_posicionesVolumenesGuardado[pos_hebra];
		if (pos >= 0) {
			W1 = d_datosVolumenes_1[pos];
			W2 = d_datosVolumenes_2[pos];
			h = W1.x;
			if (h < epsilon_h) {
				maxheps = epsilon_h;
				factor = M_SQRT2*h / sqrt(h*h*h*h + maxheps*maxheps*maxheps*maxheps);
			}
			else {
				factor = 1.0/h;
			}
			W2.x = factor*W2.x;
			W2.y = factor*W2.y;
			d_datosVolumenesGuardado_1[pos_hebra] = W1;
			d_datosVolumenesGuardado_2[pos_hebra] = W2;
		}
	}
}

#endif
