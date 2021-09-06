#ifndef _MATRIZ_H_
#define _MATRIZ_H_

#include "Constantes.hxx"

typedef struct TVec3 {
	double x, y, z;
} TVec3;

typedef double2 TVec2;

__device__ void v_copy2(TVec2 *in, TVec2 *out) {
	out->x = in->x;
	out->y = in->y;
}

__device__ void v_copy3(TVec3 *in, TVec3 *out) {
	out->x = in->x;
	out->y = in->y;
	out->z = in->z;
}

#endif
