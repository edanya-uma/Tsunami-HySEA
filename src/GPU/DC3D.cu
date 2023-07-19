#ifndef DC3D_H
#define DC3D_H

/* Common Block Declarations */

typedef struct {
	double dummy[5], sd, cd, sdsd, cdcd, sdcd, s2d, c2d;
} Tipo_1;

typedef struct {
	double alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, 
		s2d, c2d;
} Tipo_2;

typedef struct {
    double xi2, et2, q2, r__, r2, r3, r5, y, d__, tt, alx, ale, x11, y11, 
	    x32, y32, ey, ez, fy, fz, gy, gz, hy, hz;
} Tipo_c2_1;

__device__ int dccon0_(double *alpha, double *dip, Tipo_2 *c0_2)
{
    /* Initialized data */

    double f0 = 0.0;
    double f1 = 1.0;
    double f2 = 2.0;
    double pi2 = 6.283185307179586;
    double eps = 1e-6;

    /* Local variables */
    double p18;


/* ******************************************************************* */
/* *****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    ***** */
/* ******************************************************************* */

/* ***** INPUT */
/* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU) */
/* *****   DIP   : DIP-ANGLE (DEGREE) */
/* ### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO */

/* ----- */
    c0_2->alp1 = (f1 - *alpha) / f2;
    c0_2->alp2 = *alpha / f2;
    c0_2->alp3 = (f1 - *alpha) / *alpha;
    c0_2->alp4 = f1 - *alpha;
    c0_2->alp5 = *alpha;
/* ----- */
    p18 = pi2 / 360.0;
    c0_2->sd = sin(*dip * p18);
    c0_2->cd = cos(*dip * p18);
    if (abs(c0_2->cd) < eps) {
	c0_2->cd = f0;
	if (c0_2->sd > f0) {
	    c0_2->sd = f1;
	}
	if (c0_2->sd < f0) {
	    c0_2->sd = -f1;
	}
    }
    c0_2->sdsd = c0_2->sd * c0_2->sd;
    c0_2->cdcd = c0_2->cd * c0_2->cd;
    c0_2->sdcd = c0_2->sd * c0_2->cd;
    c0_2->s2d = f2 * c0_2->sdcd;
    c0_2->c2d = c0_2->cdcd - c0_2->sdsd;
    return 0;
} /* dccon0_ */

__device__ int dccon2_(double *xi, double *et, double *q, 
	double *sd, double *cd, int *kxi, int *ket, Tipo_c2_1 *c2_1)
{
    /* Initialized data */

    double f0 = 0.0;
    double f1 = 1.0;
    double f2 = 2.0;
    double eps = 1e-6;

    /* Local variables */
    double ret, rxi;


/* ********************************************************************** */
/* *****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** */
/* ********************************************************************** */

/* ***** INPUT */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
/* *****   SD,CD   : SIN, COS OF DIP-ANGLE */
/* *****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY */

/* ### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0 */

/* ----- */
    if (abs(*xi) < eps) {
	*xi = f0;
    }
    if (abs(*et) < eps) {
	*et = f0;
    }
    if (abs(*q) < eps) {
	*q = f0;
    }
    c2_1->xi2 = *xi * *xi;
    c2_1->et2 = *et * *et;
    c2_1->q2 = *q * *q;
    c2_1->r2 = c2_1->xi2 + c2_1->et2 + c2_1->q2;
    c2_1->r__ = sqrt(c2_1->r2);
    if (c2_1->r__ == f0) {
	return 0;
    }
    c2_1->r3 = c2_1->r__ * c2_1->r2;
    c2_1->r5 = c2_1->r3 * c2_1->r2;
    c2_1->y = *et * *cd + *q * *sd;
    c2_1->d__ = *et * *sd - *q * *cd;
/* ----- */
    if (*q == f0) {
	c2_1->tt = f0;
    } else {
	c2_1->tt = atan(*xi * *et / (*q * c2_1->r__));
    }
/* ----- */
    if (*kxi == 1) {
	c2_1->alx = -log(c2_1->r__ - *xi);
	c2_1->x11 = f0;
	c2_1->x32 = f0;
    } else {
	rxi = c2_1->r__ + *xi;
	c2_1->alx = log(rxi);
	c2_1->x11 = f1 / (c2_1->r__ * rxi);
	c2_1->x32 = (c2_1->r__ + rxi) * c2_1->x11 * c2_1->x11 / c2_1->r__;
    }
/* ----- */
    if (*ket == 1) {
	c2_1->ale = -log(c2_1->r__ - *et);
	c2_1->y11 = f0;
	c2_1->y32 = f0;
    } else {
	ret = c2_1->r__ + *et;
	c2_1->ale = log(ret);
	c2_1->y11 = f1 / (c2_1->r__ * ret);
	c2_1->y32 = (c2_1->r__ + ret) * c2_1->y11 * c2_1->y11 / c2_1->r__;
    }
/* ----- */
    c2_1->ey = *sd / c2_1->r__ - c2_1->y * *q / c2_1->r3;
    c2_1->ez = *cd / c2_1->r__ + c2_1->d__ * *q / c2_1->r3;
    c2_1->fy = c2_1->d__ / c2_1->r3 + c2_1->xi2 * c2_1->y32 * *sd;
    c2_1->fz = c2_1->y / c2_1->r3 + c2_1->xi2 * c2_1->y32 * *cd;
    c2_1->gy = f2 * c2_1->x11 * *sd - c2_1->y * *q * c2_1->x32;
    c2_1->gz = f2 * c2_1->x11 * *cd + c2_1->d__ * *q * c2_1->x32;
    c2_1->hy = c2_1->d__ * *q * c2_1->x32 + *xi * *q * c2_1->y32 * *sd;
    c2_1->hz = c2_1->y * *q * c2_1->x32 + *xi * *q * c2_1->y32 * *cd;
    return 0;
} /* dccon2_ */

__device__ int ua_(double *xi, double *et, double *q, double *disl1,
    double *disl2, double *disl3, double *u, Tipo_2 *c0_2, Tipo_c2_1 *c2_1)
{
    /* Initialized data */

    double f0 = 0.0;
    double f2 = 2.0;
    double pi2 = 6.283185307179586;

    int i__;
    double du[12], qx, qy, xy;


/* ******************************************************************** */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             ***** */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
/* ******************************************************************** */

/* ***** INPUT */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
/* ***** OUTPUT */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /* Parameter adjustments */
    --u;

    /* Function Body */
/* ----- */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
	u[i__] = f0;
    }
    xy = *xi * c2_1->y11;
    qx = *q * c2_1->x11;
    qy = *q * c2_1->y11;
/* ====================================== */
/* =====  STRIKE-SLIP CONTRIBUTION  ===== */
/* ====================================== */
    if (*disl1 != f0) {
	du[0] = c2_1->tt / f2 + c0_2->alp2 * *xi * qy;
	du[1] = c0_2->alp2 * *q / c2_1->r__;
	du[2] = c0_2->alp1 * c2_1->ale - c0_2->alp2 * *q * qy;
	du[3] = -c0_2->alp1 * qy - c0_2->alp2 * c2_1->xi2 * *q * c2_1->y32;
	du[4] = -c0_2->alp2 * *xi * *q / c2_1->r3;
	du[5] = c0_2->alp1 * xy + c0_2->alp2 * *xi * c2_1->q2 * c2_1->y32;
	du[6] = c0_2->alp1 * xy * c0_2->sd + c0_2->alp2 * *xi * c2_1->fy + 
		c2_1->d__ / f2 * c2_1->x11;
	du[7] = c0_2->alp2 * c2_1->ey;
	du[8] = c0_2->alp1 * (c0_2->cd / c2_1->r__ + qy * c0_2->sd) - c0_2->alp2 * 
		*q * c2_1->fy;
	du[9] = c0_2->alp1 * xy * c0_2->cd + c0_2->alp2 * *xi * c2_1->fz + c2_1->y 
		/ f2 * c2_1->x11;
	du[10] = c0_2->alp2 * c2_1->ez;
	du[11] = -c0_2->alp1 * (c0_2->sd / c2_1->r__ - qy * c0_2->cd) - c0_2->alp2 
		* *q * c2_1->fz;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
	    u[i__] += *disl1 / pi2 * du[i__ - 1];
	}
    }
/* ====================================== */
/* =====    DIP-SLIP CONTRIBUTION   ===== */
/* ====================================== */
    if (*disl2 != f0) {
	du[0] = c0_2->alp2 * *q / c2_1->r__;
	du[1] = c2_1->tt / f2 + c0_2->alp2 * *et * qx;
	du[2] = c0_2->alp1 * c2_1->alx - c0_2->alp2 * *q * qx;
	du[3] = -c0_2->alp2 * *xi * *q / c2_1->r3;
	du[4] = -qy / f2 - c0_2->alp2 * *et * *q / c2_1->r3;
	du[5] = c0_2->alp1 / c2_1->r__ + c0_2->alp2 * c2_1->q2 / c2_1->r3;
	du[6] = c0_2->alp2 * c2_1->ey;
	du[7] = c0_2->alp1 * c2_1->d__ * c2_1->x11 + xy / f2 * c0_2->sd + 
		c0_2->alp2 * *et * c2_1->gy;
	du[8] = c0_2->alp1 * c2_1->y * c2_1->x11 - c0_2->alp2 * *q * c2_1->gy;
	du[9] = c0_2->alp2 * c2_1->ez;
	du[10] = c0_2->alp1 * c2_1->y * c2_1->x11 + xy / f2 * c0_2->cd + 
		c0_2->alp2 * *et * c2_1->gz;
	du[11] = -c0_2->alp1 * c2_1->d__ * c2_1->x11 - c0_2->alp2 * *q * c2_1->gz;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
	    u[i__] += *disl2 / pi2 * du[i__ - 1];
	}
    }
/* ======================================== */
/* =====  TENSILE-FAULT CONTRIBUTION  ===== */
/* ======================================== */
    if (*disl3 != f0) {
	du[0] = -c0_2->alp1 * c2_1->ale - c0_2->alp2 * *q * qy;
	du[1] = -c0_2->alp1 * c2_1->alx - c0_2->alp2 * *q * qx;
	du[2] = c2_1->tt / f2 - c0_2->alp2 * (*et * qx + *xi * qy);
	du[3] = -c0_2->alp1 * xy + c0_2->alp2 * *xi * c2_1->q2 * c2_1->y32;
	du[4] = -c0_2->alp1 / c2_1->r__ + c0_2->alp2 * c2_1->q2 / c2_1->r3;
	du[5] = -c0_2->alp1 * qy - c0_2->alp2 * *q * c2_1->q2 * c2_1->y32;
	du[6] = -c0_2->alp1 * (c0_2->cd / c2_1->r__ + qy * c0_2->sd) - c0_2->alp2 *
		 *q * c2_1->fy;
	du[7] = -c0_2->alp1 * c2_1->y * c2_1->x11 - c0_2->alp2 * *q * c2_1->gy;
	du[8] = c0_2->alp1 * (c2_1->d__ * c2_1->x11 + xy * c0_2->sd) + c0_2->alp2 *
		 *q * c2_1->hy;
	du[9] = c0_2->alp1 * (c0_2->sd / c2_1->r__ - qy * c0_2->cd) - c0_2->alp2 * 
		*q * c2_1->fz;
	du[10] = c0_2->alp1 * c2_1->d__ * c2_1->x11 - c0_2->alp2 * *q * c2_1->gz;
	du[11] = c0_2->alp1 * (c2_1->y * c2_1->x11 + xy * c0_2->cd) + c0_2->alp2 * 
		*q * c2_1->hz;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
	    u[i__] += *disl3 / pi2 * du[i__ - 1];
	}
    }
    return 0;
} /* ua_ */

__device__ /* Subroutine */ int ub_(double *xi, double *et, double *q, double *disl1,
    double *disl2, double *disl3, double *u, Tipo_2 *c0_2, Tipo_c2_1 *c2_1)
{
    /* Initialized data */

    double f0 = 0.0;
    double f1 = 1.0;
    double f2 = 2.0;
    double pi2 = 6.283185307179586;

    /* Local variables */
    int i__;
    double x, d11, rd, du[12], qx, qy, xy, ai1, ai2, aj2, ai4, ai3,
	     aj5, ak1, ak3, aj3, aj6, ak2, ak4, aj1, rd2, aj4;


/* ******************************************************************** */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             ***** */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
/* ******************************************************************** */

/* ***** INPUT */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
/* ***** OUTPUT */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /* Parameter adjustments */
    --u;

    /* Function Body */
/* ----- */
    rd = c2_1->r__ + c2_1->d__;
    d11 = f1 / (c2_1->r__ * rd);
    aj2 = *xi * c2_1->y / rd * d11;
    aj5 = -(c2_1->d__ + c2_1->y * c2_1->y / rd) * d11;
    if (c0_2->cd != f0) {
	if (*xi == f0) {
	    ai4 = f0;
	} else {
	    x = sqrt(c2_1->xi2 + c2_1->q2);
	    ai4 = f1 / c0_2->cdcd * (*xi / rd * c0_2->sdcd + f2 * atan((*et * (
		    x + *q * c0_2->cd) + x * (c2_1->r__ + x) * c0_2->sd) / (*xi *
		     (c2_1->r__ + x) * c0_2->cd)));
	}
	ai3 = (c2_1->y * c0_2->cd / rd - c2_1->ale + c0_2->sd * log(rd)) / 
		c0_2->cdcd;
	ak1 = *xi * (d11 - c2_1->y11 * c0_2->sd) / c0_2->cd;
	ak3 = (*q * c2_1->y11 - c2_1->y * d11) / c0_2->cd;
	aj3 = (ak1 - aj2 * c0_2->sd) / c0_2->cd;
	aj6 = (ak3 - aj5 * c0_2->sd) / c0_2->cd;
    } else {
	rd2 = rd * rd;
	ai3 = (*et / rd + c2_1->y * *q / rd2 - c2_1->ale) / f2;
	ai4 = *xi * c2_1->y / rd2 / f2;
	ak1 = *xi * *q / rd * d11;
	ak3 = c0_2->sd / rd * (c2_1->xi2 * d11 - f1);
	aj3 = -(*xi) / rd2 * (c2_1->q2 * d11 - f1 / f2);
	aj6 = -c2_1->y / rd2 * (c2_1->xi2 * d11 - f1 / f2);
    }
/* ----- */
    xy = *xi * c2_1->y11;
    ai1 = -(*xi) / rd * c0_2->cd - ai4 * c0_2->sd;
    ai2 = log(rd) + ai3 * c0_2->sd;
    ak2 = f1 / c2_1->r__ + ak3 * c0_2->sd;
    ak4 = xy * c0_2->cd - ak1 * c0_2->sd;
    aj1 = aj5 * c0_2->cd - aj6 * c0_2->sd;
    aj4 = -xy - aj2 * c0_2->cd + aj3 * c0_2->sd;
/* ===== */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
	u[i__] = f0;
    }
    qx = *q * c2_1->x11;
    qy = *q * c2_1->y11;
/* ====================================== */
/* =====  STRIKE-SLIP CONTRIBUTION  ===== */
/* ====================================== */
    if (*disl1 != f0) {
	du[0] = -(*xi) * qy - c2_1->tt - c0_2->alp3 * ai1 * c0_2->sd;
	du[1] = -(*q) / c2_1->r__ + c0_2->alp3 * c2_1->y / rd * c0_2->sd;
	du[2] = *q * qy - c0_2->alp3 * ai2 * c0_2->sd;
	du[3] = c2_1->xi2 * *q * c2_1->y32 - c0_2->alp3 * aj1 * c0_2->sd;
	du[4] = *xi * *q / c2_1->r3 - c0_2->alp3 * aj2 * c0_2->sd;
	du[5] = -(*xi) * c2_1->q2 * c2_1->y32 - c0_2->alp3 * aj3 * c0_2->sd;
	du[6] = -(*xi) * c2_1->fy - c2_1->d__ * c2_1->x11 + c0_2->alp3 * (xy + 
		aj4) * c0_2->sd;
	du[7] = -c2_1->ey + c0_2->alp3 * (f1 / c2_1->r__ + aj5) * c0_2->sd;
	du[8] = *q * c2_1->fy - c0_2->alp3 * (qy - aj6) * c0_2->sd;
	du[9] = -(*xi) * c2_1->fz - c2_1->y * c2_1->x11 + c0_2->alp3 * ak1 * 
		c0_2->sd;
	du[10] = -c2_1->ez + c0_2->alp3 * c2_1->y * d11 * c0_2->sd;
	du[11] = *q * c2_1->fz + c0_2->alp3 * ak2 * c0_2->sd;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
	    u[i__] += *disl1 / pi2 * du[i__ - 1];
	}
    }
/* ====================================== */
/* =====    DIP-SLIP CONTRIBUTION   ===== */
/* ====================================== */
    if (*disl2 != f0) {
	du[0] = -(*q) / c2_1->r__ + c0_2->alp3 * ai3 * c0_2->sdcd;
	du[1] = -(*et) * qx - c2_1->tt - c0_2->alp3 * *xi / rd * c0_2->sdcd;
	du[2] = *q * qx + c0_2->alp3 * ai4 * c0_2->sdcd;
	du[3] = *xi * *q / c2_1->r3 + c0_2->alp3 * aj4 * c0_2->sdcd;
	du[4] = *et * *q / c2_1->r3 + qy + c0_2->alp3 * aj5 * c0_2->sdcd;
	du[5] = -c2_1->q2 / c2_1->r3 + c0_2->alp3 * aj6 * c0_2->sdcd;
	du[6] = -c2_1->ey + c0_2->alp3 * aj1 * c0_2->sdcd;
	du[7] = -(*et) * c2_1->gy - xy * c0_2->sd + c0_2->alp3 * aj2 * c0_2->sdcd;
	du[8] = *q * c2_1->gy + c0_2->alp3 * aj3 * c0_2->sdcd;
	du[9] = -c2_1->ez - c0_2->alp3 * ak3 * c0_2->sdcd;
	du[10] = -(*et) * c2_1->gz - xy * c0_2->cd - c0_2->alp3 * *xi * d11 * 
		c0_2->sdcd;
	du[11] = *q * c2_1->gz - c0_2->alp3 * ak4 * c0_2->sdcd;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
	    u[i__] += *disl2 / pi2 * du[i__ - 1];
	}
    }
/* ======================================== */
/* =====  TENSILE-FAULT CONTRIBUTION  ===== */
/* ======================================== */
    if (*disl3 != f0) {
	du[0] = *q * qy - c0_2->alp3 * ai3 * c0_2->sdsd;
	du[1] = *q * qx + c0_2->alp3 * *xi / rd * c0_2->sdsd;
	du[2] = *et * qx + *xi * qy - c2_1->tt - c0_2->alp3 * ai4 * c0_2->sdsd;
	du[3] = -(*xi) * c2_1->q2 * c2_1->y32 - c0_2->alp3 * aj4 * c0_2->sdsd;
	du[4] = -c2_1->q2 / c2_1->r3 - c0_2->alp3 * aj5 * c0_2->sdsd;
	du[5] = *q * c2_1->q2 * c2_1->y32 - c0_2->alp3 * aj6 * c0_2->sdsd;
	du[6] = *q * c2_1->fy - c0_2->alp3 * aj1 * c0_2->sdsd;
	du[7] = *q * c2_1->gy - c0_2->alp3 * aj2 * c0_2->sdsd;
	du[8] = -(*q) * c2_1->hy - c0_2->alp3 * aj3 * c0_2->sdsd;
	du[9] = *q * c2_1->fz + c0_2->alp3 * ak3 * c0_2->sdsd;
	du[10] = *q * c2_1->gz + c0_2->alp3 * *xi * d11 * c0_2->sdsd;
	du[11] = -(*q) * c2_1->hz + c0_2->alp3 * ak4 * c0_2->sdsd;
	for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
	    u[i__] += *disl3 / pi2 * du[i__ - 1];
	}
    }
    return 0;
} /* ub_ */

__device__ int uc_(double *xi, double *et, double *q, double *z__,
    double *disl1, double *disl2, double *disl3, double *u, Tipo_2 *c0_2, Tipo_c2_1 *c2_1)
{
    /* Initialized data */

    double f0 = 0.0;
    double f1 = 1.0;
    double f2 = 2.0;
    double f3 = 3.0;
    double pi2 = 6.283185307179586;

    double c__, h__;
    int i__;
    double y0, z0, du[12], x53, y53, z32, z53, qq, qy, qr, xy, 
	    yy0, cdr, ppy, ppz, qqy, qqz;


/* ******************************************************************** */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             ***** */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
/* ******************************************************************** */

/* ***** INPUT */
/* *****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
/* ***** OUTPUT */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /* Parameter adjustments */
    --u;

    /* Function Body */
/* ----- */
    c__ = c2_1->d__ + *z__;
    x53 = (c2_1->r2 * 8. + c2_1->r__ * 9. * *xi + f3 * c2_1->xi2) * c2_1->x11 * 
	    c2_1->x11 * c2_1->x11 / c2_1->r2;
    y53 = (c2_1->r2 * 8. + c2_1->r__ * 9. * *et + f3 * c2_1->et2) * c2_1->y11 * 
	    c2_1->y11 * c2_1->y11 / c2_1->r2;
    h__ = *q * c0_2->cd - *z__;
    z32 = c0_2->sd / c2_1->r3 - h__ * c2_1->y32;
    z53 = f3 * c0_2->sd / c2_1->r5 - h__ * y53;
    y0 = c2_1->y11 - c2_1->xi2 * c2_1->y32;
    z0 = z32 - c2_1->xi2 * z53;
    ppy = c0_2->cd / c2_1->r3 + *q * c2_1->y32 * c0_2->sd;
    ppz = c0_2->sd / c2_1->r3 - *q * c2_1->y32 * c0_2->cd;
    qq = *z__ * c2_1->y32 + z32 + z0;
    qqy = f3 * c__ * c2_1->d__ / c2_1->r5 - qq * c0_2->sd;
    qqz = f3 * c__ * c2_1->y / c2_1->r5 - qq * c0_2->cd + *q * c2_1->y32;
    xy = *xi * c2_1->y11;
    qy = *q * c2_1->y11;
    qr = f3 * *q / c2_1->r5;
    cdr = (c__ + c2_1->d__) / c2_1->r3;
    yy0 = c2_1->y / c2_1->r3 - y0 * c0_2->cd;
/* ===== */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
	u[i__] = f0;
    }
/* ====================================== */
/* =====  STRIKE-SLIP CONTRIBUTION  ===== */
/* ====================================== */
    if (*disl1 != f0) {
	du[0] = c0_2->alp4 * xy * c0_2->cd - c0_2->alp5 * *xi * *q * z32;
	du[1] = c0_2->alp4 * (c0_2->cd / c2_1->r__ + f2 * qy * c0_2->sd) - 
		c0_2->alp5 * c__ * *q / c2_1->r3;
	du[2] = c0_2->alp4 * qy * c0_2->cd - c0_2->alp5 * (c__ * *et / c2_1->r3 - 
		*z__ * c2_1->y11 + c2_1->xi2 * z32);
	du[3] = c0_2->alp4 * y0 * c0_2->cd - c0_2->alp5 * *q * z0;
	du[4] = -c0_2->alp4 * *xi * (c0_2->cd / c2_1->r3 + f2 * *q * c2_1->y32 * 
		c0_2->sd) + c0_2->alp5 * c__ * *xi * qr;
	du[5] = -c0_2->alp4 * *xi * *q * c2_1->y32 * c0_2->cd + c0_2->alp5 * *xi *
		 (f3 * c__ * *et / c2_1->r5 - qq);
	du[6] = -c0_2->alp4 * *xi * ppy * c0_2->cd - c0_2->alp5 * *xi * qqy;
	du[7] = c0_2->alp4 * f2 * (c2_1->d__ / c2_1->r3 - y0 * c0_2->sd) * 
		c0_2->sd - c2_1->y / c2_1->r3 * c0_2->cd - c0_2->alp5 * (cdr * 
		c0_2->sd - *et / c2_1->r3 - c__ * c2_1->y * qr);
	du[8] = -c0_2->alp4 * *q / c2_1->r3 + yy0 * c0_2->sd + c0_2->alp5 * (cdr *
		 c0_2->cd + c__ * c2_1->d__ * qr - (y0 * c0_2->cd + *q * z0) * 
		c0_2->sd);
	du[9] = c0_2->alp4 * *xi * ppz * c0_2->cd - c0_2->alp5 * *xi * qqz;
	du[10] = c0_2->alp4 * f2 * (c2_1->y / c2_1->r3 - y0 * c0_2->cd) * c0_2->sd 
		+ c2_1->d__ / c2_1->r3 * c0_2->cd - c0_2->alp5 * (cdr * c0_2->cd + 
		c__ * c2_1->d__ * qr);
	du[11] = yy0 * c0_2->cd - c0_2->alp5 * (cdr * c0_2->sd - c__ * c2_1->y * 
		qr - y0 * c0_2->sdsd + *q * z0 * c0_2->cd);
	for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
	    u[i__] += *disl1 / pi2 * du[i__ - 1];
	}
    }
/* ====================================== */
/* =====    DIP-SLIP CONTRIBUTION   ===== */
/* ====================================== */
    if (*disl2 != f0) {
	du[0] = c0_2->alp4 * c0_2->cd / c2_1->r__ - qy * c0_2->sd - c0_2->alp5 * 
		c__ * *q / c2_1->r3;
	du[1] = c0_2->alp4 * c2_1->y * c2_1->x11 - c0_2->alp5 * c__ * *et * *q * 
		c2_1->x32;
	du[2] = -c2_1->d__ * c2_1->x11 - xy * c0_2->sd - c0_2->alp5 * c__ * (
		c2_1->x11 - c2_1->q2 * c2_1->x32);
	du[3] = -c0_2->alp4 * *xi / c2_1->r3 * c0_2->cd + c0_2->alp5 * c__ * *xi *
		 qr + *xi * *q * c2_1->y32 * c0_2->sd;
	du[4] = -c0_2->alp4 * c2_1->y / c2_1->r3 + c0_2->alp5 * c__ * *et * qr;
	du[5] = c2_1->d__ / c2_1->r3 - y0 * c0_2->sd + c0_2->alp5 * c__ / c2_1->r3 
		* (f1 - f3 * c2_1->q2 / c2_1->r2);
	du[6] = -c0_2->alp4 * *et / c2_1->r3 + y0 * c0_2->sdsd - c0_2->alp5 * (
		cdr * c0_2->sd - c__ * c2_1->y * qr);
	du[7] = c0_2->alp4 * (c2_1->x11 - c2_1->y * c2_1->y * c2_1->x32) - 
		c0_2->alp5 * c__ * ((c2_1->d__ + f2 * *q * c0_2->cd) * c2_1->x32 
		- c2_1->y * *et * *q * x53);
	du[8] = *xi * ppy * c0_2->sd + c2_1->y * c2_1->d__ * c2_1->x32 + 
		c0_2->alp5 * c__ * ((c2_1->y + f2 * *q * c0_2->sd) * c2_1->x32 - 
		c2_1->y * c2_1->q2 * x53);
	du[9] = -(*q) / c2_1->r3 + y0 * c0_2->sdcd - c0_2->alp5 * (cdr * c0_2->cd 
		+ c__ * c2_1->d__ * qr);
	du[10] = c0_2->alp4 * c2_1->y * c2_1->d__ * c2_1->x32 - c0_2->alp5 * c__ * 
		((c2_1->y - f2 * *q * c0_2->sd) * c2_1->x32 + c2_1->d__ * *et * *
		q * x53);
	du[11] = -(*xi) * ppz * c0_2->sd + c2_1->x11 - c2_1->d__ * c2_1->d__ * 
		c2_1->x32 - c0_2->alp5 * c__ * ((c2_1->d__ - f2 * *q * c0_2->cd) *
		 c2_1->x32 - c2_1->d__ * c2_1->q2 * x53);
	for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
	    u[i__] += *disl2 / pi2 * du[i__ - 1];
	}
    }
/* ======================================== */
/* =====  TENSILE-FAULT CONTRIBUTION  ===== */
/* ======================================== */
    if (*disl3 != f0) {
	du[0] = -c0_2->alp4 * (c0_2->sd / c2_1->r__ + qy * c0_2->cd) - c0_2->alp5 *
		 (*z__ * c2_1->y11 - c2_1->q2 * z32);
	du[1] = c0_2->alp4 * f2 * xy * c0_2->sd + c2_1->d__ * c2_1->x11 - 
		c0_2->alp5 * c__ * (c2_1->x11 - c2_1->q2 * c2_1->x32);
	du[2] = c0_2->alp4 * (c2_1->y * c2_1->x11 + xy * c0_2->cd) + c0_2->alp5 * *
		q * (c__ * *et * c2_1->x32 + *xi * z32);
	du[3] = c0_2->alp4 * *xi / c2_1->r3 * c0_2->sd + *xi * *q * c2_1->y32 * 
		c0_2->cd + c0_2->alp5 * *xi * (f3 * c__ * *et / c2_1->r5 - f2 * 
		z32 - z0);
	du[4] = c0_2->alp4 * f2 * y0 * c0_2->sd - c2_1->d__ / c2_1->r3 + 
		c0_2->alp5 * c__ / c2_1->r3 * (f1 - f3 * c2_1->q2 / c2_1->r2);
	du[5] = -c0_2->alp4 * yy0 - c0_2->alp5 * (c__ * *et * qr - *q * z0);
	du[6] = c0_2->alp4 * (*q / c2_1->r3 + y0 * c0_2->sdcd) + c0_2->alp5 * (*
		z__ / c2_1->r3 * c0_2->cd + c__ * c2_1->d__ * qr - *q * z0 * 
		c0_2->sd);
	du[7] = -c0_2->alp4 * f2 * *xi * ppy * c0_2->sd - c2_1->y * c2_1->d__ * 
		c2_1->x32 + c0_2->alp5 * c__ * ((c2_1->y + f2 * *q * c0_2->sd) * 
		c2_1->x32 - c2_1->y * c2_1->q2 * x53);
	du[8] = -c0_2->alp4 * (*xi * ppy * c0_2->cd - c2_1->x11 + c2_1->y * 
		c2_1->y * c2_1->x32) + c0_2->alp5 * (c__ * ((c2_1->d__ + f2 * *q *
		 c0_2->cd) * c2_1->x32 - c2_1->y * *et * *q * x53) + *xi * qqy);
	du[9] = -(*et) / c2_1->r3 + y0 * c0_2->cdcd - c0_2->alp5 * (*z__ / 
		c2_1->r3 * c0_2->sd - c__ * c2_1->y * qr - y0 * c0_2->sdsd + *q * 
		z0 * c0_2->cd);
	du[10] = c0_2->alp4 * f2 * *xi * ppz * c0_2->sd - c2_1->x11 + c2_1->d__ * 
		c2_1->d__ * c2_1->x32 - c0_2->alp5 * c__ * ((c2_1->d__ - f2 * *q *
		 c0_2->cd) * c2_1->x32 - c2_1->d__ * c2_1->q2 * x53);
	du[11] = c0_2->alp4 * (*xi * ppz * c0_2->cd + c2_1->y * c2_1->d__ * 
		c2_1->x32) + c0_2->alp5 * (c__ * ((c2_1->y - f2 * *q * c0_2->sd) *
		 c2_1->x32 + c2_1->d__ * *et * *q * x53) + *xi * qqz);
	for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
	    u[i__] += *disl3 / pi2 * du[i__ - 1];
	}
    }
    return 0;
} /* uc_ */

__device__ int dc3d_(double *alpha, double *x, double *y, 
	double *z__, double *depth, double *dip, double *al1, 
	double *al2, double *aw1, double *aw2, double *disl1, 
	double *disl2, double *disl3, double *ux, double *uy, 
	double *uz, double *uxx, double *uyx, double *uzx, 
	double *uxy, double *uyy, double *uzy, double *uxz, 
	double *uyz, double *uzz, int *iret)
{
    /* Initialized data */

    double f0 = 0.0;
    double eps = 1e-6;

    /* Local variables */
    double d__;
    int i__, j, k;
    double p, q, u[12], r12, r21, r22, et[2], du[12];
    double xi[2], zz, dd1, dd2, dd3, dua[12], dub[12], duc[12];
    int ket[2], kxi[2];
    double ddip;
    double aalpha;
	Tipo_c2_1 c2_1;
	union {
		Tipo_1 _1;
		Tipo_2 _2;
	} c0_;


/* ******************************************************************** */
/* *****                                                          ***** */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH                      ***** */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
/* *****              CODED BY  Y.OKADA ... SEP.1991              ***** */
/* *****              REVISED ... NOV.1991, APR.1992, MAY.1993,   ***** */
/* *****                          JUL.1993, MAY.2002              ***** */
/* ******************************************************************** */

/* ***** INPUT */
/* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU) */
/* *****   X,Y,Z : COORDINATE OF OBSERVING POINT */
/* *****   DEPTH : DEPTH OF REFERENCE POINT */
/* *****   DIP   : DIP-ANGLE (DEGREE) */
/* *****   AL1,AL2   : FAULT LENGTH RANGE */
/* *****   AW1,AW2   : FAULT WIDTH RANGE */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */

/* ***** OUTPUT */
/* *****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL) */
/* *****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) / */
/* *****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) ) */
/* *****   UXZ,UYZ,UZZ : Z-DERIVATIVE */
/* *****   IRET        : RETURN CODE */
/* *****               :   =0....NORMAL */
/* *****               :   =1....SINGULAR */
/* *****               :   =2....POSITIVE Z WAS GIVEN */

/* ----- */
    *iret = 0;
    if (*z__ > 0.0) {
	*iret = 2;
	goto L99;
    }
/* ----- */
    for (i__ = 1; i__ <= 12; ++i__) {
	u[i__ - 1] = f0;
	dua[i__ - 1] = f0;
	dub[i__ - 1] = f0;
	duc[i__ - 1] = f0;
/* L111: */
    }
    aalpha = *alpha;
    ddip = *dip;
    dccon0_(&aalpha, &ddip, &c0_._2);
/* ----- */
    zz = *z__;
    dd1 = *disl1;
    dd2 = *disl2;
    dd3 = *disl3;
    xi[0] = *x - *al1;
    xi[1] = *x - *al2;
    if (abs(xi[0]) < eps) {
	xi[0] = f0;
    }
    if (abs(xi[1]) < eps) {
	xi[1] = f0;
    }
/* ====================================== */
/* =====  REAL-SOURCE CONTRIBUTION  ===== */
/* ====================================== */
    d__ = *depth + *z__;
    p = *y * c0_._1.cd + d__ * c0_._1.sd;
    q = *y * c0_._1.sd - d__ * c0_._1.cd;
    et[0] = p - *aw1;
    et[1] = p - *aw2;
    if (abs(q) < eps) {
	q = f0;
    }
    if (abs(et[0]) < eps) {
	et[0] = f0;
    }
    if (abs(et[1]) < eps) {
	et[1] = f0;
    }
/* -------------------------------- */
/* ----- REJECT SINGULAR CASE ----- */
/* -------------------------------- */
/* ----- ON FAULT EDGE */
    if (q == f0 && (xi[0] * xi[1] <= f0 && et[0] * et[1] == f0 || et[0] * et[
	    1] <= f0 && xi[0] * xi[1] == f0)) {
	*iret = 1;
	goto L99;
    }
/* ----- ON NEGATIVE EXTENSION OF FAULT EDGE */
    kxi[0] = 0;
    kxi[1] = 0;
    ket[0] = 0;
    ket[1] = 0;
    r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
    r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
    r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
    if (xi[0] < f0 && r21 + xi[1] < eps) {
	kxi[0] = 1;
    }
    if (xi[0] < f0 && r22 + xi[1] < eps) {
	kxi[1] = 1;
    }
    if (et[0] < f0 && r12 + et[1] < eps) {
	ket[0] = 1;
    }
    if (et[0] < f0 && r22 + et[1] < eps) {
	ket[1] = 1;
    }
/* ===== */
    for (k = 1; k <= 2; ++k) {
	for (j = 1; j <= 2; ++j) {
	    dccon2_(&xi[j - 1], &et[k - 1], &q, &c0_._1.sd, &c0_._1.cd, &kxi[k - 
		    1], &ket[j - 1], &c2_1);
	    ua_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dua, &c0_._2, &c2_1);
/* ----- */
	    for (i__ = 1; i__ <= 10; i__ += 3) {
		du[i__ - 1] = -dua[i__ - 1];
		du[i__] = -dua[i__] * c0_._1.cd + dua[i__ + 1] * c0_._1.sd;
		du[i__ + 1] = -dua[i__] * c0_._1.sd - dua[i__ + 1] * c0_._1.cd;
		if (i__ < 10) {
		    goto L220;
		}
		du[i__ - 1] = -du[i__ - 1];
		du[i__] = -du[i__];
		du[i__ + 1] = -du[i__ + 1];
L220:
		;
	    }
	    for (i__ = 1; i__ <= 12; ++i__) {
		if (j + k != 3) {
		    u[i__ - 1] += du[i__ - 1];
		}
		if (j + k == 3) {
		    u[i__ - 1] -= du[i__ - 1];
		}
/* L221: */
	    }
/* ----- */
/* L222: */
	}
/* L223: */
    }
/* ======================================= */
/* =====  IMAGE-SOURCE CONTRIBUTION  ===== */
/* ======================================= */
    d__ = *depth - *z__;
    p = *y * c0_._1.cd + d__ * c0_._1.sd;
    q = *y * c0_._1.sd - d__ * c0_._1.cd;
    et[0] = p - *aw1;
    et[1] = p - *aw2;
    if (abs(q) < eps) {
	q = f0;
    }
    if (abs(et[0]) < eps) {
	et[0] = f0;
    }
    if (abs(et[1]) < eps) {
	et[1] = f0;
    }
/* -------------------------------- */
/* ----- REJECT SINGULAR CASE ----- */
/* -------------------------------- */
/* ----- ON FAULT EDGE */
    if (q == f0 && (xi[0] * xi[1] <= f0 && et[0] * et[1] == f0 || et[0] * et[
	    1] <= f0 && xi[0] * xi[1] == f0)) {
	*iret = 1;
	goto L99;
    }
/* ----- ON NEGATIVE EXTENSION OF FAULT EDGE */
    kxi[0] = 0;
    kxi[1] = 0;
    ket[0] = 0;
    ket[1] = 0;
    r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
    r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
    r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
    if (xi[0] < f0 && r21 + xi[1] < eps) {
	kxi[0] = 1;
    }
    if (xi[0] < f0 && r22 + xi[1] < eps) {
	kxi[1] = 1;
    }
    if (et[0] < f0 && r12 + et[1] < eps) {
	ket[0] = 1;
    }
    if (et[0] < f0 && r22 + et[1] < eps) {
	ket[1] = 1;
    }
/* ===== */
    for (k = 1; k <= 2; ++k) {
	for (j = 1; j <= 2; ++j) {
	    dccon2_(&xi[j - 1], &et[k - 1], &q, &c0_._1.sd, &c0_._1.cd, &kxi[k - 
		    1], &ket[j - 1], &c2_1);
	    ua_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dua, &c0_._2, &c2_1);
	    ub_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dub, &c0_._2, &c2_1);
	    uc_(&xi[j - 1], &et[k - 1], &q, &zz, &dd1, &dd2, &dd3, duc, &c0_._2, &c2_1);
/* ----- */
	    for (i__ = 1; i__ <= 10; i__ += 3) {
		du[i__ - 1] = dua[i__ - 1] + dub[i__ - 1] + *z__ * duc[i__ - 
			1];
		du[i__] = (dua[i__] + dub[i__] + *z__ * duc[i__]) * c0_._1.cd - 
			(dua[i__ + 1] + dub[i__ + 1] + *z__ * duc[i__ + 1]) * 
			c0_._1.sd;
		du[i__ + 1] = (dua[i__] + dub[i__] - *z__ * duc[i__]) * 
			c0_._1.sd + (dua[i__ + 1] + dub[i__ + 1] - *z__ * duc[
			i__ + 1]) * c0_._1.cd;
		if (i__ < 10) {
		    goto L330;
		}
		du[9] += duc[0];
		du[10] = du[10] + duc[1] * c0_._1.cd - duc[2] * c0_._1.sd;
		du[11] = du[11] - duc[1] * c0_._1.sd - duc[2] * c0_._1.cd;
L330:
		;
	    }
	    for (i__ = 1; i__ <= 12; ++i__) {
		if (j + k != 3) {
		    u[i__ - 1] += du[i__ - 1];
		}
		if (j + k == 3) {
		    u[i__ - 1] -= du[i__ - 1];
		}
/* L331: */
	    }
/* ----- */
/* L333: */
	}
/* L334: */
    }
/* ===== */
    *ux = u[0];
    *uy = u[1];
    *uz = u[2];
    *uxx = u[3];
    *uyx = u[4];
    *uzx = u[5];
    *uxy = u[6];
    *uyy = u[7];
    *uzy = u[8];
    *uxz = u[9];
    *uyz = u[10];
    *uzz = u[11];
    return 0;
/* =========================================== */
/* =====  IN CASE OF SINGULAR (ON EDGE)  ===== */
/* =========================================== */
L99:
    *ux = f0;
    *uy = f0;
    *uz = f0;
    *uxx = f0;
    *uyx = f0;
    *uzx = f0;
    *uxy = f0;
    *uyy = f0;
    *uzy = f0;
    *uxz = f0;
    *uyz = f0;
    *uzz = f0;
    return 0;
} /* dc3d_ */

#endif
