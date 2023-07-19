#ifndef PRTOXY_H
#define PRTOXY_H

__device__ /* Subroutine */ int prtoxy_(double *alatdg, double *alngdg, 
	double *alato, double *alngo, double *x, double *y, 
	int *ind)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    double r__, c1, c2, v2, al, bl, an, ph1, den, cph1, rph1, rph2,
	     tph1, clat, rlat, slat, srph1, srph2, clato, rlato, slato;

/*      C*************************************************************C */
/*      C*****                                                   *****C */
/*      C***** Conversion between (lat,lon) and (X,Y)            *****C */
/*      C*****      using Gauss-Krueger projection               *****C */
/*      C*****                                                   *****C */
/*      C*************************************************************C */
/*      C***** Input/Output */
/*      C***** ALATDG,ALNGDG : latitude, longitude (deg) */
/*      C***** X , Y : EW, NS coordinates (km) */
/*      C***** Input */
/*      C***** ALATO,ALNGO : origin of projection (deg) */
/*      C***** IND : =0 ... transform (lat,lng) to ( X , Y ) */
/*      C***** : =1 ... transform ( X , Y ) to (lat,lng) */
/*      C------------------------------------------------ */
/*      C----- IND=0 : transform (lat,lng) to (X,Y) ----- */
/*      C------------------------------------------------ */
    if (*ind == 0) {
	rlat = *alatdg * 0.017453292371619689;
	slat = sin(rlat);
	clat = cos(rlat);
	v2 = clat * 0.0067397251 * clat + 1.0;
	al = *alngdg - *alngo;
	ph1 = *alatdg + v2 * 0.5 * al * al * slat * clat * 
		0.017453292371619689;
	rph1 = ph1 * 0.017453292371619689;
	rph2 = (ph1 + *alato) * 0.5 * 0.017453292371619689;
	srph1 = sin(rph1);
	srph2 = sin(rph2);
/*      C----- */
/* Computing 3rd power */
	d__1 = sqrt(1.0 - srph2 * 0.0066946053 * srph2);
	r__ = 6335.4607362597517 / (d__1 * (d__1 * d__1));
	an = 6378.16 / sqrt(1.0 - srph1 * 0.0066946053 * srph1);
	c1 = 57.29578 / r__;
	c2 = 57.29578 / an;
	*y = (ph1 - *alato) / c1;
	*x = al * clat / c2 * (al * al * cos(rlat * 2.0) / 19696.838434850401 
		+ 1.0);
/*      C------------------------------------------------ */
/*      C----- IND=1 : transform (X,Y) to (LAT,LNG) ----- */
/*      C------------------------------------------------ */
    } else if (*ind == 1) {
	rlato = *alato * 0.017453292371619689;
	slato = sin(rlato);
	clato = cos(rlato);
	den = sqrt(1.0 - slato * 0.0066946053 * slato);
/* Computing 3rd power */
	d__1 = den;
	r__ = 6335.4607362597517 / (d__1 * (d__1 * d__1));
	an = 6378.16 / den;
	v2 = clato * 0.0067397251 * clato + 1.0;
/*      C----- */
	c1 = 57.29578 / r__;
	c2 = 57.29578 / an;
	ph1 = *alato + c1 * *y;
	rph1 = ph1 * 0.017453292371619689;
	tph1 = tan(rph1);
	cph1 = cos(rph1);
	bl = c2 * *x;
	*alatdg = ph1 - bl * 0.5 * bl * v2 * tph1 * 0.017453292371619689;
	*alngdg = *alngo + bl / cph1 * (1.0 - bl * bl * (tph1 * 2.0 * tph1 + 
		1.0) / 19696.838434850401);
    }
    return 0;
} /* prtoxy_ */

#endif
