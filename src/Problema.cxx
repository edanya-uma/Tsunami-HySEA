#include "Constantes.hxx"
#include <sys/stat.h> 
#include <fstream>
#include <sstream>
#include <cstring>
#include <limits.h>

// NetCDF external functions
extern void abrirGRD(const char *nombre_fich, int *nvx, int *nvy);
extern void leerLongitudGRD(double *lon);
extern void leerLatitudGRD(double *lat);
extern void leerBatimetriaGRD(float *bati);
extern void leerEtaGRD(float *eta);
extern void cerrarGRD();

bool existeFichero(string fichero) 
{
	struct stat stFichInfo;
	bool existe;
	int intStat;

	intStat = stat(fichero.c_str(), &stFichInfo);
	if (intStat == 0) {
		existe = true;
	}
	else {
		existe = false;
	}

	return existe;
}

void liberarMemoria(int numNiveles, double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2,
		tipoDatosSubmalla datosNivel, int *posicionesVolumenesGuardado, double *lonPuntos,
		double *latPuntos)
{
	if (datosNivel.areaYCosPhi != NULL)			free(datosNivel.areaYCosPhi);
	if (datosNivel.anchoVolumenes != NULL)		free(datosNivel.anchoVolumenes);
	if (datosNivel.altoVolumenes != NULL)		free(datosNivel.altoVolumenes);
	if (datosNivel.longitud != NULL)			free(datosNivel.longitud);
	if (datosNivel.latitud != NULL)				free(datosNivel.latitud);
	if (datosVolumenesNivel_1 != NULL)			free(datosVolumenesNivel_1);
	if (datosVolumenesNivel_2 != NULL)			free(datosVolumenesNivel_2);
	if (posicionesVolumenesGuardado != NULL)	free(posicionesVolumenesGuardado);
	if (lonPuntos != NULL)						free(lonPuntos);
	if (latPuntos != NULL)						free(latPuntos);
}

int obtenerPosEnPuntosGuardado(double *lonPuntos, double *latPuntos, int num_puntos, double lon, double lat)
{
	int i;
	bool encontrado = false;

	i = 0;
	while ((i < num_puntos) && (! encontrado)) {
		if ((fabs(lonPuntos[i] - lon) < 1e-6) && (fabs(latPuntos[i] - lat) < 1e-6))
			encontrado = true;
		else i++;
	}

	return (encontrado ? i : -1);
}

int obtenerIndicePunto(int numVolxNivel0, int numVolyNivel0, tipoDatosSubmalla *datosNivel, double lon, double lat)
{
	int i, j;
	int pos;
	double *longitud, *latitud;
	bool encontrado;

	pos = -1;
	encontrado = false;
	longitud = datosNivel->longitud;
	latitud = datosNivel->latitud;
	if (lon >= longitud[0]) {
		i = 0;
		while ((i < numVolxNivel0) && (! encontrado))  {
			if (longitud[i] >= lon) {
				encontrado = true;
				if (fabs(lon - longitud[max(i-1,0)]) < fabs(lon - longitud[i]))
					i = max(i-1,0);
			}
			else i++;
		}

		if (encontrado) {
			encontrado = false;
			if (lat >= latitud[0]) {
				j = 0;
				while ((j < numVolyNivel0) && (! encontrado))  {
					if (latitud[j] >= lat) {
						encontrado = true;
						if (fabs(lat - latitud[max(j-1,0)]) < fabs(lat - latitud[j]))
							j = max(j-1,0);
					}
					else j++;
				}
				if (encontrado)
					pos = j*numVolxNivel0 + i;
			}
		}
	}

	return pos;
}

void obtenerDatosFallaOkadaStandard(ifstream &fich, double &LON_C, double &LAT_C, double &DEPTH_C, double &FAULT_L,
		double &FAULT_W, double &STRIKE, double &DIP, double &RAKE, double &SLIP)
{
	string linea;
	istringstream iss;
	size_t pos;
	bool repetir = true;

	while (repetir) {
		getline(fich,linea);
		pos = linea.find_first_not_of(" \t\r\n");
		if (pos != string::npos) {
			linea.erase(0,pos);
			if (linea[0] != '#')
				repetir = false;
		}
	}
	iss.str(linea);
	iss >> LON_C >> LAT_C >> DEPTH_C >> FAULT_L >> FAULT_W >> STRIKE >> DIP >> RAKE >> SLIP;
}

template <class T>
void obtenerSiguienteDato(ifstream &fich, T &dato)
{
 	string linea;
	istringstream iss;
	size_t pos;
	bool repetir = true;

	while (repetir) {
		getline(fich,linea);
		pos = linea.find_first_not_of(" \t\r\n");
		if (pos != string::npos) {
			linea.erase(0,pos);
			if (linea[0] != '#')
				repetir = false;
		}
	}
	iss.str(linea);
	iss >> dato;
}

template void obtenerSiguienteDato<int>(ifstream &fich, int &dato);
template void obtenerSiguienteDato<double>(ifstream &fich, double &dato);
template void obtenerSiguienteDato<string>(ifstream &fich, string &dato);

int cargarDatosProblema(string fich_ent, string &nombre_bati, string &prefijo, int *numNiveles, int *okada_flag,
		double *LON_C, double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE,
		double *DIP, double *RAKE, double *SLIP, double2 **datosVolumenesNivel_1, double2 **datosVolumenesNivel_2,
		tipoDatosSubmalla *datosNivel, int *numVolxNivel0, int *numVolyNivel0, int64_t *numVolumenesNivel,
		int *leer_fichero_puntos, int *numPuntosGuardar, int **posicionesVolumenesGuardado, double **lonPuntos,
		double **latPuntos, double *Hmin, double *borde_sup, double *borde_inf, double *borde_izq, double *borde_der,
		int *tam_spongeSup, int *tam_spongeInf, int *tam_spongeIzq, int *tam_spongeDer, double *tiempo_tot,
		double *tiempoGuardarNetCDF, double *tiempoGuardarSeries, double *CFL, double *mf0, double *vmax,
		double *epsilon_h, double *hpos, double *cvis, double *L, double *H, double *Q, double *T)
{
	int i, j;
	int pos, posPunto;
	double val, delta_lon, delta_lat;
	double radio_tierra = EARTH_RADIUS;
	double grad2rad = M_PI/180.0;
	ifstream fich2;
	ifstream fich_ptos;
	string directorio;
	string fich_topo, fich_est;
	string fich_puntos;
	double lon, lat;
	double *datosGRD;
	float *datosGRD_float;
	int64_t tam_datosGRD, tam;

	i = fich_ent.find_last_of("/");
	if (i > -1) {
		directorio = fich_ent.substr(0,i)+"/";
	}
	else {
		i = fich_ent.find_last_of("\\");
		if (i > -1) {
			directorio = fich_ent.substr(0,i)+"\\";
		}
		else {
			directorio = "";
		}
	}

	*tiempoGuardarSeries = -1.0;
	ifstream fich(fich_ent.c_str());
	obtenerSiguienteDato<string>(fich, nombre_bati);
	obtenerSiguienteDato<string>(fich, fich_topo);
	if (! existeFichero(directorio+fich_topo)) {
		cerr << "Error: File '" << directorio+fich_topo << "' not found" << endl;
		fich.close();
		return 1;
	}
	obtenerSiguienteDato<int>(fich, *okada_flag);
	if ((*okada_flag != SEA_SURFACE_FROM_FILE) && (*okada_flag != OKADA_STANDARD)) {
		cerr << "Error: The initialization flag should be " << SEA_SURFACE_FROM_FILE << " or " << OKADA_STANDARD << endl;
		fich.close();
		return 1;
	}
	if (*okada_flag == SEA_SURFACE_FROM_FILE) {
		obtenerSiguienteDato<string>(fich, fich_est);
		if (! existeFichero(directorio+fich_est)) {
			cerr << "Error: File '" << directorio+fich_est << "' not found" << endl;
			fich.close();
			return 1;
		}
	}
	else if (*okada_flag == OKADA_STANDARD) {
		obtenerDatosFallaOkadaStandard(fich, *LON_C, *LAT_C, *DEPTH_C, *FAULT_L, *FAULT_W, *STRIKE, *DIP, *RAKE, *SLIP);
	}
	obtenerSiguienteDato<string>(fich, prefijo);
	abrirGRD((directorio+fich_topo).c_str(), numVolxNivel0, numVolyNivel0);
	if ((*numVolxNivel0 < 2) || (*numVolyNivel0 < 2)) {
		cerr << "Error: Mesh size too small. The number of rows and columns should be >= 2" << endl;
		fich.close();
		cerrarGRD();
		return 1;
	}


	obtenerSiguienteDato<int>(fich, *numNiveles);
	if (*numNiveles != 1) {
		cerr << "Error: The number of levels should be 1" << endl;
		fich.close();
		cerrarGRD();
		return 1;
	}
	obtenerSiguienteDato<double>(fich, *borde_sup);
	obtenerSiguienteDato<double>(fich, *borde_inf);
	obtenerSiguienteDato<double>(fich, *borde_izq);
	obtenerSiguienteDato<double>(fich, *borde_der);
	obtenerSiguienteDato<double>(fich, *tiempo_tot);
	obtenerSiguienteDato<double>(fich, *tiempoGuardarNetCDF);
	obtenerSiguienteDato<int>(fich, *leer_fichero_puntos);
	if (*leer_fichero_puntos == 1) {
		obtenerSiguienteDato<string>(fich, fich_puntos);
		obtenerSiguienteDato<double>(fich, *tiempoGuardarSeries);
		if (! existeFichero(directorio+fich_puntos)) {
			cerr << "Error: File '" << directorio+fich_puntos << "' not found" << endl;
			fich.close();
			cerrarGRD();
			return 1;
		}
	}
	obtenerSiguienteDato<double>(fich, *CFL);
	obtenerSiguienteDato<double>(fich, *epsilon_h);
	obtenerSiguienteDato<double>(fich, *hpos);
	obtenerSiguienteDato<double>(fich, *cvis);
	obtenerSiguienteDato<double>(fich, *mf0);
	obtenerSiguienteDato<double>(fich, *vmax);
	obtenerSiguienteDato<double>(fich, *L);
	obtenerSiguienteDato<double>(fich, *H);
	*Q = sqrt(9.81*pow((*H),3.0));
	*T = (*L)*(*H)/(*Q);
	*tiempo_tot /= *T;
	*tiempoGuardarNetCDF /= *T;
	*tiempoGuardarSeries /= *T;
	*mf0 *= 9.81*(*mf0)*(*L)/pow((*H),4.0/3.0);
	*vmax /= (*Q)/(*H);
	*epsilon_h /= *H;
	*hpos /= *H;
	*cvis = 1.0 - (*cvis);
	radio_tierra /= *L;
	fich.close();

	*tam_spongeIzq = ((fabs(*borde_izq-1.0) < EPSILON) ? SPONGE_SIZE : 0);
	*tam_spongeDer = ((fabs(*borde_der-1.0) < EPSILON) ? SPONGE_SIZE : 0);
	*tam_spongeSup = ((fabs(*borde_sup-1.0) < EPSILON) ? SPONGE_SIZE : 0);
	*tam_spongeInf = ((fabs(*borde_inf-1.0) < EPSILON) ? SPONGE_SIZE : 0);

	tam = 0;
	if (*leer_fichero_puntos == 1) {
		fich_ptos.open((directorio+fich_puntos).c_str());
		fich_ptos >> tam;
		fich_ptos.close();
	}

	*numVolumenesNivel = (int64_t) ((*numVolxNivel0)*(*numVolyNivel0));
	tam = max(tam, *numVolumenesNivel);
	tam_datosGRD = (*numVolumenesNivel)*sizeof(float);
	datosGRD_float = (float *) malloc(tam_datosGRD);
	datosGRD = (double *) datosGRD_float;
	datosNivel->areaYCosPhi = (double2 *) malloc((*numVolyNivel0)*sizeof(double2));
	datosNivel->anchoVolumenes = (double *) malloc(((*numVolyNivel0)+1)*sizeof(double));
	datosNivel->altoVolumenes  = (double *) malloc((*numVolyNivel0)*sizeof(double));
	datosNivel->longitud = (double *) malloc((*numVolxNivel0)*sizeof(double));
	datosNivel->latitud  = (double *) malloc((*numVolyNivel0)*sizeof(double));
	*datosVolumenesNivel_1 = (double2 *) malloc(tam*((int64_t) sizeof(double2)));
	*datosVolumenesNivel_2 = (double2 *) malloc(tam*((int64_t) sizeof(double2)));
	if (*datosVolumenesNivel_2 == NULL) {
		cerr << "Error: Not enough CPU memory" << endl;
		liberarMemoria(*numNiveles, *datosVolumenesNivel_1, *datosVolumenesNivel_2, *datosNivel,
			*posicionesVolumenesGuardado, *lonPuntos, *latPuntos);
		free(datosGRD);
		cerrarGRD();
		return 1;
	}

	leerLongitudGRD(datosNivel->longitud);
	leerLatitudGRD(datosNivel->latitud);
	leerBatimetriaGRD(datosGRD_float);
	cerrarGRD();
	*Hmin = 1e30;
	for (j=(*numVolyNivel0)-1; j>=0; j--) {
		for (i=0; i<(*numVolxNivel0); i++) {
			pos = j*(*numVolxNivel0) + i;
			val = -1.0*datosGRD_float[pos];
			val /= *H;
			(*datosVolumenesNivel_1)[pos].y = val;
			if (*okada_flag == OKADA_STANDARD) {
				(*datosVolumenesNivel_1)[pos].x = ((val > 0.0) ? val : 0.0);
				(*datosVolumenesNivel_2)[pos].x = 0.0;
				(*datosVolumenesNivel_2)[pos].y = 0.0;
			}
			if (val < *Hmin)
				*Hmin = val;
		}
	}

	*numPuntosGuardar = 0;
	if (*leer_fichero_puntos == 1) {
		fich_ptos.open((directorio+fich_puntos).c_str());
		fich_ptos >> j;
		*lonPuntos = (double *) malloc(j*sizeof(double));
		*latPuntos = (double *) malloc(j*sizeof(double));
		*posicionesVolumenesGuardado = (int *) malloc(j*sizeof(int));
		if (*posicionesVolumenesGuardado == NULL) {
			cerr << "Error: Not enough CPU memory" << endl;
			liberarMemoria(*numNiveles, *datosVolumenesNivel_1, *datosVolumenesNivel_2, *datosNivel,
				*posicionesVolumenesGuardado, *lonPuntos, *latPuntos);
			free(datosGRD);
			fich_ptos.close();
			return 1;
		}

		for (i=0; i<j; i++) {
			fich_ptos >> lon;
			fich_ptos >> lat;
			pos = obtenerPosEnPuntosGuardado(*lonPuntos, *latPuntos, *numPuntosGuardar, lon, lat);
			if (pos == -1) {
				(*lonPuntos)[*numPuntosGuardar] = lon;
				(*latPuntos)[*numPuntosGuardar] = lat;
				posPunto = obtenerIndicePunto(*numVolxNivel0, *numVolyNivel0, datosNivel, lon, lat);
				(*posicionesVolumenesGuardado)[*numPuntosGuardar] = posPunto;
				(*numPuntosGuardar)++;
			}
		}
		fich_ptos.close();
	}

	delta_lon = fabs(datosNivel->longitud[1] - datosNivel->longitud[0]);
	for (i=0; i<(*numVolyNivel0)-1; i++) {
		delta_lat = fabs(datosNivel->latitud[i+1] - datosNivel->latitud[i]);
		datosNivel->anchoVolumenes[i] = delta_lon;
		datosNivel->altoVolumenes[i] = delta_lat;
	}
	datosNivel->anchoVolumenes[i] = delta_lon;
	datosNivel->altoVolumenes[i] = datosNivel->altoVolumenes[i-1];
	datosNivel->anchoVolumenes[i+1] = delta_lon;

	delta_lon = fabs(datosNivel->longitud[1] - datosNivel->longitud[0]);
	for (i=0; i<(*numVolyNivel0)-1; i++) {
		delta_lat = fabs(datosNivel->latitud[i+1] - datosNivel->latitud[i]);
		datosNivel->areaYCosPhi[i].x = radio_tierra*delta_lon*delta_lat*grad2rad;
		datosNivel->areaYCosPhi[i].y = cos(datosNivel->latitud[i]*grad2rad);
	}
	datosNivel->areaYCosPhi[i].x = radio_tierra*delta_lon*delta_lat*grad2rad;
	datosNivel->areaYCosPhi[i].y = cos(datosNivel->latitud[i]*grad2rad);

	if (*Hmin < 0.0) {
		for (i=0; i<(*numVolumenesNivel); i++)
			(*datosVolumenesNivel_1)[i].y -= *Hmin;
	}
	else {
		*Hmin = 0.0;
	}

	if (*okada_flag == SEA_SURFACE_FROM_FILE) {
		abrirGRD((directorio+fich_est).c_str(), &i, &j);
		leerEtaGRD(datosGRD_float);
		cerrarGRD();
		for (j=(*numVolyNivel0)-1; j>=0; j--) {
			for (i=0; i<(*numVolxNivel0); i++) {
				pos = j*(*numVolxNivel0) + i;
				val = (double) datosGRD_float[pos];
				val = val/(*H) + (*datosVolumenesNivel_1)[pos].y + (*Hmin);
				val = max(val, 0.0);                                
				(*datosVolumenesNivel_1)[pos].x = val;
				(*datosVolumenesNivel_2)[pos].x = 0.0;
				(*datosVolumenesNivel_2)[pos].y = 0.0;
			}
		}
	}
	free(datosGRD);

	return 0;
}

void mostrarDatosProblema(int numNiveles, int okada_flag, tipoDatosSubmalla datosNivel, int numVolxNivel0, int numVolyNivel0,
				double tiempo_tot, double tiempoGuardarNetCDF, int leer_fichero_puntos, double tiempoGuardarSeries, double CFL,
				double mf0, double vmax, double epsilon_h, double hpos, double cvis, double L, double H, double Q, double T)
{
	cout << "/**********************************************************************" << endl;
	cout << " Tsunami-HySEA numerical model open source v1.1.1                      " << endl;
	cout << " developed by the EDANYA Research Group, University of Malaga (Spain). " << endl;
	cout << " https://www.uma.es/edanya                                             " << endl;
	cout << " https://edanya.uma.es/hysea/                                          " << endl;
	cout << " Contact and support: hysea@uma.es                                     " << endl;
	cout << " The core version of Tsunami-HySEA is distributed under                " << endl;
	cout << " GNU General Public License, Version 2                                 " << endl;
	cout << "**********************************************************************/" << endl;
	cout << endl;

	cout << "Problem data" << endl;
	if (okada_flag == SEA_SURFACE_FROM_FILE) {
		cout << "Initialization: Sea surface displacement from file" << endl;
	}
	else if (okada_flag == OKADA_STANDARD) {
		cout << "Initialization: Standard Okada" << endl;
	}
	cout << "CFL: " << CFL << endl;
	cout << "Water-bottom friction: " << sqrt(mf0*pow(H,4.0/3.0)/(9.81*L)) << endl;
	cout << "Maximum allowed velocity of water: " << vmax*(Q/H) << endl;
	cout << "Epsilon h: " << epsilon_h*H << " m" << endl;
	cout << "Threshold for the 2s+WAF scheme: " << hpos*H << " m" << endl;
	cout << "Stability coefficient: " << 1.0-cvis << endl;
	cout << "Simulation time: " << tiempo_tot*T << " sec" << endl;
	cout << "Saving time of NetCDF files: " << tiempoGuardarNetCDF*T << " sec" << endl;
	if (leer_fichero_puntos) {
		cout << "Time series: yes (saving time: " << tiempoGuardarSeries*T << " sec)" << endl;
	}
	else {
		cout << "Time series: no" << endl;
	}
	cout << "Number of levels: " << numNiveles << endl;
	cout << "Level 0" << endl;
	cout << "  Volumes: " << numVolxNivel0 << " x " << numVolyNivel0 << " = " << numVolxNivel0*numVolyNivel0 << endl;
	cout << "  Longitude: [" << datosNivel.longitud[0] << ", " << datosNivel.longitud[numVolxNivel0-1] << "]" << endl;
	cout << "  Latitude: [" << datosNivel.latitud[0] << ", " << datosNivel.latitud[numVolyNivel0-1] << "]" << endl;
	cout << endl;
}

