#include "Constantes.hxx"
#include <sys/stat.h> 
#include <fstream>
#include <sstream>
#include <cstring>
#include <limits.h>

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
		tipoDatosSubmalla datosNivel)
{
	if (datosNivel.areaYCosPhi != NULL)		free(datosNivel.areaYCosPhi);
	if (datosNivel.anchoVolumenes != NULL)	free(datosNivel.anchoVolumenes);
	if (datosNivel.altoVolumenes != NULL)	free(datosNivel.altoVolumenes);
	if (datosNivel.longitud != NULL)		free(datosNivel.longitud);
	if (datosNivel.latitud != NULL)			free(datosNivel.latitud);
	if (datosVolumenesNivel_1 != NULL)		free(datosVolumenesNivel_1);
	if (datosVolumenesNivel_2 != NULL)		free(datosVolumenesNivel_2);
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

int cargarDatosProblema(string fich_ent, string &nombre_bati, string &prefijo, int *numNiveles,
			double2 **datosVolumenesNivel_1, double2 **datosVolumenesNivel_2, tipoDatosSubmalla *datosNivel,
			int *numVolxNivel0, int *numVolyNivel0, int *numVolumenesNivel, int *leer_fichero_puntos,
			double *Hmin, double *borde_sup, double *borde_inf, double *borde_izq, double *borde_der,
			double *tiempo_tot, double *tiempoGuardarNetCDF, double *CFL, double *mf0, double *vmax,
			double *epsilon_h, double *hpos, double *cvis, double *L, double *H, double *Q, double *T)
{
	int i, j, pos;
	double val, delta_lon, delta_lat;
	double radio_tierra = EARTH_RADIUS;
	double grad2rad = M_PI/180.0;
	ifstream fich2;
	string directorio;
	string fich_topo, fich_est, fich_est2;
	int sized = sizeof(double);
	double lon, lat;

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

	ifstream fich(fich_ent.c_str());
	obtenerSiguienteDato<string>(fich, nombre_bati);
	obtenerSiguienteDato<int>(fich, *numVolxNivel0);
	obtenerSiguienteDato<int>(fich, *numVolyNivel0);
	if ((*numVolxNivel0 < 2) || (*numVolyNivel0 < 2)) {
		cerr << "Error: Mesh size too small. The number of rows and columns should be >= 2" << endl;
		fich.close();
		return 1;
	}
	obtenerSiguienteDato<string>(fich, fich_topo);
	obtenerSiguienteDato<string>(fich, fich_est);
	obtenerSiguienteDato<string>(fich, prefijo);
	if (! existeFichero(directorio+fich_topo)) {
		cerr << "Error: File '" << directorio+fich_topo << "' not found" << endl;
		fich.close();
		return 1;
	}
	if (! existeFichero(directorio+fich_est)) {
		cerr << "Error: File '" << directorio+fich_est << "' not found" << endl;
		fich.close();
		return 1;
	}
	fich2.open((directorio+fich_topo).c_str(), ios::in | ios::binary);

	obtenerSiguienteDato<int>(fich, *numNiveles);
	if (*numNiveles != 1) {
		cerr << "Error: The number of levels should be 1" << endl;
		fich.close();
		fich2.close();
		return 1;
	}
	obtenerSiguienteDato<double>(fich, *borde_sup);
	obtenerSiguienteDato<double>(fich, *borde_inf);
	obtenerSiguienteDato<double>(fich, *borde_izq);
	obtenerSiguienteDato<double>(fich, *borde_der);
	obtenerSiguienteDato<double>(fich, *tiempo_tot);
	obtenerSiguienteDato<double>(fich, *tiempoGuardarNetCDF);
	obtenerSiguienteDato<int>(fich, *leer_fichero_puntos);
	if (*leer_fichero_puntos != 0) {
		cerr << "Error: Saving of time series is not supported" << endl;
		fich.close();
		fich2.close();
		return 1;
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
	*mf0 *= 9.81*(*mf0)*(*L)/pow((*H),4.0/3.0);
	*vmax /= (*Q)/(*H);
	*epsilon_h /= *H;
	*hpos /= *H;
	*cvis = 1.0 - (*cvis);
	radio_tierra /= *L;
	fich.close();

	*numVolumenesNivel = (*numVolxNivel0)*(*numVolyNivel0);
	datosNivel->areaYCosPhi = (double2 *) malloc((*numVolyNivel0)*sizeof(double2));
	datosNivel->anchoVolumenes = (double *) malloc(((*numVolyNivel0)+1)*sizeof(double));
	datosNivel->altoVolumenes  = (double *) malloc((*numVolyNivel0)*sizeof(double));
	datosNivel->longitud = (double *) malloc((*numVolxNivel0)*sizeof(double));
	datosNivel->latitud  = (double *) malloc((*numVolyNivel0)*sizeof(double));
	*datosVolumenesNivel_1 = (double2 *) malloc((*numVolumenesNivel)*sizeof(double2));
	*datosVolumenesNivel_2 = (double2 *) malloc((*numVolumenesNivel)*sizeof(double2));
	if (*datosVolumenesNivel_2 == NULL) {
		cerr << "Error: Not enough CPU memory" << endl;
		liberarMemoria(*numNiveles, *datosVolumenesNivel_1, *datosVolumenesNivel_2, *datosNivel);
		fich2.close();
		return 1;
	}

	*Hmin = 1e30;
	for (j=(*numVolyNivel0)-1; j>=0; j--) {
		for (i=0; i<(*numVolxNivel0); i++) {
			fich2.read((char *)(datosNivel->longitud+i), sized);
			fich2.read((char *)(datosNivel->latitud+j), sized);
			fich2.read((char *)&val, sized);
			val /= *H;
			pos = j*(*numVolxNivel0) + i;
			(*datosVolumenesNivel_1)[pos].y = val;
			if (val < *Hmin)
				*Hmin = val;
		}
	}
	fich2.close();

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

	fich2.open((directorio+fich_est).c_str(), ios::in | ios::binary);
	for (j=(*numVolyNivel0)-1; j>=0; j--) {
		for (i=0; i<(*numVolxNivel0); i++) {
			pos = j*(*numVolxNivel0) + i;
			fich2.read((char *)&val, sized);
			val = val/(*H) + (*datosVolumenesNivel_1)[pos].y + (*Hmin);
			val = max(val, 0.0);
			(*datosVolumenesNivel_1)[pos].x = val;
			fich2.read((char *)&val, sized);  (*datosVolumenesNivel_2)[pos].x = val/(*Q);
			fich2.read((char *)&val, sized);  (*datosVolumenesNivel_2)[pos].y = val/(*Q);
		}
	}
	fich2.close();

	return 0;
}

void mostrarDatosProblema(int numNiveles, tipoDatosSubmalla datosNivel, int numVolxNivel0, int numVolyNivel0,
				double tiempo_tot, double tiempoGuardarNetCDF, double CFL, double mf0, double vmax,
				double epsilon_h, double hpos, double cvis, double L, double H, double Q, double T)
{
	cout << "/**********************************************************************" << endl;
	cout << " Tsunami-HySEA numerical model                                         " << endl;
	cout << " developed by the EDANYA Research Group, University of Malaga (Spain). " << endl;
	cout << " https://www.uma.es/edanya                                             " << endl;
	cout << " https://edanya.uma.es/hysea/                                          " << endl;
	cout << " Contact and support: hysea@uma.es                                     " << endl;
	cout << " The core version of Tsunami-HySEA is distributed under                " << endl;
	cout << " Creative Commons Attribution-NonCommercial-NoDerivatives 4.0          " << endl;
	cout << " International Public License (CC BY-NC-ND 4.0)                        " << endl;
	cout << "**********************************************************************/" << endl;
	cout << endl;

	cout << "Problem data" << endl;
	cout << "CFL: " << CFL << endl;
	cout << "Water-bottom friction: " << sqrt(mf0*pow(H,4.0/3.0)/(9.81*L)) << endl;
	cout << "Maximum allowed velocity of water: " << vmax*(Q/H) << endl;
	cout << "Epsilon h: " << epsilon_h*H << " m" << endl;
	cout << "Threshold for the 2s+WAF scheme: " << hpos*H << " m" << endl;
	cout << "Stability coefficient: " << 1.0-cvis << endl;
	cout << "Simulation time: " << tiempo_tot*T << " sec" << endl;
	cout << "Saving time of NetCDF files: " << tiempoGuardarNetCDF*T << " sec" << endl;
	cout << "Number of levels: " << numNiveles << endl;
	cout << "Level 0" << endl;
	cout << "  Volumes: " << numVolxNivel0 << " x " << numVolyNivel0 << " = " << numVolxNivel0*numVolyNivel0 << endl;
	cout << "  Longitude: [" << datosNivel.longitud[0] << ", " << datosNivel.longitud[numVolxNivel0-1] << "]" << endl;
	cout << "  Latitude: [" << datosNivel.latitud[0] << ", " << datosNivel.latitud[numVolyNivel0-1] << "]" << endl;
	cout << endl;
}

