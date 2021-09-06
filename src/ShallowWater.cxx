/**********************************************************************
 Tsunami-HySEA numerical model
 developed by the EDANYA Research Group, University of Malaga (Spain).
 https://www.uma.es/edanya
 https://edanya.uma.es/hysea/
 Contact and support: hysea@uma.es
 The core version of Tsunami-HySEA is distributed under
 Creative Commons Attribution-NonCommercial-NoDerivatives 4.0
 International Public License (CC BY-NC-ND 4.0)
**********************************************************************/

#include <cstring>
#include "Constantes.hxx"
#include "Problema.cxx"

/*****************/
/* Funciones GPU */
/*****************/

extern "C" int comprobarSoporteCUDA();
extern "C" int shallowWater(int numNiveles, double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2,
			tipoDatosSubmalla datosNivel, int leer_fichero_puntos, int numVolxNivel0, int numVolyNivel0,
			int numVolumenesNivel, double Hmin, char *nombre_bati, string prefijo, double borde_sup, double borde_inf,
			double borde_izq, double borde_der, double tiempo_tot, double tiempoGuardarNetCDF, double CFL, double mf0,
			double vmax, double epsilon_h, double hpos, double cvis, double L, double H, double Q, double T, double *tiempo);

/*********************/
/* Fin funciones GPU */
/*********************/

int main(int argc, char *argv[])
{
	double2 *datosVolumenesNivel_1;
	double2 *datosVolumenesNivel_2;
	tipoDatosSubmalla datosNivel;
	int leer_fichero_puntos;
	int soporteCUDA, err;
	int numVolxNivel0, numVolyNivel0;
	int numVolumenesNivel;
	double borde_sup, borde_inf, borde_izq, borde_der;
	double tiempo_tot;
	double tiempoGuardarNetCDF;
	double Hmin, epsilon_h;
	double mf0, vmax;
	double CFL, hpos, cvis;
	double L, H, Q, T;
	string fich_ent;
	string nombre_bati;
	string prefijo;
	int numNiveles;
	double tiempo_gpu;

	if (argc < 2) {
		cerr << "/**********************************************************************" << endl;
		cerr << " Tsunami-HySEA numerical model                                         " << endl;
		cerr << " developed by the EDANYA Research Group, University of Malaga (Spain). " << endl;
		cerr << " https://www.uma.es/edanya                                             " << endl;
		cerr << " https://edanya.uma.es/hysea/                                          " << endl;
		cerr << " Contact and support: hysea@uma.es                                     " << endl;
		cerr << " The core version of Tsunami-HySEA is distributed under                " << endl;
		cerr << " Creative Commons Attribution-NonCommercial-NoDerivatives 4.0          " << endl;
		cerr << " International Public License (CC BY-NC-ND 4.0)                        " << endl;
		cerr << "**********************************************************************/" << endl;
		cerr << endl;
		cerr << "Use: " << endl;
		cerr << argv[0] << " dataFile" << endl << endl; 
		cerr << "dataFile format:" << endl;
		cerr << "  Bathymetry name" << endl;
		cerr << "  Number of columns" << endl;
		cerr << "  Number of rows" << endl;
		cerr << "  Bathymetry file" << endl;
		cerr << "  Initial state file" << endl;
		cerr << "  NetCDF file prefix" << endl;
		cerr << "  Number of levels (should be 1)" << endl;
		cerr << "  Upper border condition (1: open, -1: wall)" << endl;
		cerr << "  Lower border condition" << endl;
		cerr << "  Left border condition" << endl;
		cerr << "  Right border condition" << endl;
		cerr << "  Simulation time (sec)" << endl;
		cerr << "  Saving time of NetCDF files (sec) (-1: do not save)" << endl;
		cerr << "  Read points from file (should be 0)" << endl;
		cerr << "  CFL" << endl;
		cerr << "  Epsilon h (m)" << endl;
		cerr << "  Threshold for the 2s+WAF scheme (m)" << endl;
		cerr << "  Stability coefficient" << endl;
		cerr << "  Water-bottom friction (Manning coefficient)" << endl;
		cerr << "  Maximum allowed velocity of water" << endl;
		cerr << "  L (typical length)" << endl;
		cerr << "  H (typical depth)" << endl;
		return EXIT_FAILURE;
	}

	fich_ent = argv[1];
	if (! existeFichero(fich_ent)) {
		cerr << "Error: File '" << fich_ent << "' not found" << endl;
		return EXIT_FAILURE;
	}

	soporteCUDA = comprobarSoporteCUDA();
	if (soporteCUDA == 1) {
		cerr << "Error: There is no graphics card" << endl;
		return EXIT_FAILURE;
	}
	else if (soporteCUDA == 2) {
		cerr << "Error: There is no graphics card supporting CUDA" << endl;
		return EXIT_FAILURE;
	}

	err = cargarDatosProblema(fich_ent, nombre_bati, prefijo, &numNiveles, &datosVolumenesNivel_1, &datosVolumenesNivel_2,
			&datosNivel, &numVolxNivel0, &numVolyNivel0, &numVolumenesNivel, &leer_fichero_puntos, &Hmin, &borde_sup,
			&borde_inf, &borde_izq, &borde_der, &tiempo_tot, &tiempoGuardarNetCDF, &CFL, &mf0, &vmax, &epsilon_h,
			&hpos, &cvis, &L, &H, &Q, &T);
	if (err > 0)
		return EXIT_FAILURE;

	mostrarDatosProblema(numNiveles, datosNivel, numVolxNivel0, numVolyNivel0, tiempo_tot, tiempoGuardarNetCDF,
		CFL, mf0, vmax, epsilon_h, hpos, cvis, L, H, Q, T);

	cout << scientific;
	err = shallowWater(numNiveles, datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel, leer_fichero_puntos,
			numVolxNivel0, numVolyNivel0, numVolumenesNivel, Hmin, (char *) nombre_bati.c_str(), prefijo,
			borde_sup, borde_inf, borde_izq, borde_der, tiempo_tot, tiempoGuardarNetCDF, CFL, mf0, vmax,
			epsilon_h, hpos, cvis, L, H, Q, T, &tiempo_gpu);
	if (err > 0) {
		if (err == 1)
			cerr << "Error: Not enough GPU memory" << endl;
		else if (err == 2)
			cerr << "Error: Not enough CPU memory" << endl;
		liberarMemoria(numNiveles, datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel);
		return EXIT_FAILURE;
	}
	cout << endl << "Runtime: " << tiempo_gpu << " sec" << endl;
	liberarMemoria(numNiveles, datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel);

	return EXIT_SUCCESS;
}

