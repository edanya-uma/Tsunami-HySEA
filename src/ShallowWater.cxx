/**********************************************************************
 Tsunami-HySEA numerical model open source v1.3.0
 developed by the EDANYA Research Group, University of Malaga (Spain).
 https://www.uma.es/edanya
 https://edanya.uma.es/hysea/
 Contact and support: hysea@uma.es
 The core version of Tsunami-HySEA is distributed under
 GNU General Public License, Version 2
**********************************************************************/

#include <string>
#include "Constantes.hxx"
#include "Problema.cxx"
#include "GPU/SoporteCUDA.cu"

using namespace std;
using namespace SoporteCUDA;

extern "C" int shallowWater(int numNiveles, int okada_flag, int numFaults, int crop_flag, double crop_value, double *LON_C,
			double *LAT_C, double *DEPTH_C, double *FAULT_L, double *FAULT_W, double *STRIKE, double *DIP, double *RAKE,
			double *SLIP, double *defTime, double2 *datosVolumenesNivel_1, double2 *datosVolumenesNivel_2, tipoDatosSubmalla datosNivel,
			int leer_fichero_puntos, int numPuntosGuardar, int *posicionesVolumenesGuardado, double *lonPuntos, double *latPuntos,
			int numVolxNivel0, int numVolyNivel0, int64_t numVolumenesNivel, double Hmin, char *nombre_bati, string prefijo,
			double borde_sup, double borde_inf, double borde_izq, double borde_der, int tam_spongeSup, int tam_spongeInf,
			int tam_spongeIzq, int tam_spongeDer, double tiempo_tot, double tiempoGuardarNetCDF, double tiempoGuardarSeries,
			double CFL, double mf0, double vmax, double epsilon_h, double hpos, double cvis, double dif_at, double L, double H,
			double Q, double T, char *version, double *tiempo);

int main(int argc, char *argv[])
{
	string version = "1.3.0";
	double2 *datosVolumenesNivel_1;
	double2 *datosVolumenesNivel_2;
	tipoDatosSubmalla datosNivel;
	int *posicionesVolumenesGuardado = nullptr;
	double *lonPuntos = nullptr;
	double *latPuntos = nullptr;
	int leer_fichero_puntos;
	int numPuntosGuardar;
	int soporte_CUDA, err;
	int numVolxNivel0, numVolyNivel0;
	int64_t numVolumenesNivel;
	double borde_sup, borde_inf, borde_izq, borde_der;
	int tam_spongeSup, tam_spongeInf, tam_spongeIzq, tam_spongeDer;
	double tiempo_tot;
	double tiempoGuardarNetCDF, tiempoGuardarSeries;
	double Hmin, epsilon_h, dif_at;
	double mf0, vmax;
	double CFL, hpos, cvis;
	double L, H, Q, T;
	string fich_ent;
	string nombre_bati;
	string prefijo;
	int numNiveles;
	double tiempo_gpu;
	int okada_flag, numFaults;
	int crop_flag;
	double crop_value;
	double LON_C[MAX_FAULTS], LAT_C[MAX_FAULTS], DEPTH_C[MAX_FAULTS], FAULT_L[MAX_FAULTS], FAULT_W[MAX_FAULTS];
	double STRIKE[MAX_FAULTS], DIP[MAX_FAULTS], RAKE[MAX_FAULTS], SLIP[MAX_FAULTS], defTime[MAX_FAULTS];

	if (argc < 2) {
		cerr << "/**********************************************************************" << endl;
		cerr << " Tsunami-HySEA numerical model open source v" << version << endl;
		cerr << " developed by the EDANYA Research Group, University of Malaga (Spain). " << endl;
		cerr << " https://www.uma.es/edanya                                             " << endl;
		cerr << " https://edanya.uma.es/hysea/                                          " << endl;
		cerr << " Contact and support: hysea@uma.es                                     " << endl;
		cerr << " The core version of Tsunami-HySEA is distributed under                " << endl;
		cerr << " GNU General Public License, Version 2                                 " << endl;
		cerr << "**********************************************************************/" << endl;
		cerr << endl;
		cerr << "Use: " << endl;
		cerr << argv[0] << " dataFile" << endl << endl;
		cerr << "dataFile format:" << endl;
		cerr << "  Problem name" << endl;
		cerr << "  Bathymetry file" << endl;
		cerr << "  Initialization of states (" << SEA_SURFACE_FROM_FILE << ": Sea surface displacement from file," << endl;
		cerr << "                            " << OKADA_STANDARD << ": Standard Okada)" << endl;
		cerr << "  If " << SEA_SURFACE_FROM_FILE << ":" << endl;
		cerr << "    Initial state file" << endl;
		cerr << "  Else if " << OKADA_STANDARD << ":" << endl;
		cerr << "    Number of faults (>= 1)" << endl;
		cerr << "    For every fault, a line containing:" << endl;
		cerr << "      Time(sec) Lon_epicenter Lat_epicenter Depth_hypocenter(km) Fault_length(km) Fault_width(km) Strike Dip Rake Slip(m)" << endl;
		cerr << "    Crop deformations (" << NO_CROP << ": No, " << CROP_RELATIVE << ": Relative cropping, " << CROP_ABSOLUTE << ": Absolute cropping)" << endl;
		cerr << "    If " << CROP_RELATIVE << ":" << endl;
		cerr << "      Percentage [0,1]" << endl;
		cerr << "    Else if " << CROP_ABSOLUTE << ":" << endl;
		cerr << "      Threshold value (m)" << endl;
		cerr << "  NetCDF file prefix" << endl;
		cerr << "  Number of levels (should be 1)" << endl;
		cerr << "  Upper border condition (1: open, -1: wall)" << endl;
		cerr << "  Lower border condition (1: open, -1: wall)" << endl;
		cerr << "  Left border condition (1: open, -1: wall, 2: periodic)" << endl;
		cerr << "  Right border condition (1: open, -1: wall, 2: periodic)" << endl;
		cerr << "  Simulation time (sec)" << endl;
		cerr << "  Saving time of NetCDF files (sec) (-1: do not save)" << endl;
		cerr << "  Read points from file (0: no, 1: yes)" << endl;
		cerr << "  If 1:" << endl;
		cerr << "    File with points" << endl;
		cerr << "    Saving time of time series (sec)" << endl;
		cerr << "  CFL" << endl;
		cerr << "  Epsilon h (m)" << endl;
		cerr << "  Threshold for the 2s+WAF scheme (m)" << endl;
		cerr << "  Stability coefficient" << endl;
		cerr << "  Water-bottom friction (Manning coefficient)" << endl;
		cerr << "  Maximum allowed velocity of water" << endl;
		cerr << "  L (typical length)" << endl;
		cerr << "  H (typical depth)" << endl;
		cerr << "  Threshold for arrival times (m) (epsilon h if not specified)" << endl;
		return EXIT_FAILURE;
	}

	fich_ent = argv[1];
	if (! existeFichero(fich_ent)) {
		cerr << "Error: File '" << fich_ent << "' not found" << endl;
		return EXIT_FAILURE;
	}

	soporte_CUDA = comprobarSoporteCUDA();
	if (soporte_CUDA == 1) {
		cerr << "Error: There is no graphics card" << endl;
		return EXIT_FAILURE;
	}
	else if (soporte_CUDA == 2) {
		cerr << "Error: There is no graphics card supporting CUDA" << endl;
		return EXIT_FAILURE;
	}

	err = cargarDatosProblema(fich_ent, nombre_bati, prefijo, &numNiveles, &okada_flag, &numFaults, &crop_flag,
			&crop_value, LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W, STRIKE, DIP, RAKE, SLIP, defTime, &datosVolumenesNivel_1,
			&datosVolumenesNivel_2, &datosNivel, &numVolxNivel0, &numVolyNivel0, &numVolumenesNivel, &leer_fichero_puntos,
			&numPuntosGuardar, &posicionesVolumenesGuardado, &lonPuntos, &latPuntos, &Hmin, &borde_sup, &borde_inf,
			&borde_izq, &borde_der, &tam_spongeSup, &tam_spongeInf, &tam_spongeIzq, &tam_spongeDer, &tiempo_tot,
			&tiempoGuardarNetCDF, &tiempoGuardarSeries, &CFL, &mf0, &vmax, &epsilon_h, &hpos, &cvis, &dif_at,
			&L, &H, &Q, &T);
	if (err > 0)
		return EXIT_FAILURE;

	mostrarDatosProblema(version, numNiveles, okada_flag, numFaults, crop_flag, crop_value, datosNivel, numVolxNivel0,
		numVolyNivel0, tiempo_tot, tiempoGuardarNetCDF, leer_fichero_puntos, tiempoGuardarSeries, CFL, mf0, vmax,
		epsilon_h, hpos, cvis, dif_at, L, H, Q, T);

	cout << scientific;
	err = shallowWater(numNiveles, okada_flag, numFaults, crop_flag, crop_value, LON_C, LAT_C, DEPTH_C, FAULT_L, FAULT_W,
			STRIKE, DIP, RAKE, SLIP, defTime, datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel, leer_fichero_puntos,
			numPuntosGuardar, posicionesVolumenesGuardado, lonPuntos, latPuntos, numVolxNivel0, numVolyNivel0, numVolumenesNivel,
			Hmin, const_cast<char *>(nombre_bati.c_str()), prefijo, borde_sup, borde_inf, borde_izq, borde_der, tam_spongeSup,
			tam_spongeInf, tam_spongeIzq, tam_spongeDer, tiempo_tot, tiempoGuardarNetCDF, tiempoGuardarSeries, CFL, mf0,
			vmax, epsilon_h, hpos, cvis, dif_at, L, H, Q, T, const_cast<char *>(version.c_str()), &tiempo_gpu);
	if (err > 0) {
		if (err == 1)
			cerr << "Error: Not enough GPU memory" << endl;
		else if (err == 2)
			cerr << "Error: Not enough CPU memory" << endl;
		liberarMemoria(datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel,
			posicionesVolumenesGuardado, lonPuntos, latPuntos);
		return EXIT_FAILURE;
	}
	cout << endl << "Runtime: " << tiempo_gpu << " sec" << endl;
	liberarMemoria(datosVolumenesNivel_1, datosVolumenesNivel_2, datosNivel,
		posicionesVolumenesGuardado, lonPuntos, latPuntos);

	return EXIT_SUCCESS;
}

