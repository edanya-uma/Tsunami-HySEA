#include <stdio.h>

extern "C" int comprobarSoporteCUDA()
{
    int valor, dev, deviceCount;

    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
        valor = 1;
	else {
		for (dev=0; dev < deviceCount; ++dev) {
	        cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);
			if (deviceProp.major >= 1)
	            break;
	    }
		if (dev == deviceCount)
			valor = 2;
		else
			valor = 0;
	}
	return valor;
}
