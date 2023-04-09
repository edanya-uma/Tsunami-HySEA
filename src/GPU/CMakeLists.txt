INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

CUDA_ADD_LIBRARY(GPU
			ComprobarSoporteCUDA.cu
			ShallowWater.cu
			Arista_kernel.cu
			Reduccion_kernel.cu
			Volumen_kernel.cu
			netcdf.cu
			netcdfSeries.cu
			Matriz.cu
			Deformacion.cu
			DC3D.cu
			prtoxy.cu)


