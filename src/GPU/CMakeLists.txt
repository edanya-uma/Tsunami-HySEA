INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

CUDA_ADD_LIBRARY(GPU
			SoporteCUDA.cu
			ShallowWater.cu
			Arista_kernel.cu
			Reduccion_kernel.cu
			Volumen_kernel.cu
			Timestep.cu
			netcdfTH.cu
			netcdfSeries.cu
			Metrics.cu
			netcdf_kernel.cu
			Deformacion.cu
			DC3D.cu
			prtoxy.cu)

TARGET_LINK_LIBRARIES(GPU ${NETCDF_LIBRARIES})

