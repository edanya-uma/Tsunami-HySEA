#ifndef REDUCCION_KERNEL_H
#define REDUCCION_KERNEL_H

#include "sharedmem.cuh"

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)
*/
template <class T, unsigned int blockSize>
__global__ void reduce6_min(T *g_idata, T *g_odata, unsigned int n)
{
    SharedMemory<T> smem;
    T *sdata = smem.getPointer();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
	T thMin, val;

	if (i + blockSize >= n)
		thMin = g_idata[i];
	else
		thMin = min(g_idata[i], g_idata[i + blockSize]);
	i += gridSize;
	while (i < n) {
		if (i + blockSize >= n)
			val = g_idata[i];
		else
			val = min(g_idata[i], g_idata[i + blockSize]);
		thMin = min(thMin, val);
		i += gridSize;
	}
	sdata[tid] = thMin;

	__syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) {
		if (tid < 256) { sdata[tid] = thMin = min(sdata[tid], sdata[tid + 256]); }
		__syncthreads();
	}
    if (blockSize >= 256) {
		if (tid < 128) { sdata[tid] = thMin = min(sdata[tid], sdata[tid + 128]); }
		__syncthreads();
	}
    if (blockSize >= 128) {
		if (tid <  64) { sdata[tid] = thMin = min(sdata[tid], sdata[tid +  64]); }
		__syncthreads();
	}

/*#if (__CUDA_ARCH__ >= 300)
	if (tid < 32) {
		int warpSize = 32;

		// Fetch final intermediate sum from 2nd warp
		if (blockSize >=  64) thMin = min(thMin, sdata[tid + 32]);
		// Reduce final warp using shuffle
		for (int offset = warpSize/2; offset > 0; offset /= 2) {
			thMin = min(thMin, __shfl_down(thMin, offset));
		}
	}
#else*/
	// fully unroll reduction within a single warp
	if (tid < 32) {
		volatile T* sdatav = sdata;
		if (blockSize >=  64) sdatav[tid] = thMin = min(thMin, sdatav[tid + 32]);
		if (blockSize >=  32) sdatav[tid] = thMin = min(thMin, sdatav[tid + 16]);
		if (blockSize >=  16) sdatav[tid] = thMin = min(thMin, sdatav[tid +  8]);
		if (blockSize >=   8) sdatav[tid] = thMin = min(thMin, sdatav[tid +  4]);
		if (blockSize >=   4) sdatav[tid] = thMin = min(thMin, sdatav[tid +  2]);
		if (blockSize >=   2) sdatav[tid] = thMin = min(thMin, sdatav[tid +  1]);
	}
//#endif

	// write result for this block to global mem 
	if (tid == 0) g_odata[blockIdx.x] = thMin;
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void reduce_min(int blocks, int threads, T *d_idata, T *d_odata, int size)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);
	int smemSize = threads * sizeof(T);

	switch (threads) {
	case 512:
		reduce6_min<T, 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 256:
		reduce6_min<T, 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 128:
		reduce6_min<T, 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 64:
		reduce6_min<T,  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 32:
		reduce6_min<T,  32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 16:
		reduce6_min<T,  16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  8:
		reduce6_min<T,   8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  4:
		reduce6_min<T,   4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  2:
		reduce6_min<T,   2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  1:
		reduce6_min<T,   1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	}
}

template
void reduce_min<float>(int blocks, int threads, float *d_idata, float *d_odata, int size);

template
void reduce_min<double>(int blocks, int threads, double *d_idata, double *d_odata, int size);

// Función que pone en blocks y threads el número de bloques y de hebras por bloque,
// respectivamente, para el kernel de reducción. Tenemos en cuenta el máximo número de
// bloques porque cada hebra puede procesar un número variable de elementos.
void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
/*	if (n == 1) 
		threads = 1;
	else*/
		threads = (n < maxThreads*2) ? n / 2 : maxThreads;
	blocks = n / (threads * 2);
	blocks = min(maxBlocks, blocks);
}

template <class T>
T obtenerMinimoReduccion(T *d_data, int size)
{
	int maxBlocks = 128, maxThreads = 256;
	int numBlocks, numThreads;
	T h_data[16];
	T minimo;

	if (size > 1) {
		// Ejecutamos el kernel
		getNumBlocksAndThreads(size, maxBlocks, maxThreads, numBlocks, numThreads);
		reduce_min<T>(numBlocks, numThreads, d_data, d_data, size);

		// Obtenemos el mínimo de los resultados parciales de los bloques en GPU
		// hasta que el número de elementos sea 1
		int s = numBlocks;
		while (s > 1) {
			getNumBlocksAndThreads(s, maxBlocks, maxThreads, numBlocks, numThreads);
			reduce_min<T>(numBlocks, numThreads, d_data, d_data, s);
			s = s / (numThreads*2);
		}

		// Copiamos el resultado de la reducción de GPU a CPU
		cudaMemcpy(h_data, d_data, numBlocks*sizeof(T), cudaMemcpyDeviceToHost);
		minimo = (T) h_data[0];
	}
	else {
		// Copiamos el resultado de la reducción de GPU a CPU
		cudaMemcpy(h_data, d_data, size*sizeof(T), cudaMemcpyDeviceToHost);
		minimo = (T) h_data[0];
	}

	return minimo;
}

#endif
