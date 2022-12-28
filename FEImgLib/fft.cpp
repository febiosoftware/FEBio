#ifdef HAVE_MKL
#include <mkl.h>
#include <complex>

//---------------------------------------------------------------------------------------
// calculate the DFT of a source image x. The source image is assumed to be stored
// in row-major order, indexed by (x,y), where x is the column index, and y the row index. 
// The complex Fourier coefficients are returned in c. The buffer must be pre-allocated and must have
// the size (nx,ny), the same as the source image.
bool mkl_dft2(int nx, int ny, float* x, MKL_Complex8* c)
{
	DFTI_DESCRIPTOR_HANDLE my_desc_handle;
	MKL_LONG status;

	MKL_LONG sizes[2] = { nx, ny };
	MKL_LONG is[3] = { 0, 1, nx };
	MKL_LONG os[3] = { 0, 1, nx }; // { 0, 1, N2/2+1}
	status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 2, sizes);
	status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, is);
	status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, os);
	status = DftiCommitDescriptor(my_desc_handle);
	status = DftiComputeForward(my_desc_handle, x, c);
	status = DftiFreeDescriptor(&my_desc_handle);

	return true;
}

//---------------------------------------------------------------------------------------
// calculate the inverse DFT. The fourier coefficients are passed and assumed to be stored
// in row-major order, indexed by (x,y), where x is the column index, and y the row index. 
// The reconstructed image is returned in x, which must have the same size and storage. 
bool mkl_idft2(int nx, int ny, MKL_Complex8* c, float* x)
{
	DFTI_DESCRIPTOR_HANDLE my_desc_handle;
	MKL_LONG status;

	MKL_LONG sizes[2] = { nx, ny };
	MKL_LONG is[3] = { 0, 1, nx };  // { 0, 1, Nx/2+1}
	MKL_LONG os[3] = { 0, 1, nx };
	status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 2, sizes);
	status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, is);
	status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, os);
	status = DftiSetValue(my_desc_handle, DFTI_BACKWARD_SCALE, 1.f / ((float)nx * (float)ny));
	status = DftiCommitDescriptor(my_desc_handle);
	status = DftiComputeBackward(my_desc_handle, c, x);
	status = DftiFreeDescriptor(&my_desc_handle);

	return true;
}

//---------------------------------------------------------------------------------------
// calculate the DFT of a source image x. The source image is assumed to be stored
// in row-major order, indexed by (x,y,z), where x is the column index, y the row index, z is the plane index. 
// The complex Fourier coefficients are returned in c. The buffer must be pre-allocated and must have
// the size (nx,ny,nz), the same as the source image.
bool mkl_dft3(int nx, int ny, int nz, float* x, MKL_Complex8* c)
{
	DFTI_DESCRIPTOR_HANDLE my_desc_handle;
	MKL_LONG status;

	MKL_LONG sizes[3] = { nx, ny, nz };
	MKL_LONG is[4] = { 0, 1, nx, nx*ny };
	MKL_LONG os[4] = { 0, 1, nx, nx*ny }; // { 0, 1, nx, ny, nz/2+1}
	status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, sizes);
	status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, is);
	status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, os);
	status = DftiCommitDescriptor(my_desc_handle);
	status = DftiComputeForward(my_desc_handle, x, c);
	status = DftiFreeDescriptor(&my_desc_handle);

	return true;
}

//---------------------------------------------------------------------------------------
// calculate the inverse DFT. The fourier coefficients are passed and assumed to be stored
// in row-major order, indexed by (x,y,z), where x is the column index, y the row index, and z the plane index. 
// The reconstructed image is returned in y, which must have the same size and storage. 
bool mkl_idft3(int nx, int ny, int nz, MKL_Complex8* c, float* y)
{
	DFTI_DESCRIPTOR_HANDLE my_desc_handle;
	MKL_LONG status;

	MKL_LONG sizes[3] = { nx, ny, nz };
	MKL_LONG is[4] = { 0, 1, nx, nx * ny };  // { 0, 1, nx, ny, nz/2+1}
	MKL_LONG os[4] = { 0, 1, nx, nx * ny };
	status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, sizes);
	status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, is);
	status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, os);
	status = DftiSetValue(my_desc_handle, DFTI_BACKWARD_SCALE, 1.f / ((float)nx * (float)ny * (float) nz));
	status = DftiCommitDescriptor(my_desc_handle);
	status = DftiComputeBackward(my_desc_handle, c, y);
	status = DftiFreeDescriptor(&my_desc_handle);

	return true;
}

#endif
