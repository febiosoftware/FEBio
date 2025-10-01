#include "FEImgLib.h"
#include <FECore/FECoreKernel.h>
#include "FEImageSource.h"
#include "FEImageDataMap.h"
#include "FEImageValuator.h"
#include "ImageFilter.h"

void FEImgLib::InitModule()
{
	// image sources
	REGISTER_FECORE_CLASS(FERawImage, "raw");
	REGISTER_FECORE_CLASS(FENRRDImage, "nrrd");

	// data maps
	REGISTER_FECORE_CLASS(FEImageDataMap, "image map");

	// valuator
	REGISTER_FECORE_CLASS(FEImageValuator, "image map");

	// filter classes
	REGISTER_FECORE_CLASS(IterativeBlur, "iterative blur");
	REGISTER_FECORE_CLASS(BoxBlur, "box blur");
#ifdef HAVE_FFTW
	REGISTER_FECORE_CLASS(FFTWBlur, "FFTW blur");
#endif // HAVE_FFTW
}
