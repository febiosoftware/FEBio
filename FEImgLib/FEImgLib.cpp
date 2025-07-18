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
	REGISTER_FECORE_CLASS(IterativeBlur1D, "iterative blur 1D");
	REGISTER_FECORE_CLASS(IterativeBlur3D, "iterative blur 3D");
	REGISTER_FECORE_CLASS(BoxBlur1D, "box blur 1D");
	REGISTER_FECORE_CLASS(BoxBlur3D, "box blur 3D");
}
