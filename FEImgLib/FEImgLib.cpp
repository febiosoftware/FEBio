#include "FEImgLib.h"
#include <FECore/FECoreKernel.h>
#include "FEImageSource.h"
#include "FEImageDataMap.h"
#include "FEImageValuator.h"

void FEImgLib::InitModule()
{
	// image sources
	REGISTER_FECORE_CLASS(FERawImage, "raw");

	// data maps
	REGISTER_FECORE_CLASS(FEImageDataMap, "image map");

	// valuator
	REGISTER_FECORE_CLASS(FEImageValuator, "image map");
}
