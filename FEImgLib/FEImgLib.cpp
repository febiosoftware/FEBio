#include "FEImgLib.h"
#include <FECore/FECoreKernel.h>
#include "FEImageSource.h"

void FEImgLib::InitModule()
{
	// image sources
	REGISTER_FECORE_CLASS(FERawImage, "raw");
}
