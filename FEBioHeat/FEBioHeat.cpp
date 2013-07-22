#include "FEBioHeat.h"
#include "FEIsotropicFourier.h"

namespace FEBioHeat {

void InitModule()
{
REGISTER_MATERIAL(FEIsotropicFourier, "isotropic Fourier");
}

}
