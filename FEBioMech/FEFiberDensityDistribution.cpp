//  Created by Gerard Ateshian on 11/16/13.
//

#include "FEFiberDensityDistribution.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// define the ellipsoidal fiber density distributionmaterial parameters
BEGIN_PARAMETER_LIST(FEEllipsodialFiberDensityDistribution, FEMaterial)
	ADD_PARAMETERV(m_spa , FE_PARAM_DOUBLEV, 3, "spa" );
END_PARAMETER_LIST();

void FEEllipsodialFiberDensityDistribution::Init()
{
	FEMaterial::Init();
	if (m_spa[0] < 0) throw MaterialError("spa1 must be positive.");
	if (m_spa[1] < 0) throw MaterialError("spa2 must be positive.");
	if (m_spa[2] < 0) throw MaterialError("spa3 must be positive.");
}

double FEEllipsodialFiberDensityDistribution::FiberDensity(const vec3d n0)
{
    double R = 1.0/sqrt(SQR(n0.x/m_spa[0])+SQR(n0.y/m_spa[1])+SQR(n0.z/m_spa[2]));
    return R/m_IFD;
}

//-----------------------------------------------------------------------------
// define the 3d von Mises fiber density distribution material parameters
BEGIN_PARAMETER_LIST(FEVonMises3DFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER2(m_b, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "b" );
END_PARAMETER_LIST();

double FEVonMises3DFiberDensityDistribution::FiberDensity(const vec3d n0)
{
    // The local x-direction is the principal fiber bundle direction
    // The x-component of n0 is cos(phi)
    double R = exp(m_b*(2*SQR(n0.x)-1));
    return R/m_IFD;
}

//-----------------------------------------------------------------------------
// define the ellipsoidal fiber density distributionmaterial parameters
BEGIN_PARAMETER_LIST(FEEllipticalFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER2(m_spa[0], FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "spa1" );
	ADD_PARAMETER2(m_spa[1], FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "spa2" );
END_PARAMETER_LIST();

double FEEllipticalFiberDensityDistribution::FiberDensity(const vec3d n0)
{
    // 2d fibers lie in the local x-y plane
    // n0.x = cos(theta) and n0.y = sin(theta)
    double R = 1.0/sqrt(SQR(n0.x/m_spa[0])+SQR(n0.y/m_spa[1]));
    return R/m_IFD;
}

//-----------------------------------------------------------------------------
// define the 2d von Mises fiber density distribution material parameters
BEGIN_PARAMETER_LIST(FEVonMises2DFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER2(m_b, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "b" );
END_PARAMETER_LIST();

double FEVonMises2DFiberDensityDistribution::FiberDensity(const vec3d n0)
{
    // The fiber bundle is in the x-y plane and
    // the local x-direction is the principal fiber bundle direction
    // The x-component of n0 is cos(theta)
    double R = exp(m_b*(2*SQR(n0.x)-1));
    return R/m_IFD;
}
