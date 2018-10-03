#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

BEGIN_FECORE_CLASS(FEElasticFiberMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_thd, "theta");
	ADD_PARAMETER(m_phd, "phi");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEElasticFiberMaterial::FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem) 
{
	m_thd = 0.0;
	m_phd = 90.0;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticFiberMaterial::CreateMaterialPointData()
{
	FEFiberMaterialPoint* fp = new FEFiberMaterialPoint(FEElasticMaterial::CreateMaterialPointData());

	// Some fiber materials defined the theta,phi parameters for setting the fiber vector
	// Although this is deprecated, we still support it here for backward compatibility
	if ((m_thd != 0.0) || (m_phd != 90.0))
	{
		// convert angles from degrees to radians
		double pi = 4 * atan(1.0);
		double the = m_thd*pi / 180.;
		double phi = m_phd*pi / 180.;

		// fiber direction in local coordinate system (reference configuration)
		vec3d n0;
		n0.x = cos(the)*sin(phi);
		n0.y = sin(the)*sin(phi);
		n0.z = cos(phi);
		n0.unit();
		fp->m_n0 = n0;
	}

	return fp;
}

//-----------------------------------------------------------------------------
vec3d FEElasticFiberMaterial::GetFiberVector(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	return pt.m_Q*fp.m_n0;
}
