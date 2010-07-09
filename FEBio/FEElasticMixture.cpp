#include "stdafx.h"
#include "tens4d.h"
#include "FEElasticMixture.h"

// register the material with the framework
REGISTER_MATERIAL(FEElasticMixture, "solid mixture");

// define the material parameters
BEGIN_PARAMETER_LIST(FEElasticMixture, FEElasticMaterial)
ADD_PARAMETER(m_nMat, FE_PARAM_INT, "nsolids");
ADD_PARAMETERV(m_iMat, FE_PARAM_INTV, NMAT, "solid_ids");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Mixture of elastic solids
//////////////////////////////////////////////////////////////////////

void FEElasticMixture::Init()
{
	FEElasticMaterial::Init();

	if (m_nMat < 1) throw MaterialError("nsolid must be greater than zero");
	for (int i=0; i<m_nMat-1; ++i)
		for (int j=i+1; j<m_nMat; ++j)
			if (m_iMat[i] == m_iMat[j])
				throw MaterialError("Duplicate material id specified in solid mixture.");
	
}

mat3ds FEElasticMixture::Stress(FEMaterialPoint& mp)
{
	mat3ds s;
	
	// calculate stress
	s = m_pMat[0]->Stress(mp);
	for (int i=1; i < m_nMat; ++i)
		s += m_pMat[i]->Stress(mp);

	return s;
}

tens4ds FEElasticMixture::Tangent(FEMaterialPoint& mp)
{
	tens4ds c;

	// calculate elasticity tensor
	c = m_pMat[0]->Tangent(mp);
	for (int i=1; i < m_nMat; ++i)
		c += m_pMat[i]->Tangent(mp);

	return c;
}

double FEElasticMixture::BulkModulus()
{
	double k;
	
	// calculate bulk modulus
	k = m_pMat[0]->BulkModulus();
	for (int i=1; i < m_nMat; ++i)
		k += m_pMat[i]->BulkModulus();
	
	return k;
}

