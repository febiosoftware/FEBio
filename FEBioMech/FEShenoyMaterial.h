#pragma once
#include <FEBioMech/FEElasticMaterial.h>

// This plugin implements the constitutive model by 
// Wang et al. (Biophysical Journal, 107, 2014, pp:2592 - 2603).
// This novel material proposes a mechanism for long-range force transmission in 
// fibrous matrices enabled by tension-driven alignment of fibers.
class FEShenoyMaterial : public FEElasticMaterial
{
public:
	FEShenoyMaterial(FEModel* fem);

	mat3ds Stress(FEMaterialPoint& mp) override;

	tens4ds Tangent(FEMaterialPoint& mp) override;

private:
	double fiberStress(double lam);
	double fiberTangent(double lam);

private:
	double	m_mu;
	double	m_k;
	double	m_Ef;
	double	m_lamc;
	double	m_lamt;
	double	m_n;
	double	m_m;

	DECLARE_FECORE_CLASS();
};
