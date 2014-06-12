#include "stdafx.h"
#include "FEPreStrainCoupledTransIsoMR.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEPreStrainCoupledTransIsoMR, FEElasticMaterial)
	ADD_PARAMETER(m_mat.m_c1  , FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_mat.m_c2  , FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_mat.m_c3  , FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_mat.m_c4  , FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_mat.m_c5  , FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_mat.m_flam, FE_PARAM_DOUBLE, "lambda");
	ADD_PARAMETER(m_mat.m_K   , FE_PARAM_DOUBLE, "k");
	ADD_PARAMETER(m_ltrg	  , FE_PARAM_DOUBLE, "pre_stretch");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPreStrainCoupledTransIsoMR::FEPreStrainCoupledTransIsoMR(FEModel* pfem)  : FEElasticMaterial(pfem), m_mat(pfem)
{
}

//-----------------------------------------------------------------------------
//! create material point data for this material
FEMaterialPoint* FEPreStrainCoupledTransIsoMR::CreateMaterialPointData()
{ 
	return new FEFiberPreStretchMaterialPoint(new FEElasticMaterialPoint); 
}

//-----------------------------------------------------------------------------
double FEPreStrainCoupledTransIsoMR::FiberStretch(FEMaterialPoint& mp)
{
	FEFiberPreStretchMaterialPoint& psp = *mp.ExtractData<FEFiberPreStretchMaterialPoint>();
	if (psp.m_ltrg ==0.0) return m_ltrg;
	else
	{
		double w = m_ltrg;
		double l = psp.m_ltrg;
		return w*l + (1.0-w);
	}
}

//-----------------------------------------------------------------------------
mat3d FEPreStrainCoupledTransIsoMR::PreStrainDeformationGradient(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberPreStretchMaterialPoint& psp = *mp.ExtractData<FEFiberPreStretchMaterialPoint>();

	// get the target stretch
	double ltrg = FiberStretch(mp);

	// set-up local uni-axial stretch tensor
	double l = psp.m_lam;
	double li = 1.0;
	mat3d U(l, 0.0, 0.0, 0.0, li, 0.0, 0.0, 0.0, li);

	mat3d Q = pt.m_Q;
	mat3d Qt = Q.transpose();

	mat3d F_bar = Q*U*Qt;

	return F_bar;
}

//-----------------------------------------------------------------------------
mat3ds FEPreStrainCoupledTransIsoMR::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// pre-multiply the pre-strain
	ep.m_F = F0*PreStrainDeformationGradient(mp);
	ep.m_J = ep.m_F.det();

	// evaluate the stress
	mat3ds s = m_mat.Stress(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEPreStrainCoupledTransIsoMR::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// pre-multiply the pre-strain
	ep.m_F = F0*PreStrainDeformationGradient(mp);
	ep.m_J = ep.m_F.det();

	// evaluate the tangent
	tens4ds c = m_mat.Tangent(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	// return spatial tangent
	return c;
}
