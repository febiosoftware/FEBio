// FEContactInterface.cpp: implementation of the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEContactInterface.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEM* pfem)
{
	m_pfem = pfem;
	m_ntype = 0;
	m_nID = -1;

	m_blaugon = false;
}

FEContactInterface::~FEContactInterface()
{

}

//-----------------------------------------------------------------------------
//! This function calculates a contact penalty parameter based on the 
//! material and geometrical properties of the slave and master surfaces
//!
double FEContactInterface::AutoPenalty(FESurfaceElement& el, FESurface &s)
{
	double eps = 0;

	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// get the material
		FESolidMaterial* pm = dynamic_cast<FESolidMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));

		// extract the elastic component
		FEElasticMaterial* pme = m_pfem->GetElasticMaterial(pe->GetMatID());

		// get a material point
		FEMaterialPoint& mp = *pe->m_State[0];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// setup the material point
		pt.F = mat3dd(1.0);
		pt.J = 1;
		pt.avgJ = 1;
		pt.avgp = 0;
		pt.s.zero();

		// get the tangent at this point
		tens4ds C = pme->Tangent(pt);

		// get the upper 3x3
		mat3d T(mat3ds(C(0,0), C(1,1), C(2,2), C(0,1), C(1,2), C(0,2)));

		// see if the material is incompressible
		FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pme);
		if (pmi)
		{
			// if the material is incompressible, the elasticity tensor returned by Tangent()
			// does not contain the dilatational term, se we need to it.
			// This term is simply D2U/DJ2(J=1)1x1
			double K = pmi->Upp(1);
			T += mat3ds(K, K, K, K, K, K);
		}

		// invert the matrix
		mat3d Ti = T.inverse();

		// calculate average modulus
		eps = (1/Ti(0,0) + 1/Ti(1,1) + 1/Ti(2,2))/3;
	}

	return eps;
}
