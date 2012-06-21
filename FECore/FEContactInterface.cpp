// FEContactInterface.cpp: implementation of the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEContactInterface.h"
#include "FEModel.h"
#include "FEElasticMaterial.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEModel* pfem)
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
double FEContactInterface::BulkModulus(FESurfaceElement& el, FESurface &s)
{
	double eps = 0;

	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// extract the elastic material
		FEElasticMaterial* pme = m_pfem->GetMaterial(pe->GetMatID())->GetElasticMaterial();

		// get a material point
		FEMaterialPoint& mp = *pe->m_State[0];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// setup the material point
		pt.F = mat3dd(1.0);
		pt.J = 1;
		pt.s.zero();

		// get the tangent at this point
		tens4ds C = pme->Tangent(pt);

		// get the upper 3x3
		mat3d T(mat3ds(C(0,0), C(1,1), C(2,2), C(0,1), C(1,2), C(0,2)));

		// invert the matrix
		mat3d Ti = T.inverse();

		// calculate average modulus
		eps = (1/Ti(0,0) + 1/Ti(1,1) + 1/Ti(2,2))/3;
	}

	return eps;
}

//-----------------------------------------------------------------------------
void FEContactInterface::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_blaugon;
	}
	else
	{
		ar >> m_nID;
		ar >> m_blaugon;
	}

	// store parameters
	FEParamContainer::Serialize(ar);
}
