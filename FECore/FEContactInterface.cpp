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

	m_bactive = true;
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
	FEMesh& m = m_pfem->GetMesh();

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// extract the elastic material
		FEElasticMaterial* pme = m_pfem->GetMaterial(pe->GetMatID())->GetElasticMaterial();

		if (pme)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

			// setup the material point
			pt.F = mat3dd(1.0);
			pt.J = 1;
			pt.s.zero();

			// get the tangent (stiffness) and it inverse (compliance) at this point
			tens4ds S = pme->Tangent(pt);
			tens4ds C = S.inverse();

			// evaluate element surface normal at parametric center
			vec3d t[2];
			s.CoBaseVectors0(el, 0, 0, t);
			vec3d n = t[0] ^ t[1];
			n.unit();
		
			// evaluate normal component of the compliance matrix
			// (equivalent to inverse of Young's modulus along n)
			eps = 1./(n*(vdotTdotv(n, C, n)*n));
		}
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
