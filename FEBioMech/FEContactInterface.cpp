// FEContactInterface.cpp: implementation of the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEContactInterface.h"
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FESolver.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEModel* pfem) : FESurfacePairInteraction(pfem)
{
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
	FEMesh& m = GetFEModel()->GetMesh();

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// extract the elastic material
		FEElasticMaterial* pme = GetFEModel()->GetMaterial(pe->GetMatID())->GetElasticMaterial();

		if (pme)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

			// setup the material point
			pt.m_F = mat3dd(1.0);
			pt.m_J = 1;
			pt.m_s.zero();

			// get the tangent (stiffness) and it inverse (compliance) at this point
			tens4ds S = pme->Tangent(mp);
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
