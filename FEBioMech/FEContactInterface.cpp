#include "stdafx.h"
#include "FEContactInterface.h"
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FESolver.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEModel* pfem) : FESurfacePairConstraint(pfem)
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
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();

	// get the element this surface element belongs to
	FEElement* pe = el.m_elem[0];
	if (pe == nullptr) return 0.0;

	// extract the elastic material
	FEElasticMaterial* pme = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return 0.0;

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
	double eps = 1./(n*(vdotTdotv(n, C, n)*n));

	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

	return eps*A/V;
}

//-----------------------------------------------------------------------------
void FEContactInterface::Serialize(DumpStream& ar)
{
	// store base class
	FESurfacePairConstraint::Serialize(ar);

	// save parameters
	ar & m_blaugon;
}
