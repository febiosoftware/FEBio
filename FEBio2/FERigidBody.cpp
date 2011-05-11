// FERigidBody.cpp: implementation of the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "fem.h"
#include "FERigidBody.h"
#include "FECore/FEMaterial.h"
#include "FERigid.h"
#include "FEElasticSolidDomain.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FERigidBody::FERigidBody()
{
	m_pfem = 0;
	m_bActive = true;
	for (int i=0; i<6; ++i) m_pDC[i] = 0;
	m_prb = 0;
}

FERigidBody::~FERigidBody()
{

}

//-----------------------------------------------------------------------------
//! Calculates the total mass and center of mass of a rigid body
//!
void FERigidBody::Update()
{
	// make sure the rigid body is attached to a FEM
	if (m_pfem == 0) return;

	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	m_mass = 0;			// total mass of rigid body
	vec3d rc(0,0,0);	// center of mass

	// jacobian
	double detJ;

	// shape function values
	double* H;

	// material density
	double dens;

	
	// loop over all elements
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (int iel=0; iel<pbd->Elements(); ++iel)
			{	
				FESolidElement& el = pbd->Element(iel);

				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (fem.GetMaterial(el.GetMatID()));

				// make sure this element belongs to the rigid body
				if (pm && (pm->m_nRB == m_nID))
				{
					dens = pm->m_density;

					// unpack the element
					pbd->UnpackElement(el);

					// nr of integration points
					int nint = el.GaussPoints();

					// initial coordinates
					vec3d* r0 = el.r0();

					// integration weights
					double* gw = el.GaussWeights();

					// loop over integration points
					for (int n=0; n<nint; ++n)
					{
						// calculate jacobian
						detJ = el.detJ0(n);

						// shape functions at integration point
						H = el.H(n);

						// add to total mass
						m_mass += dens*detJ*gw[n];

						// add to com
						for (int i=0; i<el.Nodes(); ++i)
						{
							rc += r0[i]*H[i]*detJ*gw[n]*dens;
						}
					}
				}
			}
		}
	}

	// normalize com
	if (m_mass != 0) rc /= m_mass;

	// store com
	m_r0 = m_rt = rc;
}

//-----------------------------------------------------------------------------

void FERigidBody::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID << m_mat << m_mass << m_Fr << m_Mr;
		ar << m_r0 << m_rt << m_rp << m_qt << m_qp;
		ar << m_bActive;
		ar.write(m_LM , sizeof(int), 6);
		ar.write(m_Up , sizeof(double), 6);
		ar.write(m_Ut , sizeof(double), 6);
		ar.write(m_du , sizeof(double), 6);
		ar.write(m_dul, sizeof(double), 6);
	}
	else
	{
		ar >> m_nID >> m_mat >> m_mass >> m_Fr >> m_Mr;
		ar >> m_r0 >> m_rt >> m_rp >> m_qt >> m_qp;
		ar >> m_bActive;
		ar.read(m_LM , sizeof(int   ), 6);
		ar.read(m_Up , sizeof(double), 6);
		ar.read(m_Ut , sizeof(double), 6);
		ar.read(m_du , sizeof(double), 6);
		ar.read(m_dul, sizeof(double), 6);
	}
}
