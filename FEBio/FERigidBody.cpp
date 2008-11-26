// FERigidBody.cpp: implementation of the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidBody.h"
#include "FEMaterial.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FERigidBody::FERigidBody()
{
	m_pfem = 0;
}

FERigidBody::~FERigidBody()
{

}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEM::Update
// Calculates the total mass and center of mass of a rigid body
//

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
	// TODO: what if the rigid body has shells ?
	for (int iel=0; iel<mesh.SolidElements(); ++iel)
	{
		FESolidElement& el = mesh.SolidElement(iel);

		FERigid* pm = dynamic_cast<FERigid*> (fem.GetMaterial(el.GetMatID()));

		// make sure this element belongs to the rigid body
		if (pm && (pm->m_nRB == m_nID))
		{
			dens = pm->m_density;

			// unpack the element
			mesh.UnpackElement(el);

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

	// normalize com
	if (m_mass != 0) rc /= m_mass;

	// store com
	m_r0 = m_rt = rc;
}
