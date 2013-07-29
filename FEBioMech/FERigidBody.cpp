// FERigidBody.cpp: implementation of the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidBody.h"
#include "FERigid.h"
#include "FECore/FEMaterial.h"
#include "FECore/FESolidDomain.h"

//-----------------------------------------------------------------------------
FERigidBody::FERigidBody(FEModel* pfem) : FEObject(pfem)
{
	m_bActive = true;
	for (int i=0; i<6; ++i) m_pDC[i] = 0;
	m_prb = 0;

	// zero total displacements
	m_Ut[0] = m_Up[0] = 0;
	m_Ut[1] = m_Up[1] = 0;
	m_Ut[2] = m_Up[2] = 0;
	m_Ut[3] = m_Up[3] = 0;
	m_Ut[4] = m_Up[4] = 0;
	m_Ut[5] = m_Up[5] = 0;

	// initialize orientation
	m_qt = quatd(0, vec3d(0,0,1));
}

//-----------------------------------------------------------------------------
FERigidBody::~FERigidBody()
{

}

//-----------------------------------------------------------------------------
//! Reset rigid body data (called from FEM::Reset)
void FERigidBody::Reset()
{
	// zero total displacements
	m_Ut[0] = m_Up[0] = 0;
	m_Ut[1] = m_Up[1] = 0;
	m_Ut[2] = m_Up[2] = 0;
	m_Ut[3] = m_Up[3] = 0;
	m_Ut[4] = m_Up[4] = 0;
	m_Ut[5] = m_Up[5] = 0;

	// initialize orientation
	m_qt = quatd(0, vec3d(0,0,1));

	// initialize center of mass
	m_rt = m_r0;

	// reset reaction force and torque
	m_Fr = vec3d(0,0,0);
	m_Mr = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
//! Calculates the total mass and center of mass of a rigid body
//!
void FERigidBody::UpdateCOM()
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// initialize some data
	m_mass = 0;			// total mass of rigid body
	vec3d rc(0,0,0);	// center of mass

	// jacobian
	double detJ;

	// shape function values
	double* H;

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	
	// loop over all elements
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		// TODO: I should convert to a FERigidSolidDomain or FERigidShellDomain
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (pbd->GetMaterial());
			// make sure this element belongs to the rigid body
			if (pm && (pm->m_nRB == m_nID))
			{
				// get the material density
				double dens = pm->Density();

				// loop over all elements
				for (int iel=0; iel<pbd->Elements(); ++iel)
				{	
					FESolidElement& el = pbd->Element(iel);

					// nr of integration points
					int nint = el.GaussPoints();

					// number of nodes
					int neln = el.Nodes();
					assert(neln <= 8);

					// initial coordinates
					for (int i=0; i<neln; ++i) r0[i] = pbd->GetMesh()->Node(el.m_node[i]).m_r0;

					// integration weights
					double* gw = el.GaussWeights();

					// loop over integration points
					for (int n=0; n<nint; ++n)
					{
						// calculate jacobian
						detJ = pbd->detJ0(el, n);

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
void FERigidBody::ShallowCopy(FEObject *po)
{
	FERigidBody& rb = dynamic_cast<FERigidBody&>(*po);
	m_mass = rb.m_mass;
	m_Fr = rb.m_Fr;
	m_Mr = rb.m_Mr;

	m_rp = rb.m_rp;
	m_rt = rb.m_rt;

	m_qp = rb.m_qp;
	m_qt = rb.m_qt;

	m_bActive = rb.m_bActive;

	for (int i=0; i<6; ++i)
	{
		m_Up[i] = rb.m_Up[i];
		m_Ut[i] = rb.m_Ut[i];
		m_du[i] = rb.m_du[i];
		m_dul[i] = rb.m_dul[i];
		m_pDC[i] = rb.m_pDC[i];
	}
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

//-----------------------------------------------------------------------------
void FERigidBody::Update(std::vector<double>& Ui, std::vector<double>& ui)
{
	int *lm = m_LM;
	double* du = m_du;

	// first do the displacements
	if (m_prb == 0)
	{
		FERigidBodyDisplacement* pdc;
		for (int j=0; j<3; ++j)
		{
			pdc = m_pDC[j];
			if (pdc)
			{
				int lc = pdc->lc;
				// TODO: do I need to take the line search step into account here?
				du[j] = (lc < 0? 0 : pdc->sf*m_fem.GetLoadCurve(lc)->Value() - m_Up[j]);
			}
			else du[j] = (lm[j] >=0 ? Ui[lm[j]] + ui[lm[j]] : 0);
		}
	}

	m_rt.x = m_rp.x + du[0];
	m_rt.y = m_rp.y + du[1];
	m_rt.z = m_rp.z + du[2];

	// next, we do the rotations. We do this seperatly since
	// they need to be interpreted differently than displacements
	if (m_prb == 0)
	{
		FERigidBodyDisplacement* pdc;
		for (int j=3; j<6; ++j)
		{
			pdc = m_pDC[j];
			if (pdc)
			{
				int lc = pdc->lc;
				// TODO: do I need to take the line search step into account here?
				du[j] = (lc < 0? 0 : pdc->sf*m_fem.GetLoadCurve(lc)->Value() - m_Up[j]);
			}
			else du[j] = (lm[j] >=0 ? Ui[lm[j]] + ui[lm[j]] : 0);
		}
	}

	vec3d r = vec3d(du[3], du[4], du[5]);
	double w = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
	quatd dq = quatd(w, r);

	m_qt = dq*m_qp;
	m_qt.MakeUnit();

	if (m_prb) du = m_dul;
	m_Ut[0] = m_Up[0] + du[0];
	m_Ut[1] = m_Up[1] + du[1];
	m_Ut[2] = m_Up[2] + du[2];
	m_Ut[3] = m_Up[3] + du[3];
	m_Ut[4] = m_Up[4] + du[4];
	m_Ut[5] = m_Up[5] + du[5];

	// update the mesh' nodes
	FEMesh& mesh = m_fem.GetMesh();
	int N = mesh.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid == m_nID)
		{
			vec3d a0 = node.m_r0 - m_r0;
			vec3d at = m_qt*a0;
			node.m_rt = m_rt + at;
		}
	}
}
