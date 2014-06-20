// FERigidBody.cpp: implementation of the FERigidBody class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidBody.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidBody, FEObject);
	ADD_PARAMETER(m_Fr.x, FE_PARAM_DOUBLE, "Fx");
	ADD_PARAMETER(m_Fr.y, FE_PARAM_DOUBLE, "Fy");
	ADD_PARAMETER(m_Fr.z, FE_PARAM_DOUBLE, "Fz");
	ADD_PARAMETER(m_Mr.x, FE_PARAM_DOUBLE, "Mx");
	ADD_PARAMETER(m_Mr.y, FE_PARAM_DOUBLE, "My");
	ADD_PARAMETER(m_Mr.z, FE_PARAM_DOUBLE, "Mz");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidBody::FERigidBody(FEModel* pfem) : FEObject(pfem)
{
	m_bActive = true;
    m_bpofr = false;
	for (int i=0; i<6; ++i)
	{
		m_BC[i] = 0;
		m_pDC[i] = 0;
	}
	m_prb = 0;

	// zero total displacements
	m_Ut[0] = m_Up[0] = 0;
	m_Ut[1] = m_Up[1] = 0;
	m_Ut[2] = m_Up[2] = 0;
	m_Ut[3] = m_Up[3] = 0;
	m_Ut[4] = m_Up[4] = 0;
	m_Ut[5] = m_Up[5] = 0;

    // initialize velocity and acceleration of center of mass
    m_vt = m_at = vec3d(0,0,0);
    
	// initialize orientation
	m_qt = quatd(0, vec3d(0,0,1));
    
    // initialize angular velocity and acceleration
    m_wt = m_alt = vec3d(0,0,0);
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

    // initialize velocity and acceleration of center of mass
    m_vp = m_vt = vec3d(0,0,0);
    m_ap = m_at = vec3d(0,0,0);
    
	// initialize orientation
	m_qp = m_qt = quatd(0, vec3d(0,0,1));

    // initialize angular velocity and acceleration
    m_wp = m_wt = vec3d(0,0,0);
    m_alp = m_alt = vec3d(0,0,0);
    
	// initialize center of mass
	m_rt = m_r0;

	// reset reaction force and torque
	m_Fr = vec3d(0,0,0);
	m_Mr = vec3d(0,0,0);
    
    // if any of the rotation dofs is prescribed, all others should be
    // prescribed or fixed.
    m_bpofr = false;
    if (m_pDC[3] || m_pDC[4] || m_pDC[5])
    {
        bool bpofr[3] = {false};
        for (int j=3; j<6; ++j) if (m_pDC[j] || (m_LM[j] < 0)) bpofr[j-3] = true;
        if (bpofr[0] && bpofr[1] && bpofr[2]) m_bpofr = true;
        else
        {
            printf("FATAL ERROR: Rigid body rotations cannot mix prescribed and free components.\n");
            printf("Rigid body: %d, Material: %d\n",m_nID, GetMaterialID());
            throw "FATAL ERROR";
        }
    }
}

//-----------------------------------------------------------------------------
//! This function is called at the start of each time step and is used to update
//! some variables.
void FERigidBody::Init()
{
	// clear reaction forces
	m_Fr = m_Mr = vec3d(0,0,0);

	// store previous state
	m_rp = m_rt;
    m_vp = m_vt;
    m_ap = m_at;
	m_qp = m_qt;
    m_wp = m_wt;
    m_alp = m_alt;
	m_Up[0] = m_Ut[0];
	m_Up[1] = m_Ut[1];
	m_Up[2] = m_Ut[2];
	m_Up[3] = m_Ut[3];
	m_Up[4] = m_Ut[4];
	m_Up[5] = m_Ut[5];

	// zero incremental displacements
	m_du[0] = m_dul[0] = 0.0;
	m_du[1] = m_dul[1] = 0.0;
	m_du[2] = m_dul[2] = 0.0;
	m_du[3] = m_dul[3] = 0.0;
	m_du[4] = m_dul[4] = 0.0;
	m_du[5] = m_dul[5] = 0.0;
    
    // if any of the rotation dofs is prescribed, all others should be
    // prescribed or fixed.
    m_bpofr = false;
    if (m_pDC[3] || m_pDC[4] || m_pDC[5])
    {
        bool bpofr[3] = {false};
        for (int j=3; j<6; ++j) if (m_pDC[j] || (m_LM[j] < 0)) bpofr[j-3] = true;
        if (bpofr[0] && bpofr[1] && bpofr[2]) m_bpofr = true;
        else
        {
            printf("FATAL ERROR: Rigid body rotations cannot mix prescribed and free components.\n");
            printf("Rigid body: %d, Material: %d\n",m_nID, GetMaterialID());
            throw "FATAL ERROR";
        }
    }
}

//-----------------------------------------------------------------------------
//! Set the rigid body's center of mass directly
void FERigidBody::SetCOM(vec3d rc)
{
	m_r0 = m_rt = rc;
}

//-----------------------------------------------------------------------------
//! Calculates the rigid body's total mass, center of mass, and mass moment
//! of inertia about the center of mass
//!
void FERigidBody::UpdateCOM()
{
	// get the mesh
	FEMesh& mesh = m_fem.GetMesh();

	// initialize some data
	m_mass = 0;			// total mass of rigid body
	vec3d rc(0,0,0);	// center of mass
    mat3d moi(0,0,0,0,0,0,0,0,0);    // mass moment of inertia about origin
    mat3dd I(1);        // identity tensor

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
		if (mesh.Domain(nd).Class() == FE_DOMAIN_SOLID)
		{
			FESolidDomain* pbd = static_cast<FESolidDomain*>(&mesh.Domain(nd));
			FEMaterial* pm = pbd->GetMaterial();
			// make sure this element belongs to the rigid body
			if (pm->IsRigid() && (pm->GetRigidBodyID() == m_nID))
			{
				// get the material density
				double dens = pm->Density();
				assert(dens > 0.0);
				if (dens == 0.0) dens = 1.0;

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

						// add to com and moi
						for (int i=0; i<el.Nodes(); ++i)
						{
							rc += r0[i]*H[i]*detJ*gw[n]*dens;
                            for (int j=0; j<el.Nodes(); ++j) {
                                moi += ((r0[i]*r0[j])*I - (r0[i] & r0[j]))*H[i]*H[j]*detJ*gw[n]*dens;
                            }
						}
					}
				}
			}
		}
	}

	// normalize com
	if (m_mass != 0) rc /= m_mass;
    
    // use parallel axis theorem to transfer moi to com
    // and store moi
    m_moi = moi.sym() - m_mass*((rc*rc)*I - dyad(rc));

	// store com
	m_r0 = m_rt = rc;
}

//-----------------------------------------------------------------------------
void FERigidBody::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_mass;
        dmp << m_moi;
		dmp << m_Fr << m_Mr;
		dmp << m_rp << m_rt;
		dmp << m_vp << m_vt;
        dmp << m_ap << m_at;
		dmp << m_qp << m_qt;
		dmp << m_wp << m_wt;
		dmp << m_alp << m_alt;
		dmp << m_bActive;
		for (int i=0; i<6; ++i)
		{
			dmp << m_Up[i];
			dmp << m_Ut[i];
			dmp << m_du[i];
			dmp << m_dul[i];
		}
	}
	else
	{
		dmp >> m_mass;
        dmp >> m_moi;
		dmp >> m_Fr >> m_Mr;
		dmp >> m_rp >> m_rt;
		dmp >> m_vp >> m_vt;
        dmp >> m_ap >> m_at;
		dmp >> m_qp >> m_qt;
		dmp >> m_wp >> m_wt;
		dmp >> m_alp >> m_alt;
		dmp >> m_bActive;
		for (int i=0; i<6; ++i)
		{
			dmp >> m_Up[i];
			dmp >> m_Ut[i];
			dmp >> m_du[i];
			dmp >> m_dul[i];
		}
	}
}

//-----------------------------------------------------------------------------

void FERigidBody::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID << m_mat << m_mass << m_moi << m_Fr << m_Mr;
		ar << m_r0 << m_rt << m_rp << m_vt << m_vp << m_at << m_ap;
        ar << m_qt << m_qp << m_wt << m_wp << m_alt << m_alp;
		ar << m_bActive << m_bpofr;
		ar.write(m_LM , sizeof(int), 6);
		ar.write(m_Up , sizeof(double), 6);
		ar.write(m_Ut , sizeof(double), 6);
		ar.write(m_du , sizeof(double), 6);
		ar.write(m_dul, sizeof(double), 6);
	}
	else
	{
		ar >> m_nID >> m_mat >> m_mass >> m_moi >> m_Fr >> m_Mr;
		ar >> m_r0 >> m_rt >> m_rp >> m_vt >> m_vp >> m_at >> m_ap;
        ar >> m_qt >> m_qp >> m_wt >> m_wp >> m_alt >> m_alp;
		ar >> m_bActive >> m_bpofr;
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
				du[j] = (lc < 0? 0 : pdc->sf*m_fem.GetLoadCurve(lc)->Value() - m_Up[j] + pdc->ref);
			}
			else 
			{
				du[j] = (lm[j] >=0 ? Ui[lm[j]] + ui[lm[j]] : 0);
			}
		}
	}

	m_rt.x = m_rp.x + du[0];
	m_rt.y = m_rp.y + du[1];
	m_rt.z = m_rp.z + du[2];

	// next, we do the rotations. We do this seperatly since
	// they need to be interpreted differently than displacements
	if (m_prb == 0)
	{
        quatd qdu;          // quaternion of net increment
        
        if (m_bpofr) {
            // if all rotation components are known (prescribed or fixed)
            // evaluate net increment from load curve
            double Ut[3] = {0};
            for (int j=3; j<6; ++j) {
                if (m_pDC[j]) {
                    int lc = m_pDC[j]->lc;
                    Ut[j-3] = (lc < 0? 0 : m_pDC[j]->sf*m_fem.GetLoadCurve(lc)->Value());
                }
            }
            quatd qUt(vec3d(Ut[0],Ut[1],Ut[2]));
            qdu = qUt*m_qp.Inverse();
        }
        else
        {
            // rotation components are either free or fixed
            vec3d vUi(0,0,0);   // initialize total increment so far
            vec3d vui(0,0,0);   // initialize current increment
            if (lm[3] >= 0) { vUi.x = Ui[lm[3]]; vui.x = ui[lm[3]]; }
            if (lm[4] >= 0) { vUi.y = Ui[lm[4]]; vui.y = ui[lm[4]]; }
            if (lm[5] >= 0) { vUi.z = Ui[lm[5]]; vui.z = ui[lm[5]]; }
            quatd qUi(2*atan(vUi.norm()/2),vUi);                    // Cayley transform
            quatd qui(2*atan(vui.norm()/2),vui);                    // Cayley transform
            qdu = qui*qUi;
        }
        
        qdu.MakeUnit();                                         // clean-up roundoff errors
        vec3d vdu = qdu.GetVector()*(2*tan(qdu.GetAngle()/2));  // Cayley transform
        du[3] = vdu.x; du[4] = vdu.y; du[5] = vdu.z;
        
	}
    
    vec3d vdu(du[3],du[4],du[5]);
    quatd qdu(2*atan(vdu.norm()/2),vdu);
	m_qt = qdu*m_qp;     // update at the current time step
	m_qt.MakeUnit();
    
	if (m_prb) du = m_dul;
    // update RB center of mass translations
	m_Ut[0] = m_Up[0] + du[0];
	m_Ut[1] = m_Up[1] + du[1];
	m_Ut[2] = m_Up[2] + du[2];
    // update RB rotations
    vec3d vUt = m_qt.GetVector()*m_qt.GetAngle();
    m_Ut[3] = vUt.x;
    m_Ut[4] = vUt.y;
    m_Ut[5] = vUt.z;
    
    
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

//-----------------------------------------------------------------------------
void FERigidBody::UpdateIncrement(std::vector<double>& Ui, std::vector<double>& ui, bool emap)
{
	if (m_prb == 0)
	{
        int *lm = m_LM;
        vec3d v;
        quatd qUi;
        
        // first do the displacements
		for (int j=0; j<3; ++j)
            if (lm[j] >=0) Ui[lm[j]] += ui[lm[j]];
        
        // next, we do the rotations. We do this separately since
        // they need to be interpreted differently than displacements
        vec3d vUi(0,0,0);
        vec3d vui(0,0,0);
        if (lm[3] >= 0) { vUi.x = Ui[lm[3]]; vui.x = ui[lm[3]]; }
        if (lm[4] >= 0) { vUi.y = Ui[lm[4]]; vui.y = ui[lm[4]]; }
        if (lm[5] >= 0) { vUi.z = Ui[lm[5]]; vui.z = ui[lm[5]]; }
        if (emap) qUi = quatd(vUi);
        else qUi = quatd(2*atan(vUi.norm()/2),vUi);     // Cayley transform
        quatd qui(2*atan(vui.norm()/2),vui);            // Cayley transform
        quatd q = qui*qUi;
        q.MakeUnit();
        if (emap) v = q.GetVector()*q.GetAngle();
        else v = q.GetVector()*(2*tan(q.GetAngle()/2)); // Cayley transform
        if (lm[3] >= 0) Ui[lm[3]] = v.x;
        if (lm[4] >= 0) Ui[lm[4]] = v.y;
        if (lm[5] >= 0) Ui[lm[5]] = v.z;
	}
}
