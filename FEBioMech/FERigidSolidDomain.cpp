#include "stdafx.h"
#include "FERigidSolidDomain.h"
#include "FERigidMaterial.h"
#include <FECore/FERigidBody.h>

//-----------------------------------------------------------------------------
bool FERigidSolidDomain::Initialize(FEModel& fem)
{
	return FESolidDomain::Initialize(fem);
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidSolidDomain::Reset()
{
	// nothing to reset here
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements.
//! Rigid elements don't generate stress, so there is nothing to do here
void FERigidSolidDomain::StiffnessMatrix(FESolver* psolver)
{
	// Caught you looking!
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FESolver* psolver, vector<double>& R)
{
	// what you looking at ?!
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}


//-----------------------------------------------------------------------------
// Calculate inertial forces
void FERigidSolidDomain::InertialForces(FEGlobalVector& R, vector<double>& Q)
{
	// This is calculated elsewhere
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::MassMatrix(FESolver* psolver, double scale)
{
	// this is calculated elsewhere
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    FERigidMaterial* pmr = dynamic_cast<FERigidMaterial*>(GetMaterial());
	FEModel& fem = *pmr->GetFEModel();
	FESolver* psolver = fem.GetCurrentStep()->GetFESolver();
    
    // get the rigid body
    FERigidBody& RB = dynamic_cast<FERigidBody&>(*fem.Object(pmr->GetRigidBodyID()));
        
    // 6 dofs per rigid body
    double fe[6];
        
    // evaluate body force per mass
    FEElasticMaterialPoint mp;
    mp.m_r0 = RB.m_r0;
    // use alpha rule
    double alpha = psolver->m_alpha;
    mp.m_rt = RB.m_rt*alpha + RB.m_rp*(1-alpha);
    vec3d b = BF.force(mp);
        
    // body force = mass*(body force per mass)
    vec3d F = b*RB.m_mass;
        
    fe[0] = -F.x;
    fe[1] = -F.y;
    fe[2] = -F.z;
        
    // moment of body force about center of mass is zero
    vec3d M(0,0,0);
        
    fe[3] = -M.x;
    fe[4] = -M.y;
    fe[5] = -M.z;
        
    R.Assemble(RB.m_LM, fe, 6);
        
    // add to rigid body force
    RB.m_Fr += F;
        
    // add to rigid body torque
    RB.m_Mr += M;
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
	FEModel& fem = psolver->GetFEModel();
    
    FERigidMaterial* pmr = dynamic_cast<FERigidMaterial*>(GetMaterial()); assert(pmr);
    
    vector<int> lm;
        
    // get the rigid body
    FERigidBody& RB = dynamic_cast<FERigidBody&>(*fem.Object(pmr->GetRigidBodyID()));
        
    // 6 dofs per rigid body
    matrix ke(6, 6);
    ke.zero();
        
    // evaluate body force stiffness per mass
    FEElasticMaterialPoint mp;
    mp.m_r0 = RB.m_r0;
    // use alpha rule
    double alpha = psolver->m_alpha;
    mp.m_rt = RB.m_rt*alpha + RB.m_rp*(1-alpha);
    mat3ds k = bf.stiffness(mp);
        
    // body force stiffness = mass*(body force stiffness per mass)
    // multiply by alpha because of alpha rule
    mat3ds K = k*(RB.m_mass*alpha);
        
    // since moment of body force about center of mass is zero
    // there is no moment contribution to the stiffness
    ke[0][0] = K.xx(); ke[0][1] = K.xy(); ke[0][2] = K.xz();
    ke[1][0] = K.xy(); ke[1][1] = K.yy(); ke[1][2] = K.yz();
    ke[2][0] = K.xz(); ke[2][1] = K.yz(); ke[2][2] = K.zz();
        
    lm.assign(RB.m_LM, RB.m_LM+6);
        
    psolver->AssembleStiffness(lm, ke);
}
