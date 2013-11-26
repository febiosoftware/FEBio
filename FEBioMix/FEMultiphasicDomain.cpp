#include "FEMultiphasicDomain.h"
#include "FEMultiphasic.h"
#include "FECore/log.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FEMultiphasicDomain::FEMultiphasicDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat)
{ 
	m_ntype = FE_MULTIPHASIC_DOMAIN; 
}

//-----------------------------------------------------------------------------
//! clone domain
FEDomain* FEMultiphasicDomain::Clone()
{
	FEMultiphasicDomain* pd = new FEMultiphasicDomain(m_pMesh, m_pMat);
	pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
	return pd;
}

//-----------------------------------------------------------------------------
bool FEMultiphasicDomain::Initialize(FEModel &mdl)
{
	// initialize base class
	FEElasticSolidDomain::Initialize(mdl);

	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());
		
	// get the multiphasic material
	FEMultiphasic* pmb = dynamic_cast<FEMultiphasic*>(pm); assert(pmb);
	const int nsol = pmb->Solutes();
	const int nsbm = pmb->SBMs();

	// extract the initial concentrations of the solid-bound molecules
	vector<double> sbmr(nsbm,0);
	for (int i=0; i<nsbm; ++i) {
		sbmr[i] = pmb->GetSBM(i)->m_rho0;
	}

	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
			
			// initialize multiphasic solutes
			ps.m_nsol = nsol;
			ps.m_c.assign(nsol,0);
			ps.m_ca.assign(nsol,0);
			ps.m_gradc.assign(nsol,0);
			ps.m_k.assign(nsol, 0);
			ps.m_dkdJ.assign(nsol, 0);
			ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
			ps.m_j.assign(nsol,0);
			ps.m_nsbm = nsbm;
			ps.m_sbmr = sbmr;
			ps.m_sbmrp = sbmr;
			ps.m_sbmrhat.assign(nsbm,0);
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::Reset()
{
	// reset base class
	FEElasticSolidDomain::Reset();
	
	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());
    
	// get the multiphasic material
	FEMultiphasic* pmb = dynamic_cast<FEMultiphasic*>(pm); assert(pmb);
	const int nsol = pmb->Solutes();
	const int nsbm = pmb->SBMs();
	
	// extract the initial concentrations of the solid-bound molecules
	vector<double> sbmr(nsbm,0);
	for (int i=0; i<nsbm; ++i) {
		sbmr[i] = pmb->GetSBM(i)->m_rho0;
	}
	
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
			
			// initialize multiphasic solutes
			ps.m_nsol = nsol;
			ps.m_c.assign(nsol,0);
			ps.m_ca.assign(nsol,0);
			ps.m_gradc.assign(nsol,0);
			ps.m_k.assign(nsol, 0);
			ps.m_dkdJ.assign(nsol, 0);
			ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
			ps.m_j.assign(nsol,0);
			ps.m_nsbm = nsbm;
			ps.m_sbmr = sbmr;
			ps.m_sbmrp = sbmr;
			ps.m_sbmrhat.assign(nsbm,0);
		}
	}
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InitElements()
{
	FEElasticSolidDomain::InitElements();
	
	// store previous mesh state
	// we need it for receptor-ligand complex calculations
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// reset referential solid volume fraction at previous time
			pt.m_phi0p = pt.m_phi0;

			// reset referential solid-bound molecule concentrations at previous time
			for (int j=0; j<ps.m_nsbm; ++j) {
				ps.m_sbmrp[j] = ps.m_sbmr[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
/*
void FEMultiphasicDomain::Residual(FESolidSolver* psolver, vector<double>& R)
{
	int i, j, dofc;
	
	FEM& fem = psolver->m_fem;
	
	// make sure we are in multiphasic mode
	assert(fem.m_pStep->m_nModule == FE_POROSOLUTE);
	
	// element force vector
	vector<double> fe;
	
	vector<int> elm;
	
	int NE = m_Elem.size();
	if (fem.m_pStep->m_nanalysis == FE_STEADY_STATE) {
		for (i=0; i<NE; ++i)
		{
			// get the element
			FESolidElement& el = m_Elem[i];
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the elements material
			FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (fem.GetMaterial(el.GetMatID()));
			assert(pm);
			const int nsol = pm->m_pSolute.size();
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::InternalForces(el, fe);
			
			// apply body forces
			// TODO: can we calculate body-forces with our formulation
			//       of biphasic theory
			//
			//if (fem.UseBodyForces())
			// {
			// BodyForces(fem, el, fe);
			// }
						
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			// calculate fluid internal work
			InternalFluidWorkSS(fem, el, fe);
			
			// add fluid work to global residual
			int neln = el.Nodes();
			int J;
			for (j=0; j<neln; ++j)
			{
				J = elm[3*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			for (int isol=0; isol<nsol; ++isol) {
				// calculate solute internal work
				InternalSoluteWorkSS(fem, el, fe, isol);
				
				// add solute work to global residual
				dofc = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
				for (j=0; j<neln; ++j)
				{
					J = elm[dofc*neln+j];
					if (J >= 0) R[J] += fe[j];
				}
			}
		}
	} else {
		for (i=0; i<NE; ++i)
		{
			// get the element
			FESolidElement& el = m_Elem[i];
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the elements material
			FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (fem.GetMaterial(el.GetMatID()));
			assert(pm);
			const int nsol = pm->m_pSolute.size();
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::InternalForces(el, fe);
			
			// apply body forces
			// TODO: can we calculate body-forces with our formulation
			//       of biphasic theory
			//
			// if (fem.UseBodyForces())
			// {
			// BodyForces(fem, el, fe);
			// }
			//
			
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			// calculate fluid internal work
			InternalFluidWork(fem, el, fe);
			
			// add fluid work to global residual
			int neln = el.Nodes();
			int J;
			for (j=0; j<neln; ++j)
			{
				J = elm[3*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			for (int isol=0; isol<nsol; ++isol) {
				// calculate solute internal work
				InternalSoluteWork(fem, el, fe, isol);
				
				// add solute work to global residual
				dofc = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
				for (j=0; j<neln; ++j)
				{
					J = elm[dofc*neln+j];
					if (J >= 0) R[J] += fe[j];
				}
				
			}
		}
	}
}
*/

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the solute work
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes
void FEMultiphasicDomain::InternalSoluteWorkSS(FESolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	const int nsol = pm->Solutes();

	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int neln = el.Nodes();
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
			
		for (int isol=0; isol<nsol; ++isol) 
		{
			// calculate solute internal work
			ElementInternalSoluteWorkSS(el, fe, dt, isol);
				
			// add solute work to global residual
			int dofc = DOF_C + pm->GetSolute(isol)->GetSoluteID();
			for (int j=0; j<neln; ++j)
			{
				int J = elm[dofc*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEMultiphasicDomain::ElementInternalSoluteWorkSS(FESolidElement& el, vector<double>& fe, double dt, const int sol)
{
	int i, isol, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;
	
	// gradient of shape functions
	vector<vec3d> gradN(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);

	const int nsol = pm->Solutes();
	vector<int> z(nsol);
	
	const int nreact = pm->Reactions();
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			// save spatial gradient of shape functions
			gradN[i] = vec3d(Gx,Gy,Gz);
		}
		
		vector<vec3d> j(spt.m_j);
		vec3d je(0);
		
		for (isol=0; isol<nsol; ++isol) {
			// get the charge number
			z[isol] = pm->GetSolute(isol)->ChargeNumber();
			// current density (flux units)
			je += j[isol]*z[isol];
		}

		// chemical reactions
		double phiw = pm->Porosity(mp);
		double chat = 0;
		for (i=0; i<nreact; ++i)
			chat += pm->GetReaction(i)->m_v[sol]*pm->GetReaction(i)->ReactionSupply(mp);

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[sol]+je*pm->m_penalty)
						 + H[i]*phiw*chat
						 )*detJ*wg[n];
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalSoluteWork(FESolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	const int nsol = pm->Solutes();

	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int neln = el.Nodes();
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
			
		for (int isol=0; isol<nsol; ++isol) 
		{
			// calculate solute internal work
			ElementInternalSoluteWork(el, fe, dt, isol);
				
			// add solute work to global residual
			int dofc = DOF_C + pm->GetSolute(isol)->GetSoluteID();
			for (int j=0; j<neln; ++j)
			{
				int J = elm[dofc*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes
bool FEMultiphasicDomain::ElementInternalSoluteWork(FESolidElement& el, vector<double>& fe, double dt, const int sol)
{
	int i, isol, jsol, n, isbm;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// gradient of shape functions
	vector<vec3d> gradN(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	FEMesh& mesh = *m_pMesh;
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);

	const int nsol = pm->Solutes();
    const int nsbm = pm->SBMs();

	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->GetSolute(isol)->GetSoluteID();
	
	const int nreact = pm->Reactions();
	
	const int NE = FEElement::MAX_NODES;
	vec3d r0[NE], rt[NE], rp[NE], vt[NE];
	vector< vector<double> > cp(nsol, vector<double>(8));
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		vt[i] = mesh.Node(el.m_node[i]).m_vt;
		for (isol=0; isol<nsol; ++isol)
			cp[isol][i] = mesh.Node(el.m_node[i]).m_cp[sid[isol]];
	}
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
		vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
		vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		invjac0(el, J0i, n);
		vec3d G1(J0i[0][0],J0i[0][1],J0i[0][2]);
		vec3d G2(J0i[1][0],J0i[1][1],J0i[1][2]);
		vec3d G3(J0i[2][0],J0i[2][1],J0i[2][2]);
		
		// next we calculate the deformation gradient and the solid velocity
		mat3d Fp;
		Fp.zero();
		vec3d vs(0);
		vector<double> cprev(nsol,0);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			GX = J0i[0][0]*Gr[i]+J0i[1][0]*Gs[i]+J0i[2][0]*Gt[i];
			GY = J0i[0][1]*Gr[i]+J0i[1][1]*Gs[i]+J0i[2][1]*Gt[i];
			GZ = J0i[0][2]*Gr[i]+J0i[1][2]*Gs[i]+J0i[2][2]*Gt[i];
			
			Fp[0][0] += rp[i].x*GX; Fp[1][0] += rp[i].y*GX; Fp[2][0] += rp[i].z*GX;
			Fp[0][1] += rp[i].x*GY; Fp[1][1] += rp[i].y*GY; Fp[2][1] += rp[i].z*GY;
			Fp[0][2] += rp[i].x*GZ; Fp[1][2] += rp[i].y*GZ; Fp[2][2] += rp[i].z*GZ;
			
			// calculate solid velocity
			vs += vt[i]*H[i];
			
			// save spatial gradient of shape functions
			gradN[i] = vec3d(Gx,Gy,Gz);
			
			// calculate effective concentration at previous time step
			for (isol=0; isol<nsol; ++isol)
				cprev[isol] += cp[isol][i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.m_J;
		double dJdt = (J-Jp)/dt;
		
		vector<vec3d> j(spt.m_j);
		vector<double> c(spt.m_c);
		vector<int> z(nsol);
		vector<double> kappa(spt.m_k);
		vec3d je(0);
		vector<double> dcdt(nsol);
		
		for (isol=0; isol<nsol; ++isol) {
			// get the time derivative of the effective concentration
			dcdt[isol] = (c[isol] - cprev[isol])/dt;
			// get the charge number
			z[isol] = pm->GetSolute(isol)->ChargeNumber();
			je += j[isol]*z[isol];
		}
		
		double dkdJ = spt.m_dkdJ[sol];
		vector<double> dkdc(spt.m_dkdc[sol]);
		vector<double> dkdr(spt.m_dkdr[sol]);
		
		// evaluate the porosity, its derivative w.r.t. J, and its gradient
		double phiw = pm->Porosity(mp);
		// evaluate time derivatives of solubility and porosity
		double dkdt = dkdJ*dJdt;
		for (jsol=0; jsol<nsol; ++jsol)
			dkdt += dkdc[jsol]*dcdt[jsol];
        double dphi0dt = 0;
        double chat = 0;
		
		// chemical reactions
        for (i=0; i<nreact; ++i) {
			FEChemicalReaction* pri = pm->GetReaction(i);
            double zhat = pri->ReactionSupply(mp);
            chat += phiw*zhat*pri->m_v[sol];
            for (isbm=0; isbm<nsbm; ++isbm) {
                double M = pm->SBMMolarMass(isbm);
                double rhoT = pm->SBMDensity(isbm);
                double v = pri->m_v[nsol+isbm];
                dkdt += J*phiw*zhat*M*v*dkdr[isbm];
                dphi0dt += J*phiw*zhat*v*M/rhoT;
            }
        }

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[sol]+je*pm->m_penalty)
						 + H[i]*(chat - (dJdt-dphi0dt)*kappa[sol]*c[sol]/J
                                 - phiw*dkdt*c[sol]
                                 - phiw*kappa[sol]*dcdt[sol]
                                 )
						 )*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalFluidWork(FESolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;
	
	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
		
		// calculate fluid internal work
		ElementInternalFluidWork(el, fe, dt);
			
		// add fluid work to global residual
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j)
		{
			int J = elm[3*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FEMultiphasicDomain::ElementInternalFluidWork(FESolidElement& el, vector<double>& fe, double dt)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	FEMesh& mesh = *m_pMesh;
	
	vec3d rp[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
	}
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);

	const int nreact = pm->Reactions();
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(el.m_State[n]->ExtractData<FEBiphasicMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		invjac0(el, J0i, n);
		
		// next we calculate the deformation gradient
		mat3d Fp;
		Fp.zero();
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			GX = J0i[0][0]*Gr[i]+J0i[1][0]*Gs[i]+J0i[2][0]*Gt[i];
			GY = J0i[0][1]*Gr[i]+J0i[1][1]*Gs[i]+J0i[2][1]*Gt[i];
			GZ = J0i[0][2]*Gr[i]+J0i[1][2]*Gs[i]+J0i[2][2]*Gt[i];
			
			Fp[0][0] += rp[i].x*GX; Fp[1][0] += rp[i].y*GX; Fp[2][0] += rp[i].z*GX;
			Fp[0][1] += rp[i].x*GY; Fp[1][1] += rp[i].y*GY; Fp[2][1] += rp[i].z*GY;
			Fp[0][2] += rp[i].x*GZ; Fp[1][2] += rp[i].y*GZ; Fp[2][2] += rp[i].z*GZ;
			
			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.m_J;
		
		// and then finally
		double divv = ((J-Jp)/dt)/J;
		
		// get the flux
		vec3d& w = ppt.m_w;

		// get the solvent supply
		double phiwhat = 0;
		if (pm->GetSolventSupply()) phiwhat = pm->GetSolventSupply()->Supply(mp);
	
		// chemical reactions
        double phiw = pm->Porosity(mp);
		for (i=0; i<nreact; ++i)
			phiwhat += phiw*pm->GetReaction(i)->m_Vbar*pm->GetReaction(i)->ReactionSupply(mp);

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z + (phiwhat - divv)*H[i])*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalFluidWorkSS(FESolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;
	
	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
		
		// calculate fluid internal work
		ElementInternalFluidWorkSS(el, fe, dt);
			
		// add fluid work to global residual
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j)
		{
			int J = elm[3*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! for a steady-state analysis (zero solid velocity)
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FEMultiphasicDomain::ElementInternalFluidWorkSS(FESolidElement& el, vector<double>& fe, double dt)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);

	const int nreact = pm->Reactions();

	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEBiphasicMaterialPoint& ppt = *(el.m_State[n]->ExtractData<FEBiphasicMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
		}
		
		// get the solvent flux
		vec3d& w = ppt.m_w;

		// get the solvent supply
		double phiwhat = 0;
		if (pm->GetSolventSupply()) phiwhat = pm->GetSolventSupply()->Supply(mp);
		
		// chemical reactions
        double phiw = pm->Porosity(mp);
		for (i=0; i<nreact; ++i)
			phiwhat += phiw*pm->GetReaction(i)->m_Vbar*pm->GetReaction(i)->ReactionSupply(mp);

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z + H[i]*phiwhat)*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
/*
void FEMultiphasicDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;
	
	// element stiffness matrix
	matrix ke;
	
	vector<int> elm;
	
	// repeat over all solid elements
	int NE = m_Elem.size();
	
	if (fem.m_pStep->m_nanalysis == FE_STEADY_STATE) {
		for (int iel=0; iel<NE; ++iel)
		{
			FESolidElement& el = m_Elem[iel];
			
			UnpackLM(el, elm);
			
			// get the elements material
			FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (fem.GetMaterial(el.GetMatID()));
			assert(pm);
			const int nsol = pm->m_pSolute.size();
			
			// allocate stiffness matrix
			int neln = el.Nodes();
			int ndpn = 4+nsol;
			int ndof = neln*ndpn;
			ke.resize(ndof, ndof);
			
			// calculate the element stiffness matrix
			ElementMultiphasicStiffnessSS(fem, el, ke);
			
			// TODO: the problem here is that the LM array that is returned by the UnpackLM
			// function does not give the equation numbers in the right order. For this reason we
			// have to create a new lm array and place the equation numbers in the right order.
			// What we really ought to do is fix the UnpackLM function so that it returns
			// the LM vector in the right order for solute-solid elements.
			vector<int> lm(ndof);
			int isol;
			vector<int> dofc(nsol);
			for (isol=0; isol<nsol; ++isol)
				dofc[isol] = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
			for (int i=0; i<neln; ++i)
			{
				lm[ndpn*i  ] = elm[3*i];
				lm[ndpn*i+1] = elm[3*i+1];
				lm[ndpn*i+2] = elm[3*i+2];
				lm[ndpn*i+3] = elm[3*neln+i];
				for (isol=0; isol<nsol; ++isol)
					lm[ndpn*i+4+isol] = elm[dofc[isol]*neln+i];
			}
			
			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	} else {
		for (int iel=0; iel<NE; ++iel)
		{
			FESolidElement& el = m_Elem[iel];
			
			UnpackLM(el, elm);
			
			// get the elements material
			FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (fem.GetMaterial(el.GetMatID()));
			assert(pm);
			const int nsol = pm->m_pSolute.size();
			
			// allocate stiffness matrix
			int neln = el.Nodes();
			int ndpn = 4+nsol;
			int ndof = neln*ndpn;
			ke.resize(ndof, ndof);
			
			// calculate the element stiffness matrix
			ElementMultiphasicStiffness(fem, el, ke);
			
			// TODO: the problem here is that the LM array that is returned by the UnpackLM
			// function does not give the equation numbers in the right order. For this reason we
			// have to create a new lm array and place the equation numbers in the right order.
			// What we really ought to do is fix the UnpackLM function so that it returns
			// the LM vector in the right order for solute-solid elements.
			vector<int> lm(ndof);
			int isol;
			vector<int> dofc(nsol);
			for (isol=0; isol<nsol; ++isol)
				dofc[isol] = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
			for (int i=0; i<neln; ++i)
			{
				lm[ndpn*i  ] = elm[3*i];
				lm[ndpn*i+1] = elm[3*i+1];
				lm[ndpn*i+2] = elm[3*i+2];
				lm[ndpn*i+3] = elm[3*neln+i];
				for (isol=0; isol<nsol; ++isol)
					lm[ndpn*i+4+isol] = elm[dofc[isol]*neln+i];
			}
			
			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}
*/

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp)
{
	// element stiffness matrix
	matrix ke;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->Solutes();

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		UnpackLM(el, elm);
			
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 4+nsol;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementMultiphasicStiffness(el, ke, bsymm, tp.dt);
			
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for solute-solid elements.
		vector<int> lm(ndof);
		int isol;
		vector<int> dofc(nsol);
		for (isol=0; isol<nsol; ++isol)
			dofc[isol] = DOF_C + pm->GetSolute(isol)->GetSoluteID();
		for (int i=0; i<neln; ++i)
		{
			lm[ndpn*i  ] = elm[3*i];
			lm[ndpn*i+1] = elm[3*i+1];
			lm[ndpn*i+2] = elm[3*i+2];
			lm[ndpn*i+3] = elm[3*neln+i];
			for (isol=0; isol<nsol; ++isol)
				lm[ndpn*i+4+isol] = elm[dofc[isol]*neln+i];
		}
			
		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

void FEMultiphasicDomain::StiffnessMatrixSS(FESolver* psolver, bool bsymm, const FETimePoint& tp)
{
	// element stiffness matrix
	matrix ke;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->Solutes();

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		UnpackLM(el, elm);
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 4+nsol;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
			
		// calculate the element stiffness matrix
		ElementMultiphasicStiffnessSS(el, ke, bsymm, tp.dt);
			
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for solute-solid elements.
		vector<int> lm(ndof);
		int isol;
		vector<int> dofc(nsol);
		for (isol=0; isol<nsol; ++isol)
			dofc[isol] = DOF_C + pm->GetSolute(isol)->GetSoluteID();
		for (int i=0; i<neln; ++i)
		{
			lm[ndpn*i  ] = elm[3*i];
			lm[ndpn*i+1] = elm[3*i+1];
			lm[ndpn*i+2] = elm[3*i+2];
			lm[ndpn*i+3] = elm[3*neln+i];
			for (isol=0; isol<nsol; ++isol)
				lm[ndpn*i+4+isol] = elm[dofc[isol]*neln+i];
		}
			
		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEMultiphasicDomain::ElementMultiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt)
{
	int i, j, isol, jsol, ksol, n, ireact, isbm;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	// Gradient of shape functions
	vector<vec3d> gradN(neln);
	double tmp;
	
	// gauss-weights
	double* gw = el.GaussWeights();
	
	FEMesh& mesh = *m_pMesh;
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	
	const int nsol = pm->Solutes();
	int ndpn = 4+nsol;
	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->GetSolute(isol)->GetSoluteID();
	
	const int nsbm   = pm->SBMs();
	const int nreact = pm->Reactions();

	vec3d rp[FEElement::MAX_NODES], v[FEElement::MAX_NODES];
	vector< vector<double> > cp(nsol, vector<double>(FEElement::MAX_NODES));
	for (i=0; i<neln; ++i)
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		v[i]  = mesh.Node(el.m_node[i]).m_vt;
		for (isol=0; isol<nsol; ++isol)
			cp[isol][i] = mesh.Node(el.m_node[i]).m_cp[sid[isol]];
	}
	
	// zero stiffness matrix
	ke.zero();
	
	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(el, ks);
	
	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[ndpn*i  ][ndpn*j] = ks[3*i  ][3*j  ]; ke[ndpn*i  ][ndpn*j+1] = ks[3*i  ][3*j+1]; ke[ndpn*i  ][ndpn*j+2] = ks[3*i  ][3*j+2];
			ke[ndpn*i+1][ndpn*j] = ks[3*i+1][3*j  ]; ke[ndpn*i+1][ndpn*j+1] = ks[3*i+1][3*j+1]; ke[ndpn*i+1][ndpn*j+2] = ks[3*i+1][3*j+2];
			ke[ndpn*i+2][ndpn*j] = ks[3*i+2][3*j  ]; ke[ndpn*i+2][ndpn*j+1] = ks[3*i+2][3*j+1]; ke[ndpn*i+2][ndpn*j+2] = ks[3*i+2][3*j+2];
		}
	

	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
		vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
		vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		invjac0(el, J0i, n);
		vec3d G1(J0i[0][0],J0i[0][1],J0i[0][2]);
		vec3d G2(J0i[1][0],J0i[1][1],J0i[1][2]);
		vec3d G3(J0i[2][0],J0i[2][1],J0i[2][2]);
		
		// next we calculate the deformation gradient and the solid velocity
		mat3d Fp, gradv;
		Fp.zero();
		gradv.zero();
		vec3d vs(0);
		vector<double> cprev(nsol,0);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			GX = J0i[0][0]*Gr[i]+J0i[1][0]*Gs[i]+J0i[2][0]*Gt[i];
			GY = J0i[0][1]*Gr[i]+J0i[1][1]*Gs[i]+J0i[2][1]*Gt[i];
			GZ = J0i[0][2]*Gr[i]+J0i[1][2]*Gs[i]+J0i[2][2]*Gt[i];
			
			Fp[0][0] += rp[i].x*GX; Fp[1][0] += rp[i].y*GX; Fp[2][0] += rp[i].z*GX;
			Fp[0][1] += rp[i].x*GY; Fp[1][1] += rp[i].y*GY; Fp[2][1] += rp[i].z*GY;
			Fp[0][2] += rp[i].x*GZ; Fp[1][2] += rp[i].y*GZ; Fp[2][2] += rp[i].z*GZ;
			
			// calculate solid velocity and its gradient
			vs += v[i]*H[i];
			gradv[0][0] += v[i].x*Gx; gradv[1][0] += v[i].y*Gx; gradv[2][0] += v[i].z*Gx;
			gradv[0][1] += v[i].x*Gy; gradv[1][1] += v[i].y*Gy; gradv[2][1] += v[i].z*Gy;
			gradv[0][2] += v[i].x*Gz; gradv[1][2] += v[i].y*Gz; gradv[2][2] += v[i].z*Gz;
			
			// calculate Bp matrix
			gradN[i] = vec3d(Gx,Gy,Gz);
			
			// calculate effective concentration at previous time step
			for (isol=0; isol<nsol; ++isol)
				cprev[isol] += cp[isol][i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.m_J;
		
		// and then finally
		double dJdt = (J-Jp)/dt;
		double divv = dJdt/J;
		
		// get the fluid flux and pressure gradient
		vec3d w = ppt.m_w;
		vec3d gradp = ppt.m_gradp;
		
		vector<double> c(spt.m_c);
		vector<vec3d> gradc(spt.m_gradc);
		vector<double> dcdt(nsol);
		vector<int> z(nsol);
		
		vector<double> zz(nsol);
		vector<double> kappa(spt.m_k);
		
		for (isol=0; isol<nsol; ++isol) {
			// get the time derivative of the effective concentration
			dcdt[isol] = (c[isol] - cprev[isol])/dt;
			// get the charge number
			z[isol] = pm->GetSolute(isol)->ChargeNumber();
		}
		
		vector<double> dkdJ(spt.m_dkdJ);
		vector< vector<double> > dkdc(spt.m_dkdc);
		vector<double> dkdJJ(spt.m_dkdJJ);
		vector< vector<double> > dkdJc(spt.m_dkdJc);
		vector< vector< vector<double> > > dkdcc(spt.m_dkdcc);
		vector< vector<double> > dkdr(spt.m_dkdr);
		vector< vector<double> > dkdJr(spt.m_dkdJr);
		vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
        double phi0 = ppt.m_phi0;
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		double dpdJJ = -2*phis/(J*J);
		
		// evaluate the osmotic coefficient
		double osmc = pm->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
		
		// evaluate the permeability
		mat3ds K = pm->GetPermeability()->Permeability(mp);
		tens4ds dKdE = pm->GetPermeability()->Tangent_Permeability_Strain(mp);
		
		vector<mat3ds> dKdc(nsol);
		vector<mat3ds> D(nsol);
		vector<tens4ds> dDdE(nsol);
		vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
		vector<double> D0(nsol);
		vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
		vector<double> dodc(nsol);
		vector<mat3ds> dTdc(nsol);
		vector<mat3ds> ImD(nsol);
		mat3dd I(1);
		
		// evaluate the solvent supply and its derivatives
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		vector<double> Phic(nsol,0);
        vector<mat3ds> dchatde(nsol);
		if (pm->GetSolventSupply()) {
			Phie = pm->GetSolventSupply()->Tangent_Supply_Strain(mp);
			Phip = pm->GetSolventSupply()->Tangent_Supply_Pressure(mp);
		}
        
        // chemical reactions
		for (i=0; i<nreact; ++i)
			Phie += pm->GetReaction(i)->m_Vbar*(I*pm->GetReaction(i)->ReactionSupply(mp)
              +pm->GetReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
		
		for (isol=0; isol<nsol; ++isol) {
			// evaluate the permeability derivatives
			dKdc[isol] = pm->GetPermeability()->Tangent_Permeability_Concentration(mp,isol);
			
			// evaluate the diffusivity tensor and its derivatives
			D[isol] = pm->GetSolute(isol)->m_pDiff->Diffusivity(mp);
			dDdE[isol] = pm->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Strain(mp);
			
			// evaluate the solute free diffusivity
			D0[isol] = pm->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
			
			// evaluate the derivative of the osmotic coefficient
			dodc[isol] = pm->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
			
			// evaluate the stress tangent with concentration
			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
			}
			
			// evaluate the solvent supply tangent with concentration
			if (pm->GetSolventSupply()) Phic[isol] = pm->GetSolventSupply()->Tangent_Supply_Concentration(mp,isol);

            // chemical reactions
            dchatde[isol].zero();
            for (ireact=0; ireact<nreact; ++ireact) {
                dchatde[isol] += pm->GetReaction(ireact)->m_v[isol]
                *(I*pm->GetReaction(ireact)->ReactionSupply(mp)
                  +pm->GetReaction(ireact)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
                Phic[isol] += phiw*pm->GetReaction(ireact)->m_Vbar
                *pm->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, isol);
            }
		}
		
		// Miscellaneous constants
		double R = pm->m_Rgas;
		double T = pm->m_Tabs;
		double penalty = pm->m_penalty;
		
		// evaluate the effective permeability and its derivatives
		mat3ds Ki = K.inverse();
		mat3ds Ke(0,0,0,0,0,0);
		tens4ds G = dyad1s(Ki,I) - dyad4s(Ki,I)*2 - ddots(dyad2s(Ki),dKdE)*0.5;
		vector<mat3ds> Gc(nsol);
		vector<mat3ds> dKedc(nsol);
		for (isol=0; isol<nsol; ++isol) {
			Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
			G += dyad1s(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/2/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
			+(dyad1s(I) - dyad4s(I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
			Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
			for (jsol=0; jsol<nsol; ++jsol) {
				Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
				-(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
			}
			Gc[isol] *= R*T/phiw;
		}
		Ke = (Ki + Ke*(R*T/phiw)).inverse();
		tens4ds dKedE = dyad1s(Ke,I) - 2*dyad4s(Ke,I) - ddots(dyad2s(Ke),G)*0.5;
		for (isol=0; isol<nsol; ++isol)
			dKedc[isol] = -Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke;
		
		// calculate all the matrices
		vec3d vtmp,gp,qpu;
		vector<vec3d> gc(nsol),qcu(nsol),wc(nsol),jce(nsol);
		vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
		mat3d wu, jue;
		vector<mat3d> ju(nsol);
		vector< vector<double> > qcc(nsol, vector<double>(nsol));
		vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
		double sum;
		mat3ds De;
		tmp = detJ*gw[n];
		for (i=0; i<neln; ++i)
		{
			for (j=0; j<neln; ++j)
			{
				// calculate the kpu matrix
				gp = vec3d(0,0,0);
				for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
				gp = gradp+gp*(R*T);
				wu = vdotTdotv(-gp, dKedE, gradN[j]);
				for (isol=0; isol<nsol; ++isol) {
					wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
							+Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
					+ (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
				}
//				qpu = -gradN[j]*(divv+1.0/dt)-gradv.transpose()*gradN[j];
				qpu = -gradN[j]*(divv+1.0/dt)+gradv.transpose()*gradN[j];
				vtmp = (wu.transpose()*gradN[i] + (qpu + Phie*gradN[j])*H[i])*(tmp*dt);
				ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
				
				// calculate the kup matrix
				vtmp = -gradN[i]*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
				
				// calculate the kpp matrix
				ke[ndpn*i+3][ndpn*j+3] += (H[i]*H[j]*Phip - gradN[i]*(Ke*gradN[j]))*(tmp*dt);
				
				// calculate kcu matrix data
				jue.zero();
				De.zero();
				for (isol=0; isol<nsol; ++isol) {
					gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
					ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol]) 
					+ vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
					+ (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
					   +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
					   )*kappa[isol]
					+D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
					jue += ju[isol]*z[isol];
					De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
					sum = 0;
					for (jsol=0; jsol<nsol; ++jsol)
						sum += dcdt[jsol]*((phiw+J*dpdJ)*dkdc[isol][jsol]+J*phiw*dkdJc[isol][jsol]);
					qcu[isol] = -gradN[j]*(c[isol]*(dJdt*(2*(dpdJ*kappa[isol]+phiw*dkdJ[isol]+J*dpdJ*dkdJ[isol])
														  +J*(dpdJJ*kappa[isol]+phiw*dkdJJ[isol])) + sum)
										   +dcdt[isol]*((phiw+J*dpdJ)*kappa[isol]+J*phiw*dkdJ[isol]))
					+qpu*(c[isol]*((phiw+J*dpdJ)*kappa[isol]+J*phiw*dkdJ[isol]));

                    // chemical reactions
                    for (ireact=0; ireact<nreact; ++ireact) {
                        double sum1 = 0;
                        double sum2 = 0;
                        for (isbm=0; isbm<nsbm; ++isbm) {
                            sum1 += pm->SBMMolarMass(isbm)*pm->GetReaction(ireact)->m_v[nsol+isbm]*
                            ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/pm->SBMDensity(isbm));
                            sum2 += pm->SBMMolarMass(isbm)*pm->GetReaction(ireact)->m_v[nsol+isbm]*
                            (dkdr[isol][isbm]+(J-phi0)*dkdJr[isol][isbm]-dkdJ[isol]/pm->SBMDensity(isbm));
                        }
                        double zhat = pm->GetReaction(ireact)->ReactionSupply(mp);
                        mat3dd zhatI(zhat);
                        mat3ds dzde = pm->GetReaction(ireact)->Tangent_ReactionSupply_Strain(mp);
                        qcu[isol] -= ((zhatI+dzde*(J-phi0))*gradN[j])*(sum1*c[isol])
                        +gradN[j]*(c[isol]*(J-phi0)*sum2*zhat);
                    }
				}
				
				for (isol=0; isol<nsol; ++isol) {
					
					// calculate the kcu matrix
					vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i]
							+ (qcu[isol] + dchatde[isol]*gradN[j])*H[i])*(tmp*dt);
					ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
					ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
					ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
					
					// calculate the kcp matrix
					ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
															  (D[isol]*(kappa[isol]*c[isol]/D0[isol])
															   +De*penalty)
															  *(Ke*gradN[j])
															  ))*(tmp*dt);
				
					// calculate the kuc matrix
					sum = 0;
					for (jsol=0; jsol<nsol; ++jsol)
						sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
					vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*tmp;
					ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
					ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
					ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
					
					// calculate the kpc matrix
					vtmp = vec3d(0,0,0);
					for (jsol=0; jsol<nsol; ++jsol)
						vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
								+dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
					wc[isol] = (dKedc[isol]*gp)*(-H[j])
					-Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
					ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol]+H[i]*H[j]*Phic[isol])*(tmp*dt);
					
				}
				
				// calculate data for the kcc matrix
				jce.assign(nsol, vec3d(0,0,0));
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						if (jsol != isol) {
							jc[isol][jsol] = 
							((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
							+(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
							
							sum = dkdc[isol][jsol]*dcdt[isol];
							for (ksol=0; ksol<nsol; ++ksol)
								sum += c[isol]*dkdcc[isol][jsol][ksol]*dcdt[ksol];
							
							qcc[isol][jsol] = -H[j]*((c[isol]*((phiw+J*dpdJ)*dkdc[isol][jsol]+J*phiw*dkdJc[isol][jsol]))*divv
											   +phiw*(c[isol]*dkdc[isol][jsol])/dt
											   +phiw*sum);
						}
						else {
							jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
							+((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
							+(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
							
							sum = dkdc[isol][jsol]*dcdt[isol];
							for (ksol=0; ksol<nsol; ++ksol)
								sum += (dkdc[isol][ksol]+c[isol]*dkdcc[isol][jsol][ksol])*dcdt[ksol];
							
							qcc[isol][jsol] = -H[j]*(((phiw+J*dpdJ)*kappa[isol]+J*phiw*dkdJ[isol]
													  +c[isol]*((phiw+J*dpdJ)*dkdc[isol][jsol]+J*phiw*dkdJc[isol][jsol]))*divv
													 +phiw*(kappa[isol] + c[isol]*dkdc[isol][jsol])/dt
													 +phiw*sum);
						}
						jce[jsol] += jc[isol][jsol]*z[isol];
						
						// chemical reactions
						dchatdc[isol][jsol] = 0;
						for (ireact=0; ireact<nreact; ++ireact) {
							dchatdc[isol][jsol] += pm->GetReaction(ireact)->m_v[isol]
							*pm->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                            double sum1 = 0;
                            double sum2 = 0;
                            for (isbm=0; isbm<nsbm; ++isbm) {
                                sum1 += pm->SBMMolarMass(isbm)*pm->GetReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/pm->SBMDensity(isbm));
                                sum2 += pm->SBMMolarMass(isbm)*pm->GetReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdrc[isol][isbm][jsol]-dkdc[isol][jsol]/pm->SBMDensity(isbm));
                            }
                            double zhat = pm->GetReaction(ireact)->ReactionSupply(mp);
                            double dzdc = pm->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, jsol);
                            if (jsol != isol) {
                                qcc[isol][jsol] -= H[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                            }
                            else {
                                qcc[isol][jsol] -= H[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                            }
                        }
					}
				}

				// calculate the kcc matrix
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
															 + H[i]*(qcc[isol][jsol]
																	 + H[j]*phiw*dchatdc[isol][jsol]))*(tmp*dt);
					}
				}
			}
		}
	}
	
	// Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
	if (bsymm) {
		for (i=0; i<ndpn*neln; ++i)
			for (j=i+1; j<ndpn*neln; ++j) {
				tmp = 0.5*(ke[i][j]+ke[j][i]);
				ke[i][j] = ke[j][i] = tmp;
			}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//!
bool FEMultiphasicDomain::ElementMultiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt)
{
	int i, j, isol, jsol, n, ireact;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;
	
	// jacobian
	double Ji[3][3], detJ;
	
	// Gradient of shape functions
	vector<vec3d> gradN(neln);
	double tmp;
	
	// gauss-weights
	double* gw = el.GaussWeights();
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	
	const int nsol = pm->Solutes();
	int ndpn = 4+nsol;
	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->GetSolute(isol)->GetSoluteID();
	
	const int nreact = pm->Reactions();

	// zero stiffness matrix
	ke.zero();
	
	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(el, ks);
	
	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[ndpn*i  ][ndpn*j] = ks[3*i  ][3*j  ]; ke[ndpn*i  ][ndpn*j+1] = ks[3*i  ][3*j+1]; ke[ndpn*i  ][ndpn*j+2] = ks[3*i  ][3*j+2];
			ke[ndpn*i+1][ndpn*j] = ks[3*i+1][3*j  ]; ke[ndpn*i+1][ndpn*j+1] = ks[3*i+1][3*j+1]; ke[ndpn*i+1][ndpn*j+2] = ks[3*i+1][3*j+2];
			ke[ndpn*i+2][ndpn*j] = ks[3*i+2][3*j  ]; ke[ndpn*i+2][ndpn*j+1] = ks[3*i+2][3*j+1]; ke[ndpn*i+2][ndpn*j+2] = ks[3*i+2][3*j+2];
		}
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
		vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
		vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			// calculate Bp matrix
			gradN[i] = vec3d(Gx,Gy,Gz);
			
		}
		
		// next we get the determinant
		double J = ept.m_J;
		
		// get the fluid flux and pressure gradient
		vec3d w = ppt.m_w;
		vec3d gradp = ppt.m_gradp;
		
		vector<double> c(spt.m_c);
		vector<vec3d> gradc(spt.m_gradc);
		vector<int> z(nsol);
		
		vector<double> zz(nsol);
		vector<double> kappa(spt.m_k);
		
		for (isol=0; isol<nsol; ++isol) {
			// get the charge number
			z[isol] = pm->GetSolute(isol)->ChargeNumber();
		}
		
		vector<double> dkdJ(spt.m_dkdJ);
		vector< vector<double> > dkdc(spt.m_dkdc);
		vector<double> dkdJJ(spt.m_dkdJJ);
		vector< vector<double> > dkdJc(spt.m_dkdJc);
		vector< vector< vector<double> > > dkdcc(spt.m_dkdcc);
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		
		// evaluate the osmotic coefficient
		double osmc = pm->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
		
		// evaluate the permeability
		mat3ds K = pm->GetPermeability()->Permeability(mp);
		tens4ds dKdE = pm->GetPermeability()->Tangent_Permeability_Strain(mp);
		
		vector<mat3ds> dKdc(nsol);
		vector<mat3ds> D(nsol);
		vector<tens4ds> dDdE(nsol);
		vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
		vector<double> D0(nsol);
		vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
		vector<double> dodc(nsol);
		vector<mat3ds> dTdc(nsol);
		vector<mat3ds> ImD(nsol);
		mat3dd I(1);
		
		// evaluate the solvent supply and its derivatives
		double phiwhat = 0;
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		vector<double> Phic(nsol,0);
		if (pm->GetSolventSupply()) {
			phiwhat = pm->GetSolventSupply()->Supply(mp);
			Phie = pm->GetSolventSupply()->Tangent_Supply_Strain(mp);
			Phip = pm->GetSolventSupply()->Tangent_Supply_Pressure(mp);
		}
		
        // chemical reactions
		for (i=0; i<nreact; ++i)
			Phie += pm->GetReaction(i)->m_Vbar*(I*pm->GetReaction(i)->ReactionSupply(mp)
                                             +pm->GetReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
		
		for (isol=0; isol<nsol; ++isol) {
			// evaluate the permeability derivatives
			dKdc[isol] = pm->GetPermeability()->Tangent_Permeability_Concentration(mp,isol);
			
			// evaluate the diffusivity tensor and its derivatives
			D[isol] = pm->GetSolute(isol)->m_pDiff->Diffusivity(mp);
			dDdE[isol] = pm->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Strain(mp);
			
			// evaluate the solute free diffusivity
			D0[isol] = pm->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
			
			// evaluate the derivative of the osmotic coefficient
			dodc[isol] = pm->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
			
			// evaluate the stress tangent with concentration
			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
			}

			// evaluate the solvent supply tangent with concentration
			if (pm->GetSolventSupply()) Phic[isol] = pm->GetSolventSupply()->Tangent_Supply_Concentration(mp,isol);
			
		}
		
		// Miscellaneous constants
		double R = pm->m_Rgas;
		double T = pm->m_Tabs;
		double penalty = pm->m_penalty;
		
		// evaluate the effective permeability and its derivatives
		mat3ds Ki = K.inverse();
		mat3ds Ke(0,0,0,0,0,0);
		tens4ds G = dyad1s(Ki,I) - dyad4s(Ki,I)*2 - ddots(dyad2s(Ki),dKdE)*0.5;
		vector<mat3ds> Gc(nsol);
		vector<mat3ds> dKedc(nsol);
		for (isol=0; isol<nsol; ++isol) {
			Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
			G += dyad1s(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/2/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
			+(dyad1s(I) - dyad4s(I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
			Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
			for (jsol=0; jsol<nsol; ++jsol) {
				Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
				-(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
			}
			Gc[isol] *= R*T/phiw;
		}
		Ke = (Ki + Ke*(R*T/phiw)).inverse();
		tens4ds dKedE = dyad1s(Ke,I) - 2*dyad4s(Ke,I) - ddots(dyad2s(Ke),G)*0.5;
		for (isol=0; isol<nsol; ++isol)
			dKedc[isol] = -Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke;
		
		// calculate all the matrices
		vec3d vtmp,gp,qpu;
		vector<vec3d> gc(nsol),wc(nsol),jce(nsol);
		vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
		mat3d wu, jue;
		vector<mat3d> ju(nsol);
		vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
		double sum;
		mat3ds De;
		tmp = detJ*gw[n];
		for (i=0; i<neln; ++i)
		{
			for (j=0; j<neln; ++j)
			{
				// calculate the kpu matrix
				gp = vec3d(0,0,0);
				for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
				gp = gradp+gp*(R*T);
				wu = vdotTdotv(-gp, dKedE, gradN[j]);
				for (isol=0; isol<nsol; ++isol) {
					wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
						   +Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
					+ (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
				}
				qpu = Phie*gradN[j];
				vtmp = (wu.transpose()*gradN[i] + qpu*H[i])*(tmp*dt);
				ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
				
				// calculate the kup matrix
				vtmp = -gradN[i]*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
				
				// calculate the kpp matrix
				ke[ndpn*i+3][ndpn*j+3] += (H[i]*H[j]*Phip - gradN[i]*(Ke*gradN[j]))*(tmp*dt);
				
				// calculate kcu matrix data
				jue.zero();
				De.zero();
				for (isol=0; isol<nsol; ++isol) {
					gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
					ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol]) 
					+ vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
					+ (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
					   +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
					   )*kappa[isol]
					+D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
					jue += ju[isol]*z[isol];
					De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
				}
				
				for (isol=0; isol<nsol; ++isol) {
					
					// calculate the kcu matrix
					vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i])*(tmp*dt);
					ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
					ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
					ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
					
					// calculate the kcp matrix
					ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
															  (D[isol]*(kappa[isol]*c[isol]/D0[isol])
															   +De*penalty)
															  *(Ke*gradN[j])
															  ))*(tmp*dt);
					
					// calculate the kuc matrix
					sum = 0;
					for (jsol=0; jsol<nsol; ++jsol)
						sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
					vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*tmp;
					ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
					ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
					ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
					
					// calculate the kpc matrix
					vtmp = vec3d(0,0,0);
					for (jsol=0; jsol<nsol; ++jsol)
						vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
								 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
					wc[isol] = (dKedc[isol]*gp)*(-H[j])
					-Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
					ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol]+H[i]*H[j]*Phic[isol])*(tmp*dt);
					
				}
				
				// calculate data for the kcc matrix
				jce.assign(nsol, vec3d(0,0,0));
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						if (jsol != isol) {
							jc[isol][jsol] = 
							((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
							+(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
						}
						else {
							jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
							+((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
							+(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
						}
						jce[jsol] += jc[isol][jsol]*z[isol];
						
						// chemical reactions
						dchatdc[isol][jsol] = 0;
						for (ireact=0; ireact<nreact; ++ireact)
							dchatdc[isol][jsol] += pm->GetReaction(ireact)->m_v[isol]
							*pm->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
					}
				}
				
				// calculate the kcc matrix
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
															 + H[i]*H[j]*phiw*dchatdc[isol][jsol])*(tmp*dt);
					}
				}
			}
		}
	}
	
	// Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
	if (bsymm) {
		for (i=0; i<ndpn*neln; ++i)
			for (j=i+1; j<ndpn*neln; ++j) {
				tmp = 0.5*(ke[i][j]+ke[j][i]);
				ke[i][j] = ke[j][i] = tmp;
			}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FEMultiphasicDomain::SolidElementStiffness(FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	ElementMultiphasicMaterialStiffness(el, ke);
	
	// calculate geometrical stiffness
	ElementGeometricalStiffness(el, ke);
	
	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	int ndof = 3*el.Nodes();
	int i, j;
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEMultiphasicDomain::ElementMultiphasicMaterialStiffness(FESolidElement &el, matrix &ke)
{
	// see if this is a multiphasic material
	FEMultiphasic* pmat = dynamic_cast<FEMultiphasic*>(GetMaterial());
	assert(pmat);
	
	int i, i3, j, j3, n;
	
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	
	// global derivatives of shape functions
	const int NE = FEElement::MAX_NODES;
	double Gx[NE], Gy[NE], Gz[NE];
	
	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;
	
	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix
	
	// The 'D*BL' matrix
	double DBL[6][3];
	
	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;
	
	// jacobian
	double Ji[3][3], detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();
	
	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n)*gw[n];
		
		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);
		
		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.m_State[n];
		
		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
		C.extract(D);
		
		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];
			
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}
		
		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = Gx[i];
			Gyi = Gy[i];
			Gzi = Gz[i];
			
			for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = Gx[j];
				Gyj = Gy[j];
				Gzj = Gz[j];
				
				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);
				
				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);
				
				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);
				
				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);
				
				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);
				
				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);
				
				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;
				
				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;
				
				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::UpdateStresses(FEModel &fem)
{
	bool berr = false;
	int NE = (int) m_Elem.size();
	double dt = fem.GetCurrentStep()->m_dt;
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			UpdateElementStress(i, dt);
		}
		catch (NegativeJacobian e)
		{
			// A negative jacobian was detected
			felog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
			#pragma omp critical
			berr = true;
		}
	}

	// if we encountered an error, we request a running restart
	if (berr) throw DoRunningRestart();
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::UpdateElementStress(int iel, double dt)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];

   // extract the elastic component
    FEElasticMaterial* pme = m_pMat->GetElasticMaterial();

	// get the multiphasic material
	FEMultiphasic* pmb = dynamic_cast<FEMultiphasic*>(GetMaterial());
	assert(pmb);
	const int nsol = pmb->Solutes();
	const int nsbm = pmb->SBMs();
	vector<int> sid(nsol);
	for (int j=0; j<nsol; ++j) sid[j] = pmb->GetSolute(j)->GetSoluteID();
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	// get the number of nodes
	int neln = el.Nodes();
		
	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES];
	vector< vector<double> > ct(nsol, vector<double>(FEElement::MAX_NODES));
	for (int j=0; j<neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
		pn[j] = mesh.Node(el.m_node[j]).m_pt;
		for (int k=0; k<nsol; ++k)
			ct[k][j] = mesh.Node(el.m_node[j]).m_ct[sid[k]];
	}

	// check if this mixture includes chemical reactions
	int nreact = pmb->Reactions();
	if (nreact) 
	{
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEElasticMaterialPoint&  pt  = *(mp.ExtractData<FEElasticMaterialPoint >());
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());

			// for chemical reactions involving solid-bound molecules,
			// update their concentration
            double phi0 = ppt.m_phi0;
			for (int isbm=0; isbm<nsbm; ++isbm) {
				spt.m_sbmrhat[isbm] = 0;
				// combine the molar supplies from all the reactions
				for (int k=0; k<nreact; ++k) {
					double zetahat = pmb->GetReaction(k)->ReactionSupply(mp);
					double v = pmb->GetReaction(k)->m_v[nsol+isbm];
					// remember to convert from molar supply to referential mass supply
					spt.m_sbmrhat[isbm] += (pt.m_J-phi0)*pmb->SBMMolarMass(isbm)*v*zetahat;
				}
				// perform the time integration (Euler's method)
				spt.m_sbmr[isbm] = spt.m_sbmrp[isbm] + dt*spt.m_sbmrhat[isbm];
				// check bounds
				if (spt.m_sbmr[isbm] < pmb->GetSBM(isbm)->m_rhomin)
					spt.m_sbmr[isbm] = pmb->GetSBM(isbm)->m_rhomin;
				if ((pmb->GetSBM(isbm)->m_rhomax > 0) && (spt.m_sbmr[isbm] > pmb->GetSBM(isbm)->m_rhomax))
					spt.m_sbmr[isbm] = pmb->GetSBM(isbm)->m_rhomax;
			}
		}
	}
		
	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			
		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);
			
		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
			
		// solute-poroelastic data
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());

        // evaluate referential solid volume fraction
        ppt.m_phi0 = pmb->SolidReferentialVolumeFraction(mp);

		// evaluate fluid pressure at gauss-point
		ppt.m_p = el.Evaluate(pn, n);
			
		// calculate the gradient of p at gauss-point
		ppt.m_gradp = gradient(el, pn, n);
			
		for (int k=0; k<nsol; ++k) {
			// evaluate effective solute concentrations at gauss-point
			spt.m_c[k] = el.Evaluate(&ct[k][0], n);
			// calculate the gradient of c at gauss-point
			spt.m_gradc[k] = gradient(el, &ct[k][0], n);
		}
			
		// update the fluid and solute fluxes
		// and evaluate the actual fluid pressure and solute concentration
		ppt.m_w = pmb->FluidFlux(mp);
		spt.m_psi = pmb->ElectricPotential(mp);
		for (int k=0; k<nsol; ++k) {
			spt.m_ca[k] = pmb->Concentration(mp,k);
			spt.m_j[k] = pmb->SoluteFlux(mp,k);
		}
		ppt.m_pa = pmb->Pressure(mp);
		spt.m_cF = pmb->FixedChargeDensity(mp);
		spt.m_Ie = pmb->CurrentDensity(mp);
           pmb->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc,
                                              spt.m_dkdJJ, spt.m_dkdJc, spt.m_dkdcc,
                                              spt.m_dkdr, spt.m_dkdJr, spt.m_dkdrc);
			
		pt.m_s = pmb->Stress(mp);

		// evaluate the strain energy density
		pt.m_sed = pme->StrainEnergy(mp);
            
		// evaluate the referential solid density
		pt.m_rhor = pmb->SolidReferentialApparentDensity(mp);
	}
}
