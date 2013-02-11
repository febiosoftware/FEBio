#include "stdafx.h"
#include "FEMultiphasicDomain.h"
#include "FECore/FEMaterial.h"
#include "FEMultiphasic.h"
#include "FECore/log.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
bool FEMultiphasicDomain::Initialize(FEModel &mdl)
{
	// initialize base class
	FEElasticSolidDomain::Initialize(mdl);

	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());
		
	// get the multiphasic material
	FEMultiphasic* pmb = dynamic_cast<FEMultiphasic*>(pm); assert(pmb);
	const int nsol = pmb->m_pSolute.size();

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
			ps.m_j.assign(nsol,0);
		}
	}
	
	return true;
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
			
			// reset referential solid volume fraction at previous time
			pt.m_phi0p = pt.m_phi0;
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
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
			//! this element should not be UDG
			assert(el.Type() != FE_HEX8G1);
			
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
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
			//! this element should not be UDG
			assert(el.Type() != FE_HEX8G1);
			
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
void FEMultiphasicDomain::InternalSoluteWorkSS(FENLSolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	const int nsol = pm->m_pSolute.size();

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
			int dofc = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
			for (int j=0; j<neln; ++j)
			{
				int J = elm[dofc*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEMultiphasicDomain::ElementInternalSoluteWorkSS(FESolidElement& elem, vector<double>& fe, double dt, const int sol)
{
	// jacobian
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;
	
	// gradient of shape functions
	int neln = elem.Nodes();
	vector<vec3d> gradN(neln);
	
	// gauss-weights
	double* wg = elem.GaussWeights();
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->m_pSolute.size();
	vector<vec3d> j(nsol);
	vector<int> z(nsol);
	
	zero(fe);
	
	// loop over gauss-points
	int nint = elem.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *elem.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(elem, Ji, n);
		
		Gr = elem.Gr(n);
		Gs = elem.Gs(n);
		Gt = elem.Gt(n);
		
		H = elem.H(n);
		
		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];
			
			// save spatial gradient of shape functions
			gradN[i] = vec3d(Gx,Gy,Gz);
		}
		
		vec3d je(0);
		
		for (int isol=0; isol<nsol; ++isol) {
			// get the solute flux
			j[isol] = spt.m_j[isol];
			// get the charge number
			z[isol] = pm->m_pSolute[isol]->ChargeNumber();
			// current density (flux units)
			je += j[isol]*z[isol];
		}
		
		// update force vector
		for (int i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[sol]+je*pm->m_penalty)
						 )*detJ*wg[n];
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalSoluteWork(FENLSolver* psolver, vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial());
	assert(pm);
	const int nsol = pm->m_pSolute.size();

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
			int dofc = DOF_C + pm->m_pSolute[isol]->GetSoluteID();
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
	int i, isol, jsol;
	
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
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);

	const int nsol = pm->m_pSolute.size();
	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->m_pSolute[isol]->GetSoluteID();
	
	const int NE = FEElement::MAX_NODES;
	vec3d r0[NE], rt[NE], rp[NE], vt[NE];
	vector< vector<double> > cp(nsol, vector<double>(FEElement::MAX_NODES));
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
	for (int n=0; n<nint; ++n)
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
		double J = ept.J;
		double dJdt = (J-Jp)/dt;
		
		// and then finally
		double divv = dJdt/J;
		
		vector<vec3d> j(nsol);
		vector<double> c(nsol);
		vector<int> z(nsol);
		vector<double> khat(nsol);
		vector<double> dkhdJ(nsol);
		vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
		vector<double> zz(nsol);
		vector<double> kappa(nsol);
		double den = 0;
		vec3d je(0);
		vector<double> dcdt(nsol);
		double zeta = pm->ElectricPotential(mp, true);
		
		for (isol=0; isol<nsol; ++isol) {
			// get the solute flux
			j[isol] = spt.m_j[isol];
			// get the effective concentration and its time derivative
			c[isol] = spt.m_c[isol];
			dcdt[isol] = (c[isol] - cprev[isol])/dt;
			// get the charge number
			z[isol] = pm->m_pSolute[isol]->ChargeNumber();
			// evaluate the solubility and its derivatives w.r.t. J and c
			khat[isol] = pm->m_pSolute[isol]->m_pSolub->Solubility(mp);
			dkhdJ[isol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
			for (jsol=0; jsol<nsol; ++jsol) {
				dkhdc[isol][jsol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
			}
			// evaluate electric potential (nondimensional exponential form)
			// also evaluate partition coefficients
			zz[isol] = pow(zeta, z[isol]);
			kappa[isol] = zz[isol]*khat[isol];
			den += SQR(z[isol])*kappa[isol]*c[isol];
			je += j[isol]*z[isol];
		}
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		
		// evaluate derivatives of electric potential (nondimensional exponential form)
		// and partition coefficients
		double zidzdJ = 0;
		vector<double> zidzdc(nsol,0);
		if (den > 0) {
			for (isol=0; isol<nsol; ++isol) {
				zidzdJ += z[isol]*zz[isol]*dkhdJ[isol]*c[isol];
				for (jsol=0; jsol<nsol; ++jsol)
					zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
				zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
			}
			zidzdJ = -(dcFdJ+zidzdJ)/den;
		}
		double dkdJ;
		vector<double> dkdc(nsol);
		dkdJ = zz[sol]*dkhdJ[sol]+z[sol]*kappa[sol]*zidzdJ;
		for (jsol=0; jsol<nsol; ++jsol)
			dkdc[jsol] = zz[sol]*dkhdc[sol][jsol]+z[sol]*kappa[sol]*zidzdc[jsol];
		
		// evaluate the porosity, its derivative w.r.t. J, and its gradient
		double phiw = pm->Porosity(mp);
		double dpdJ = (1. - phiw)/J;
		// evaluate time derivatives of solubility and porosity
		double dkdt = dkdJ*dJdt;
		for (jsol=0; jsol<nsol; ++jsol)
			dkdt += dkdc[jsol]*dcdt[jsol];
		double dpdt = dpdJ*dJdt;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[sol]+je*pm->m_penalty)
						 - H[i]*(dpdt*kappa[sol]*c[sol]+phiw*dkdt*c[sol]
								 +phiw*kappa[sol]*dcdt[sol]+phiw*kappa[sol]*c[sol]*divv)
						 )*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalFluidWork(FENLSolver* psolver, vector<double>& R, double dt)
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
	int i;
	
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
	
	zero(fe);
	
	// loop over gauss-points
	for (int n=0; n<nint; ++n)
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
		double J = ept.J;
		
		// and then finally
		double divv = ((J-Jp)/dt)/J;
		
		// get the flux
		vec3d& w = ppt.m_w;

		// get the solvent supply
		double phiwhat = 0;
		if (pm->m_pSupp) phiwhat = pm->m_pSupp->Supply(mp);
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z + (phiwhat - divv)*H[i])*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicDomain::InternalFluidWorkSS(FENLSolver* psolver, vector<double>& R, double dt)
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
	int i;
	
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
	
	zero(fe);
	
	// loop over gauss-points
	for (int n=0; n<nint; ++n)
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
		
		// get the flux
		vec3d& w = ppt.m_w;

		// get the solvent supply
		double phiwhat = 0;
		if (pm->m_pSupp) phiwhat = pm->m_pSupp->Supply(mp);
		
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
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
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
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
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
void FEMultiphasicDomain::StiffnessMatrix(FENLSolver* psolver, bool bsymm, double dt)
{
	// element stiffness matrix
	matrix ke;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->m_pSolute.size();

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
		ElementMultiphasicStiffness(el, ke, bsymm, dt);
			
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

void FEMultiphasicDomain::StiffnessMatrixSS(FENLSolver* psolver, bool bsymm, double dt)
{
	// element stiffness matrix
	matrix ke;
	vector<int> elm;

	// get the elements material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->m_pSolute.size();

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
		ElementMultiphasicStiffnessSS(el, ke, bsymm, dt);
			
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

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEMultiphasicDomain::ElementMultiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt)
{
	int i, j, isol, jsol, ksol, n;
	
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
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->m_pSolute.size();
	int ndpn = 4+nsol;
	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->m_pSolute[isol]->GetSoluteID();
	
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
		double J = ept.J;
		
		// and then finally
		double dJdt = (J-Jp)/dt;
		double divv = dJdt/J;
		
		// get the fluid flux and pressure gradient
		vec3d w = ppt.m_w;
		vec3d gradp = ppt.m_gradp;
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		double dcFdJJ = 2*cF/SQR(J-phi0);
		
		vector<double> c(nsol);
		vector<vec3d> gradc(nsol);
		vector<double> dcdt(nsol);
		vector<int> z(nsol);
		vector<double> khat(nsol);
		vector<double> dkhdJ(nsol);
		vector<double> dkhdJJ(nsol);
		vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
		vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
		vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);	// use dkhdc to initialize only
		
		double zeta = pm->ElectricPotential(mp, true);
		vector<double> zz(nsol);
		vector<double> kappa(nsol);
		double den = 0;
		
		for (isol=0; isol<nsol; ++isol) {
			// get the effective concentration, its gradient and its time derivative
			c[isol] = spt.m_c[isol];
			gradc[isol] = spt.m_gradc[isol];
			dcdt[isol] = (c[isol] - cprev[isol])/dt;
			
			// get the charge number
			z[isol] = pm->m_pSolute[isol]->ChargeNumber();
			
			// evaluate the solubility and its derivatives w.r.t. J and c
			khat[isol] = pm->m_pSolute[isol]->m_pSolub->Solubility(mp);
			dkhdJ[isol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
			dkhdJJ[isol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Strain(mp);
			for (jsol=0; jsol<nsol; ++jsol) {
				dkhdc[isol][jsol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
				dkhdJc[isol][jsol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,jsol);
				for (ksol=0; ksol<nsol; ++ksol) {
					dkhdcc[isol][jsol][ksol] = 
					pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,jsol,ksol);
				}
			}
			zz[isol] = pow(zeta, z[isol]);
			kappa[isol] = zz[isol]*khat[isol];
			den += SQR(z[isol])*kappa[isol]*c[isol];
		}
		
		// evaluate electric potential (nondimensional exponential form) and its derivatives
		// also evaluate partition coefficients and their derivatives
		double zidzdJ = 0;
		double zidzdJJ = 0, zidzdJJ1 = 0, zidzdJJ2 = 0;
		vector<double> zidzdc(nsol,0);
		vector<double> zidzdJc(nsol,0), zidzdJc1(nsol,0), zidzdJc2(nsol,0);
		vector< vector<double> > zidzdcc(nsol, vector<double>(nsol,0));
		vector< vector<double> > zidzdcc1(nsol, vector<double>(nsol,0));
		vector<double> zidzdcc2(nsol,0);
		double zidzdcc3 = 0;
		if (den > 0) {
			
			for (isol=0; isol<nsol; ++isol)
				zidzdJ += z[isol]*zz[isol]*dkhdJ[isol]*c[isol];
			zidzdJ = -(dcFdJ+zidzdJ)/den;
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdJJ1 += SQR(z[jsol])*c[jsol]*(z[jsol]*zidzdJ*kappa[jsol]+zz[jsol]*dkhdJ[jsol]);
					zidzdJJ2 += z[jsol]*zz[jsol]*c[jsol]*(zidzdJ*z[jsol]*dkhdJ[jsol]+dkhdJJ[jsol]);
					zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
				}
				zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
				zidzdcc3 += pow(double(z[isol]),3)*kappa[isol]*c[isol];
			}
			zidzdJJ = zidzdJ*(zidzdJ-zidzdJJ1/den)-(dcFdJJ+zidzdJJ2)/den;
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdJc1[isol] += SQR(z[jsol])*c[jsol]*(zidzdc[isol]*z[jsol]*kappa[jsol]+zz[jsol]*dkhdc[jsol][isol]);
					zidzdJc2[isol] += z[jsol]*zz[jsol]*c[jsol]*(zidzdc[isol]*z[jsol]*dkhdJ[jsol]+dkhdJc[jsol][isol]);
					zidzdcc2[isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
					for (ksol=0; ksol<nsol; ++ksol)
						zidzdcc1[isol][jsol] += z[ksol]*zz[ksol]*c[ksol]*dkhdcc[ksol][isol][jsol];
				}
				zidzdJc[isol] = zidzdJ*(zidzdc[isol]-(SQR(z[isol])*kappa[isol] + zidzdJc1[isol])/den)
				-(z[isol]*zz[isol]*dkhdJ[isol] + zidzdJc2[isol])/den;
			}
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdcc[isol][jsol] = zidzdc[isol]*zidzdc[jsol]*(1 - zidzdcc3/den)
					- zidzdcc1[isol][jsol]/den
					- z[isol]*(z[isol]*kappa[isol]*zidzdc[jsol]+zz[isol]*dkhdc[isol][jsol])/den
					- z[jsol]*(z[jsol]*kappa[jsol]*zidzdc[isol]+zz[jsol]*dkhdc[jsol][isol])/den
					- zidzdc[jsol]*zidzdcc2[isol]/den
					- zidzdc[isol]*zidzdcc2[jsol]/den;
				}
			}
		}
		
		vector<double> dkdJ(nsol);
		vector< vector<double> > dkdc(nsol, vector<double>(nsol));
		vector<double> dkdJJ(nsol);
		vector< vector<double> > dkdJc(nsol, vector<double>(nsol));
		vector< vector< vector<double> > > dkdcc(nsol, dkdc);	// use dkhdc for initialization only
		
		for (isol=0; isol<nsol; ++isol) {
			dkdJ[isol] = zz[isol]*dkhdJ[isol]+z[isol]*kappa[isol]*zidzdJ;
			dkdJJ[isol] = zz[isol]*dkhdJJ[isol]+2*z[isol]*zz[isol]*dkhdJ[isol]*zidzdJ
			+z[isol]*kappa[isol]*((z[isol]-1)*SQR(zidzdJ)+zidzdJJ);
			for (jsol=0; jsol<nsol; ++jsol) {
				dkdc[isol][jsol] = zz[isol]*dkhdc[isol][jsol]+z[isol]*kappa[isol]*zidzdc[jsol];
				dkdJc[isol][jsol] = zz[isol]*dkhdJc[isol][jsol]
				+z[isol]*zz[isol]*(dkhdJ[isol]*zidzdc[jsol]+dkhdc[isol][jsol]*zidzdJ)
				+z[isol]*kappa[isol]*((z[isol]-1)*zidzdc[jsol]*zidzdJ+zidzdJc[jsol]);
				for (ksol=0; ksol<nsol; ++ksol) {
					dkdcc[isol][jsol][ksol] = zz[isol]*(dkhdcc[isol][jsol][ksol]
														+z[isol]*(dkhdc[isol][jsol]*zidzdc[ksol]
																  +dkhdc[isol][ksol]*zidzdc[jsol]))
					+z[isol]*kappa[isol]*((z[isol]-1)*zidzdc[ksol]*zidzdc[jsol]+zidzdcc[jsol][ksol]);
				}
			}
		}
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		double dpdJJ = -2*phis/(J*J);
		
		// evaluate the osmotic coefficient
		double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
		
		// evaluate the permeability
		mat3ds K = pm->m_pPerm->Permeability(mp);
		tens4ds dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
		
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
		if (pm->m_pSupp) {
			phiwhat = pm->m_pSupp->Supply(mp);
			Phie = pm->m_pSupp->Tangent_Supply_Strain(mp);
			Phip = pm->m_pSupp->Tangent_Supply_Pressure(mp);
		}

		for (isol=0; isol<nsol; ++isol) {
			// evaluate the permeability derivatives
			dKdc[isol] = pm->m_pPerm->Tangent_Permeability_Concentration(mp,isol);
			
			// evaluate the diffusivity tensor and its derivatives
			D[isol] = pm->m_pSolute[isol]->m_pDiff->Diffusivity(mp);
			dDdE[isol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(mp);
			
			// evaluate the solute free diffusivity
			D0[isol] = pm->m_pSolute[isol]->m_pDiff->Free_Diffusivity(mp);
			
			// evaluate the derivative of the osmotic coefficient
			dodc[isol] = pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,isol);
			
			// evaluate the stress tangent with concentration
			dTdc[isol] = pm->m_pSolid->Tangent_Concentration(mp,isol);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
			}

			// evaluate the solvent supply tangent with concentration
			if (pm->m_pSupp) Phic[isol] = pm->m_pSupp->Tangent_Supply_Concentration(mp,isol);
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
				qpu = -gradN[j]*(divv+1.0/dt)-gradv.transpose()*gradN[j];
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
				ke[ndpn*i+3][ndpn*j+3] -= gradN[i]*(Ke*gradN[j])*(tmp*dt);
				
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
				}
				
				for (isol=0; isol<nsol; ++isol) {
					
					// calculate the kcu matrix
					vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i]
							+ qcu[isol]*H[i])*(tmp*dt);
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
					ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol])*(tmp*dt);
					
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
					}
				}

				// calculate the kcc matrix
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
												   + H[i]*qcc[isol][jsol])*(tmp*dt);
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
	int i, j, isol, jsol, ksol, n;
	
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
	
	FEMesh& mesh = *m_pMesh;
	
	// get the element's material
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (GetMaterial()); assert(pm);
	const int nsol = pm->m_pSolute.size();
	int ndpn = 4+nsol;
	vector<int> sid(nsol);
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->m_pSolute[isol]->GetSoluteID();
	
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
		double J = ept.J;
		
		// get the fluid flux and pressure gradient
		vec3d w = ppt.m_w;
		vec3d gradp = ppt.m_gradp;
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		double dcFdJJ = 2*cF/SQR(J-phi0);
		
		vector<double> c(nsol);
		vector<vec3d> gradc(nsol);
		vector<int> z(nsol);
		vector<double> khat(nsol);
		vector<double> dkhdJ(nsol);
		vector<double> dkhdJJ(nsol);
		vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
		vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
		vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);	// use dkhdc for initialization only
		
		double zeta = pm->ElectricPotential(mp, true);
		vector<double> zz(nsol);
		vector<double> kappa(nsol);
		double den = 0;
		
		for (isol=0; isol<nsol; ++isol) {
			// get the effective concentration, its gradient and its time derivative
			c[isol] = spt.m_c[isol];
			gradc[isol] = spt.m_gradc[isol];
			
			// get the charge number
			z[isol] = pm->m_pSolute[isol]->ChargeNumber();
			
			// evaluate the solubility and its derivatives w.r.t. J and c
			khat[isol] = pm->m_pSolute[isol]->m_pSolub->Solubility(mp);
			dkhdJ[isol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
			dkhdJJ[isol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Strain(mp);
			for (jsol=0; jsol<nsol; ++jsol) {
				dkhdc[isol][jsol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
				dkhdJc[isol][jsol] = pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,jsol);
				for (ksol=0; ksol<nsol; ++ksol) {
					dkhdcc[isol][jsol][ksol] = 
					pm->m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,jsol,ksol);
				}
			}
			zz[isol] = pow(zeta, z[isol]);
			kappa[isol] = zz[isol]*khat[isol];
			den += SQR(z[isol])*kappa[isol]*c[isol];
		}
		
		// evaluate electric potential (nondimensional exponential form) and its derivatives
		// also evaluate partition coefficients and their derivatives
		double zidzdJ = 0;
		double zidzdJJ = 0, zidzdJJ1 = 0, zidzdJJ2 = 0;
		vector<double> zidzdc(nsol,0);
		vector<double> zidzdJc(nsol,0), zidzdJc1(nsol,0), zidzdJc2(nsol,0);
		vector< vector<double> > zidzdcc(nsol, vector<double>(nsol,0)), zidzdcc1(nsol, vector<double>(nsol,0));
		vector<double> zidzdcc2(nsol,0);
		double zidzdcc3 = 0;
		if (den > 0) {
			
			for (isol=0; isol<nsol; ++isol)
				zidzdJ += z[isol]*zz[isol]*dkhdJ[isol]*c[isol];
			zidzdJ = -(dcFdJ+zidzdJ)/den;
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdJJ1 += SQR(z[jsol])*c[jsol]*(z[jsol]*zidzdJ*kappa[jsol]+zz[jsol]*dkhdJ[jsol]);
					zidzdJJ2 += z[jsol]*zz[jsol]*c[jsol]*(zidzdJ*z[jsol]*dkhdJ[jsol]+dkhdJJ[jsol]);
					zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
				}
				zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
				zidzdcc3 += pow(double(z[isol]),3)*kappa[isol]*c[isol];
			}
			zidzdJJ = zidzdJ*(zidzdJ-zidzdJJ1/den)-(dcFdJJ+zidzdJJ2)/den;
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdJc1[isol] += SQR(z[jsol])*c[jsol]*(zidzdc[isol]*z[jsol]*kappa[jsol]+zz[jsol]*dkhdc[jsol][isol]);
					zidzdJc2[isol] += z[jsol]*zz[jsol]*c[jsol]*(zidzdc[isol]*z[jsol]*dkhdJ[jsol]+dkhdJc[jsol][isol]);
					zidzdcc2[isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
					for (ksol=0; ksol<nsol; ++ksol)
						zidzdcc1[isol][jsol] += z[ksol]*zz[ksol]*c[ksol]*dkhdcc[ksol][isol][jsol];
				}
				zidzdJc[isol] = zidzdJ*(zidzdc[isol]-(SQR(z[isol])*kappa[isol] + zidzdJc1[isol])/den)
				-(z[isol]*zz[isol]*dkhdJ[isol] + zidzdJc2[isol])/den;
			}
			
			for (isol=0; isol<nsol; ++isol) {
				for (jsol=0; jsol<nsol; ++jsol) {
					zidzdcc[isol][jsol] = zidzdc[isol]*zidzdc[jsol]*(1 - zidzdcc3/den)
					- zidzdcc1[isol][jsol]/den
					- z[isol]*(z[isol]*kappa[isol]*zidzdc[jsol]+zz[isol]*dkhdc[isol][jsol])/den
					- z[jsol]*(z[jsol]*kappa[jsol]*zidzdc[isol]+zz[jsol]*dkhdc[jsol][isol])/den
					- zidzdc[jsol]*zidzdcc2[isol]/den
					- zidzdc[isol]*zidzdcc2[jsol]/den;
				}
			}
		}
		
		vector<double> dkdJ(nsol);
		vector< vector<double> > dkdc(nsol, vector<double>(nsol));
		vector<double> dkdJJ(nsol);
		vector< vector<double> > dkdJc(nsol, vector<double>(nsol));
		vector< vector< vector<double> > > dkdcc(nsol, dkdc);	// use dkdc for initialization only
		
		for (isol=0; isol<nsol; ++isol) {
			dkdJ[isol] = zz[isol]*dkhdJ[isol]+z[isol]*kappa[isol]*zidzdJ;
			dkdJJ[isol] = zz[isol]*dkhdJJ[isol]+2*z[isol]*zz[isol]*dkhdJ[isol]*zidzdJ
			+z[isol]*kappa[isol]*((z[isol]-1)*SQR(zidzdJ)+zidzdJJ);
			for (jsol=0; jsol<nsol; ++jsol) {
				dkdc[isol][jsol] = zz[isol]*dkhdc[isol][jsol]+z[isol]*kappa[isol]*zidzdc[jsol];
				dkdJc[isol][jsol] = zz[isol]*dkhdJc[isol][jsol]
				+z[isol]*zz[isol]*(dkhdJ[isol]*zidzdc[jsol]+dkhdc[isol][jsol]*zidzdJ)
				+z[isol]*kappa[isol]*((z[isol]-1)*zidzdc[jsol]*zidzdJ+zidzdJc[jsol]);
				for (ksol=0; ksol<nsol; ++ksol) {
					dkdcc[isol][jsol][ksol] = zz[isol]*(dkhdcc[isol][jsol][ksol]
														+z[isol]*(dkhdc[isol][jsol]*zidzdc[ksol]
																  +dkhdc[isol][ksol]*zidzdc[jsol]))
					+z[isol]*kappa[isol]*((z[isol]-1)*zidzdc[ksol]*zidzdc[jsol]+zidzdcc[jsol][ksol]);
				}
			}
		}
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		
		// evaluate the osmotic coefficient
		double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
		
		// evaluate the permeability
		mat3ds K = pm->m_pPerm->Permeability(mp);
		tens4ds dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
		
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
		if (pm->m_pSupp) {
			phiwhat = pm->m_pSupp->Supply(mp);
			Phie = pm->m_pSupp->Tangent_Supply_Strain(mp);
			Phip = pm->m_pSupp->Tangent_Supply_Pressure(mp);
		}

		for (isol=0; isol<nsol; ++isol) {
			// evaluate the permeability derivatives
			dKdc[isol] = pm->m_pPerm->Tangent_Permeability_Concentration(mp,isol);
			
			// evaluate the diffusivity tensor and its derivatives
			D[isol] = pm->m_pSolute[isol]->m_pDiff->Diffusivity(mp);
			dDdE[isol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(mp);
			
			// evaluate the solute free diffusivity
			D0[isol] = pm->m_pSolute[isol]->m_pDiff->Free_Diffusivity(mp);
			
			// evaluate the derivative of the osmotic coefficient
			dodc[isol] = pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,isol);
			
			// evaluate the stress tangent with concentration
			dTdc[isol] = pm->m_pSolid->Tangent_Concentration(mp,isol);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
			}

			// evaluate the solvent supply tangent with concentration
			if (pm->m_pSupp) Phic[isol] = pm->m_pSupp->Tangent_Supply_Concentration(mp,isol);
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
		vec3d vtmp,gp;
		vector<vec3d> gc(nsol),wc(nsol),jce(nsol);
		vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
		mat3d wu, jue;
		vector<mat3d> ju(nsol);
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
				vtmp = (wu.transpose()*gradN[i])*(tmp*dt);
				ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
				
				// calculate the kup matrix
				vtmp = -gradN[i]*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
				
				// calculate the kpp matrix
				ke[ndpn*i+3][ndpn*j+3] -= gradN[i]*(Ke*gradN[j])*(tmp*dt);
				
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
					ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol])*(tmp*dt);
					
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
					}
				}
				
				// calculate the kcc matrix
				for (isol=0; isol<nsol; ++isol) {
					for (jsol=0; jsol<nsol; ++jsol) {
						ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty))*(tmp*dt);
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
	// NOTE: hard-coding of hex elements!
	// Gx = dH/dx
	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];
	
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
	int i, j, k, n;
	int nint, neln;
	double* gw;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES];
	
	FEMesh& mesh = *m_pMesh;

	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());

	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the multiphasic material
		FEMultiphasic* pmb = dynamic_cast<FEMultiphasic*>(pm);
		assert(pmb);
		const int nsol = pmb->m_pSolute.size();
		double ct[MAX_CDOFS][FEElement::MAX_NODES];
		vector<int> sid(nsol);
		for (j=0; j<nsol; ++j) sid[j] = pmb->m_pSolute[j]->GetSoluteID();
		
		// get the number of integration points
		nint = el.GaussPoints();
		
		// get the number of nodes
		neln = el.Nodes();
		assert(neln <= 8);
		
		// get the integration weights
		gw = el.GaussWeights();
		
		// get the nodal data
		for (j=0; j<neln; ++j)
		{
			r0[j] = mesh.Node(el.m_node[j]).m_r0;
			rt[j] = mesh.Node(el.m_node[j]).m_rt;
			pn[j] = mesh.Node(el.m_node[j]).m_pt;
			for (k=0; k<nsol; ++k)
				ct[k][j] = mesh.Node(el.m_node[j]).m_ct[sid[k]];
		}
		
		// loop over the integration points and calculate
		// the stress at the integration point
		for (n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// material point coordinates
			// TODO: I'm not entirly happy with this solution
			//		 since the material point coordinates are used by most materials.
			pt.r0 = el.Evaluate(r0, n);
			pt.rt = el.Evaluate(rt, n);
			
			// get the deformation gradient and determinant
			pt.J = defgrad(el, pt.F, n);
			
			// solute-poroelastic data
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// evaluate fluid pressure at gauss-point
			ppt.m_p = el.Evaluate(pn, n);
			
			// calculate the gradient of p at gauss-point
			ppt.m_gradp = gradient(el, pn, n);
			
			for (k=0; k<nsol; ++k) {
				// evaluate effective solute concentrations at gauss-point
				spt.m_c[k] = el.Evaluate(ct[k], n);
				// calculate the gradient of c at gauss-point
				spt.m_gradc[k] = gradient(el, ct[k], n);
			}
			
			// for multiphasic materials also update the porosity, fluid and solute fluxes
			// and evaluate the actual fluid pressure and solute concentration
			ppt.m_w = pmb->FluidFlux(mp);
			spt.m_psi = pmb->ElectricPotential(mp);
			for (k=0; k<nsol; ++k) {
				spt.m_ca[k] = pmb->Concentration(mp,k);
				spt.m_j[k] = pmb->SoluteFlux(mp,k);
			}
			ppt.m_pa = pmb->Pressure(mp);
			spt.m_cF = pmb->FixedChargeDensity(mp);
			spt.m_Ie = pmb->CurrentDensity(mp);
			
			pt.s = pmb->Stress(mp);
		}
	}
}
