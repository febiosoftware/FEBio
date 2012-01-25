#include "stdafx.h"
#include "FETriphasicDomain.h"
#include "FESolidSolver.h"
#include "FECore/FEMaterial.h"
#include "FEMicroMaterial.h"
#include "FEBioLib/FETriphasic.h"
#include "FEBioLib/log.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
bool FETriphasicDomain::Initialize(FEModel &mdl)
{
	// initialize base class
	FEElasticSolidDomain::Initialize(mdl);
	
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		assert(!el.IsRigid());
		
		assert(el.Type() != FE_UDGHEX);
		
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(mdl.GetMaterial(el.GetMatID()));
		
		assert(dynamic_cast<FETriphasic*>(pm) != 0);
		
		// get the triphasic material
		FETriphasic* pmb = dynamic_cast<FETriphasic*>(pm);
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InitElements()
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
//! calculate internal equivalent nodal forces
void FETriphasicDomain::InternalForces(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	// local element force vector
	vector<double> fl;
	
	vector<int> elm;
	
	// this element should not be rigid
	assert(!el.IsRigid());
	
	//! this element should not be UDG
	assert(el.Type() != FE_UDGHEX);
	
	// unpack the element
	UnpackLM(el, elm);
	
	// get the element force vector and initialize it to zero
	int neln = el.Nodes();
	int ndpn = 6;
	int ndof = ndpn*neln;
	fe.assign(ndof, 0);
	fl.assign(3*neln, 0);
	
	// calculate internal force vector
	// (This function is inherited from FEElasticSolidDomain)
	FEElasticSolidDomain::InternalForces(el, fl);
	
	// copy fl into fe
	int i;
	for (i=0; i<neln; ++i) {
		fe[ndpn*i  ] = fl[3*i  ];
		fe[ndpn*i+1] = fl[3*i+1];
		fe[ndpn*i+2] = fl[3*i+2];
	}
	
	// calculate fluid internal work
	InternalFluidWork(fem, el, fl);
	
	// copy fl into fe
	for (i=0; i<neln; ++i) {
		fe[ndpn*i+3] = fl[i];
	}
	
	// calculate cation internal work
	InternalSoluteWork(fem, el, fl, 0);
	
	// copy fl into fe
	for (i=0; i<neln; ++i) {
		fe[ndpn*i+4] = fl[i];
	}
	
	// calculate anion internal work
	InternalSoluteWork(fem, el, fl, 1);
	
	// copy fl into fe
	for (i=0; i<neln; ++i) {
		fe[ndpn*i+5] = fl[i];
	}
	
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Residual(FENLSolver* psolver, vector<double>& R)
{
	int i, j;
	
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());
	
	// make sure we are in poro-solute mode
	assert(fem.GetCurrentStep()->m_nModule == FE_TRIPHASIC);
	
	// element force vector
	vector<double> fe;
	
	vector<int> elm;
	
	int NE = m_Elem.size();
	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE) {
		for (i=0; i<NE; ++i)
		{
			// get the element
			FESolidElement& el = m_Elem[i];
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
			//! this element should not be UDG
			assert(el.Type() != FE_UDGHEX);
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::InternalForces(el, fe);
			
			// apply body forces
			// TODO: can we calculate body-forces with our formulation
			//       of biphasic theory
			/*
			 if (fem.UseBodyForces())
			 {
			 BodyForces(fem, el, fe);
			 }
			 */
			
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			FEMaterial* pm = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pm) != 0);
			
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
			
			// calculate cation internal work
			InternalSoluteWorkSS(fem, el, fe, 0);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[11*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate anion internal work
			InternalSoluteWorkSS(fem, el, fe, 1);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[12*neln+j];
				if (J >= 0) R[J] += fe[j];
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
			assert(el.Type() != FE_UDGHEX);
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::InternalForces(el, fe);
			
			// apply body forces
			// TODO: can we calculate body-forces with our formulation
			//       of biphasic theory
			/*
			 if (fem.UseBodyForces())
			 {
			 BodyForces(fem, el, fe);
			 }
			 */
			
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			FEMaterial* pm = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pm) != 0);
			
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
			
			// calculate cation internal work
			InternalSoluteWork(fem, el, fe, 0);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[11*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate anion internal work
			InternalSoluteWork(fem, el, fe, 1);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[12*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FETriphasicDomain::InternalFluidWork(FEM& fem, FESolidElement& el, vector<double>& fe)
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
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d rp[8];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
	}
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.GetCurrentStep()->m_dt;
	
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
		double J = ept.J;
		
		// and then finally
		double divv = ((J-Jp)/dt)/J;
		
		// get the flux
		vec3d& w = ppt.m_w;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z - divv*H[i])*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! for a steady-state analysis (zero solid velocity)
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FETriphasicDomain::InternalFluidWorkSS(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt;
	double Gx, Gy, Gz;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEBiphasicMaterialPoint& ppt = *(el.m_State[n]->ExtractData<FEBiphasicMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
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
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z)*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FETriphasicDomain::InternalSoluteWork(FEM& fem, FESolidElement& el, vector<double>& fe, const int ion)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	double *Gr, *Gs, *Gt, *H;
	double *Grr, *Gsr, *Gtr, *Grs, *Gss, *Gts, *Grt, *Gst, *Gtt;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d r0[8], rt[8], rp[8], vt[8];
	double cp[2][8];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[0][i] = mesh.Node(el.m_node[i]).m_cp[0];
		cp[1][i] = mesh.Node(el.m_node[i]).m_cp[1];
		vt[i] = mesh.Node(el.m_node[i]).m_vt;
	}
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESaltMaterialPoint& spt = *(el.m_State[n]->ExtractData<FESaltMaterialPoint>());
		
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
		vec3d gradJ(0);
		double cprev[2] = {0,0};
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		Grr = el.Grr(n); Grs = el.Grs(n); Grt = el.Grt(n);
		Gsr = el.Gsr(n); Gss = el.Gss(n); Gst = el.Gst(n);
		Gtr = el.Gtr(n); Gts = el.Gts(n); Gtt = el.Gtt(n);
		
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
			
			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
			
			// calculate gradJ
			gradJ += (g1*Grr[i] + g2*Grs[i] + g3*Grt[i])*(rt[i]*g1-r0[i]*G1)
			+ (g1*Gsr[i] + g2*Gss[i] + g3*Gst[i])*(rt[i]*g2-r0[i]*G2)
			+ (g1*Gtr[i] + g2*Gts[i] + g3*Gtt[i])*(rt[i]*g3-r0[i]*G3);
			
			// calculate effective concentration at previous time step
			cprev[0] += cp[0][i]*H[i];
			cprev[1] += cp[1][i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.J;
		double dJdt = (J-Jp)/dt;
		gradJ *= J;

		// and then finally
		double divv = dJdt/J;
		
		// get the solute flux
		vec3d& j = spt.m_j[ion];
		// get the effective concentration
		double c[2] = {spt.m_c[0],spt.m_c[1]};

		// get the charge number
		double z[2] = {pm->m_pSolute[0]->ChargeNumber(),
			pm->m_pSolute[1]->ChargeNumber()};
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		
		// evaluate the solubility and its derivatives w.r.t. J and c
		double khat[2] = {
			pm->m_pSolute[0]->m_pSolub->Solubility(mp),
			pm->m_pSolute[1]->m_pSolub->Solubility(mp)};
		double dkhdJ[2] = {
			pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain(mp),
			pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain(mp)};
		double dkhdc[2][2] = {
			{pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,1)},
			{pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,1)}};
		
		// evaluate electric potential (nondimensional exponential form) and its derivatives
		// also evaluate partition coefficients and their derivatives
		double zeta = pm->ElectricPotential(mp, true);
		double zz[2] = {pow(zeta, z[0]), pow(zeta, z[1])};
		double kappa[2] = {zz[0]*khat[0], zz[1]*khat[1]};
		double den = SQR(z[0])*kappa[0]*c[0]+SQR(z[1])*kappa[1]*c[1];
		double zidzdJ = 0;
		double zidzdc[2] = {0,0};
		if (den > 0) {
			zidzdJ = -(dcFdJ+z[0]*zz[0]*dkhdJ[0]*c[0]
					   +z[1]*zz[1]*dkhdJ[1]*c[1])/den;
			zidzdc[0] = -(z[0]*kappa[0]
						  +z[0]*zz[0]*dkhdc[0][0]*c[0]
						  +z[1]*zz[1]*dkhdc[1][0]*c[1])/den;
			zidzdc[1] = -(z[1]*kappa[1]
						  +z[0]*zz[0]*dkhdc[0][1]*c[0]
						  +z[1]*zz[1]*dkhdc[1][1]*c[1])/den;
		}
		double dkdJ[2] = {zz[0]*dkhdJ[0]+z[0]*kappa[0]*zidzdJ,
			zz[1]*dkhdJ[1]+z[1]*kappa[1]*zidzdJ};
		double dkdc[2][2] = {{zz[0]*dkhdc[0][0]+z[0]*kappa[0]*zidzdc[0],
			zz[0]*dkhdc[0][1]+z[0]*kappa[0]*zidzdc[1]},
			{zz[1]*dkhdc[1][0]+z[1]*kappa[1]*zidzdc[0],
				zz[1]*dkhdc[1][1]+z[1]*kappa[1]*zidzdc[1]}};
		
		// evaluate the porosity, its derivative w.r.t. J, and its gradient
		double phiw = pm->Porosity(mp);
		double dpdJ = (1. - phiw)/J;
		vec3d gradp = gradJ*dpdJ;
		// evaluate time derivatives of concentration, solubility and porosity
		double dcdt[2] = {(c[0] - cprev[0])/dt, (c[1] - cprev[1])/dt};
		double dkdt[2] = {dkdJ[0]*dJdt + dkdc[0][0]*dcdt[0] + dkdc[0][1]*dcdt[1],
			dkdJ[1]*dJdt + dkdc[1][0]*dcdt[0] + dkdc[1][1]*dcdt[1]};
		double dpdt = dpdJ*dJdt;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*j.x+B2[i]*j.y+B3[i]*j.z 
						 - H[i]*(dpdt*kappa[ion]*c[ion]+phiw*dkdt[ion]*c[ion]
								 +phiw*kappa[ion]*dcdt[ion]+phiw*kappa[ion]*c[ion]*divv)
						 )*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FETriphasicDomain::InternalSoluteWorkSS(FEM& fem, FESolidElement& el, vector<double>& fe, const int ion)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt;
	double Gx, Gy, Gz;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FESaltMaterialPoint& spt = *(el.m_State[n]->ExtractData<FESaltMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
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
		
		// get the solute flux
		vec3d& j = spt.m_j[ion];
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*j.x+B2[i]*j.y+B3[i]*j.z)*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------

void FETriphasicDomain::StiffnessMatrix(FENLSolver* psolver)
{
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());
	
	// element stiffness matrix
	matrix ke;
	
	vector<int> elm;
	
	// repeat over all solid elements
	int NE = m_Elem.size();
	
	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE) {
		for (int iel=0; iel<NE; ++iel)
		{
			FESolidElement& el = m_Elem[iel];
			
			// this element should not be rigid
			assert(!el.IsRigid());
			
			UnpackLM(el, elm);
			
			// get the elements material
			FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pmat) != 0);
			
			// allocate stiffness matrix
			int neln = el.Nodes();
			int ndpn = 6;
			int ndof = neln*ndpn;
			ke.resize(ndof, ndof);
			
			// calculate the element stiffness matrix
			ElementTriphasicStiffnessSS(fem, el, ke);
			
			// TODO: the problem here is that the LM array that is returned by the UnpackLM
			// function does not give the equation numbers in the right order. For this reason we
			// have to create a new lm array and place the equation numbers in the right order.
			// What we really ought to do is fix the UnpackLM function so that it returns
			// the LM vector in the right order for solute-solid elements.
			vector<int> lm(ndof);
			for (int i=0; i<neln; ++i)
			{
				lm[ndpn*i  ] = elm[3*i];
				lm[ndpn*i+1] = elm[3*i+1];
				lm[ndpn*i+2] = elm[3*i+2];
				lm[ndpn*i+3] = elm[3*neln+i];
				lm[ndpn*i+4] = elm[11*neln+i];
				lm[ndpn*i+5] = elm[12*neln+i];
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
			FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pmat) != 0);
			
			// allocate stiffness matrix
			int neln = el.Nodes();
			int ndpn = 6;
			int ndof = neln*ndpn;
			ke.resize(ndof, ndof);
			
			// calculate the element stiffness matrix
			ElementTriphasicStiffness(fem, el, ke);
			
			// TODO: the problem here is that the LM array that is returned by the UnpackLM
			// function does not give the equation numbers in the right order. For this reason we
			// have to create a new lm array and place the equation numbers in the right order.
			// What we really ought to do is fix the UnpackLM function so that it returns
			// the LM vector in the right order for solute-solid elements.
			vector<int> lm(ndof);
			for (int i=0; i<neln; ++i)
			{
				lm[ndpn*i  ] = elm[3*i];
				lm[ndpn*i+1] = elm[3*i+1];
				lm[ndpn*i+2] = elm[3*i+2];
				lm[ndpn*i+3] = elm[3*neln+i];
				lm[ndpn*i+4] = elm[11*neln+i];
				lm[ndpn*i+5] = elm[12*neln+i];
			}
			
			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FETriphasicDomain::ElementTriphasicStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	int ndpn = 6;
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	vector<vec3d> gradN(neln);
	double tmp;
	
	// gauss-weights
	double* gw = el.GaussWeights();
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d r0[8], rt[8], rp[8], v[8];
	double cp[2][8];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[0][i] = mesh.Node(el.m_node[i]).m_cp[0];
		cp[1][i] = mesh.Node(el.m_node[i]).m_cp[1];
		v[i]  = mesh.Node(el.m_node[i]).m_vt;
	}
	
	// zero stiffness matrix
	ke.zero();
	
	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(fem, el, ks);
	
	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[ndpn*i  ][ndpn*j] = ks[3*i  ][3*j  ]; ke[ndpn*i  ][ndpn*j+1] = ks[3*i  ][3*j+1]; ke[ndpn*i  ][ndpn*j+2] = ks[3*i  ][3*j+2];
			ke[ndpn*i+1][ndpn*j] = ks[3*i+1][3*j  ]; ke[ndpn*i+1][ndpn*j+1] = ks[3*i+1][3*j+1]; ke[ndpn*i+1][ndpn*j+2] = ks[3*i+1][3*j+2];
			ke[ndpn*i+2][ndpn*j] = ks[3*i+2][3*j  ]; ke[ndpn*i+2][ndpn*j+1] = ks[3*i+2][3*j+1]; ke[ndpn*i+2][ndpn*j+2] = ks[3*i+2][3*j+2];
		}
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	// check if we use the symmetric version of the poro-implementation
	bool bsymm = fem.m_bsym_poro;
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(el.m_State[n]->ExtractData<FEBiphasicMaterialPoint>());
		FESaltMaterialPoint& spt = *(el.m_State[n]->ExtractData<FESaltMaterialPoint>());
		
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
		double cprev[2] = {0,0};
		
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
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
			gradN[i] = vec3d(Gx,Gy,Gz);
			
			// calculate effective concentration at previous time step
			cprev[0] += cp[0][i]*H[i];
			cprev[1] += cp[1][i]*H[i];
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
		
		// get the effective concentration, its gradient and its time derivative
		double c[2] = {spt.m_c[0],spt.m_c[1]};
		vec3d gradc[2] = {spt.m_gradc[0],spt.m_gradc[1]};
		double dcdt[2] = {(c[0] - cprev[0])/dt,(c[1] - cprev[1])/dt};
		
		// get the charge number
		double z[2] = {pm->m_pSolute[0]->ChargeNumber(),
			pm->m_pSolute[1]->ChargeNumber()};
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		double dcFdJJ = 2*cF/SQR(J-phi0);
		
		// evaluate the solubility and its derivatives w.r.t. J and c
		double khat[2] = {
			pm->m_pSolute[0]->m_pSolub->Solubility(mp),
			pm->m_pSolute[1]->m_pSolub->Solubility(mp)};
		double dkhdJ[2] = {
			pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain(mp),
			pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain(mp)};
		double dkhdJJ[2] = {
			pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Strain(mp),
			pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Strain(mp)};
		double dkhdc[2][2] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,1)}};
		double dkhdJc[2][2] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,1)}};
		double dkhdcc[2][3] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,1),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,1,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,1),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,1,1)}};
		
		// evaluate electric potential (nondimensional exponential form) and its derivatives
		// also evaluate partition coefficients and their derivatives
		double zeta = pm->ElectricPotential(mp, true);
		double zz[2] = {pow(zeta, z[0]), pow(zeta, z[1])};
		double kappa[2] = {zz[0]*khat[0], zz[1]*khat[1]};
		double den = SQR(z[0])*kappa[0]*c[0]+SQR(z[1])*kappa[1]*c[1];
		double zidzdJ = 0;
		double zidzdJJ = 0;
		double zidzdc[2] = {0,0};
		double zidzdJc[2] = {0,0};
		if (den > 0) {
			zidzdJ = -(dcFdJ+z[0]*zz[0]*dkhdJ[0]*c[0]
					   +z[1]*zz[1]*dkhdJ[1]*c[1])/den;
			zidzdc[0] = -(z[0]*kappa[0]+z[0]*zz[0]*dkhdc[0][0]*c[0]+z[1]*zz[1]*dkhdc[1][0]*c[1])/den;
			zidzdc[1] = -(z[1]*kappa[1]+z[0]*zz[0]*dkhdc[0][1]*c[0]+z[1]*zz[1]*dkhdc[1][1]*c[1])/den;
			zidzdJJ = -(dcFdJJ+z[0]*zz[0]*c[0]*(z[0]*zidzdJ*dkhdJ[0]+dkhdJJ[0])
						+z[1]*zz[1]*c[1]*(z[1]*zidzdJ*dkhdJ[1]+dkhdJJ[1]))/den
			+zidzdJ*(SQR(z[0])*zz[0]*c[0]*(dkhdJ[0]+kappa[0]*(z[0]-1)*zidzdJ)
					 +SQR(z[1])*zz[1]*c[1]*(dkhdJ[1]+kappa[1]*(z[1]-1)*zidzdJ))/den;
			zidzdJc[0] = -(z[0]*zz[0]*dkhdJ[0]+z[0]*zz[0]*c[0]*(z[0]*zidzdc[0]*dkhdJ[0]+dkhdJc[0][0])
						   +z[1]*zz[1]*c[1]*(z[1]*zidzdc[0]*dkhdJ[1]+dkhdJc[1][0]))/den
			+zidzdJ*(SQR(z[0])*kappa[0]+SQR(z[0])*zz[0]*c[0]*(dkhdc[0][0]+kappa[0]*(z[0]-1)*zidzdc[0])
					 +SQR(z[1])*zz[1]*c[1]*(dkhdc[1][0]+kappa[1]*(z[1]-1)*zidzdc[0]))/den;
			zidzdJc[1] = -(z[1]*zz[1]*dkhdJ[1]+z[0]*zz[0]*c[0]*(z[0]*zidzdc[1]*dkhdJ[0]+dkhdJc[0][1])
						   +z[1]*zz[1]*c[1]*(z[1]*zidzdc[1]*dkhdJ[1]+dkhdJc[1][1]))/den
			+zidzdJ*(SQR(z[1])*kappa[1]+SQR(z[0])*zz[0]*c[0]*(dkhdc[0][1]+kappa[0]*(z[0]-1)*zidzdc[1])
					 +SQR(z[1])*zz[1]*c[1]*(dkhdc[1][1]+kappa[1]*(z[1]-1)*zidzdc[1]))/den;
		}
		double dkdJ[2] = {zz[0]*dkhdJ[0]+z[0]*kappa[0]*zidzdJ,
			zz[1]*dkhdJ[1]+z[1]*kappa[1]*zidzdJ};
		double dkdc[2][2] = {{zz[0]*dkhdc[0][0]+z[0]*kappa[0]*zidzdc[0],
			zz[0]*dkhdc[0][1]+z[0]*kappa[0]*zidzdc[1]},
			{zz[1]*dkhdc[1][0]+z[1]*kappa[1]*zidzdc[0],
				zz[1]*dkhdc[1][1]+z[1]*kappa[1]*zidzdc[1]}};
		double dkdJJ[2] = {zz[0]*dkhdJJ[0]+2*z[0]*zz[0]*dkhdJ[0]*zidzdJ
			+z[0]*kappa[0]*((z[0]-1)*SQR(zidzdJ)+zidzdJJ),
			zz[1]*dkhdJJ[1]+2*z[1]*zz[1]*dkhdJ[1]*zidzdJ
			+z[1]*kappa[1]*((z[1]-1)*SQR(zidzdJ)+zidzdJJ)};
		double dkdJc[2][2] = {
			{
				zz[0]*dkhdJc[0][0]+z[0]*zz[0]*(dkhdJ[0]*zidzdc[0]+dkhdc[0][0]*zidzdJ)
				+z[0]*kappa[0]*((z[0]-1)*zidzdc[0]*zidzdJ+zidzdJc[0]),
				zz[0]*dkhdJc[0][1]+z[0]*zz[0]*(dkhdJ[0]*zidzdc[1]+dkhdc[0][1]*zidzdJ)
				+z[0]*kappa[0]*((z[0]-1)*zidzdc[1]*zidzdJ+zidzdJc[1])},
			{
				zz[1]*dkhdJc[1][0]+z[1]*zz[1]*(dkhdJ[1]*zidzdc[0]+dkhdc[1][0]*zidzdJ)
				+z[1]*kappa[1]*((z[1]-1)*zidzdc[0]*zidzdJ+zidzdJc[0]),
				zz[1]*dkhdJc[1][1]+z[1]*zz[1]*(dkhdJ[1]*zidzdc[1]+dkhdc[1][1]*zidzdJ)
				+z[1]*kappa[1]*((z[1]-1)*zidzdc[1]*zidzdJ+zidzdJc[1])}};
		double dkdt[2] = {
			dkdJ[0]*dJdt+dkdc[0][0]*dcdt[0]+dkdc[0][1]*dcdt[1],
			dkdJ[1]*dJdt+dkdc[1][0]*dcdt[0]+dkdc[1][1]*dcdt[1]};
			
		// evaluate the permeability and its derivatives
		mat3ds K = pm->m_pPerm->Permeability(mp);
		tens4ds dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
		mat3ds dKdc[2] = {pm->m_pPerm->Tangent_Permeability_Concentration(mp,0),
			pm->m_pPerm->Tangent_Permeability_Concentration(mp,1)};
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		double dpdJJ = -2*phis/(J*J);
		
		// evaluate the diffusivity tensor and its derivatives
		mat3ds D[2] = {
			pm->m_pSolute[0]->m_pDiff->Diffusivity(mp),
			pm->m_pSolute[1]->m_pDiff->Diffusivity(mp)};
		tens4ds dDdE[2] = {
			pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Strain(mp),
			pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Strain(mp)};
		mat3ds dDdc[2][2] = {
			{
				pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Concentration(mp,1)}};
		
		// evaluate the solute free diffusivity
		double D0[2] = {
			pm->m_pSolute[0]->m_pDiff->Free_Diffusivity(mp),
			pm->m_pSolute[1]->m_pDiff->Free_Diffusivity(mp)};
		double dD0dc[2][2] = {
			{
				pm->m_pSolute[0]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[0]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[1]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,1)}};
			
		// evaluate the osmotic coefficient and its derivatives
		double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
		double dodc[2] = {
			pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,0),
			pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,1)};
		
		// evaluate the stress tangent with concentration
		mat3ds dTdc[2] = {
			pm->m_pSolid->Tangent_Concentration(mp,0), 
			pm->m_pSolid->Tangent_Concentration(mp,1)};
		
		// Miscellaneous constants
		mat3dd I(1);
		double R = pm->m_Rgas;
		double T = pm->m_Tabs;
		
		// evaluate the effective permeability and its derivatives
		mat3ds Ki = K.inverse();
		mat3ds ImD[2] = {I-D[0]/D0[0],I-D[1]/D0[1]};
		mat3ds Ke = (Ki + (ImD[0]*(kappa[0]*c[0]/D0[0])+ImD[1]*(kappa[1]*c[1]/D0[1]))*(R*T/phiw)).inverse();
		tens4ds G = dyad1s(Ki,I) - dyad4s(Ki,I)*2 - ddots(dyad2s(Ki),dKdE)*0.5
		+dyad1s(ImD[0],I)*(R*T*c[0]*J/D0[0]/2/phiw*(dkdJ[0]-kappa[0]/phiw*dpdJ))
		+(dyad1s(I) - dyad4s(I)*2 - dDdE[0]/D0[0])*(R*T*kappa[0]*c[0]/phiw/D0[0])
		+dyad1s(ImD[1],I)*(R*T*c[1]*J/D0[1]/2/phiw*(dkdJ[1]-kappa[1]/phiw*dpdJ))
		+(dyad1s(I) - dyad4s(I)*2 - dDdE[1]/D0[1])*(R*T*kappa[1]*c[1]/phiw/D0[1]);
		tens4ds dKedE = dyad1s(Ke,I) - 2*dyad4s(Ke,I) - ddots(dyad2s(Ke),G)*0.5;
		mat3ds Gc[2] = {
			(ImD[0]*(kappa[0]/D0[0])
			 +ImD[0]*(c[0]/D0[0]*(dkdc[0][0]-kappa[0]/D0[0]*dD0dc[0][0]))
			 +ImD[1]*(c[1]/D0[1]*(dkdc[1][0]-kappa[1]/D0[1]*dD0dc[1][0]))
			 -(dDdc[0][0]-D[0]*(dD0dc[0][0]/D0[0])*(kappa[0]*c[0]/SQR(D0[0])))
			 -(dDdc[1][0]-D[1]*(dD0dc[1][0]/D0[1])*(kappa[1]*c[1]/SQR(D0[1])))
			 )*(R*T/phiw),
			(ImD[1]*(kappa[1]/D0[1])
			 +ImD[0]*(c[0]/D0[0]*(dkdc[0][1]-kappa[0]/D0[0]*dD0dc[0][1]))
			 +ImD[1]*(c[1]/D0[1]*(dkdc[1][1]-kappa[1]/D0[1]*dD0dc[1][1]))
			 -(dDdc[0][1]-D[0]*(dD0dc[0][1]/D0[0])*(kappa[0]*c[0]/SQR(D0[0])))
			 -(dDdc[1][1]-D[1]*(dD0dc[1][1]/D0[1])*(kappa[1]*c[1]/SQR(D0[1])))
			 )*(R*T/phiw)};
		mat3ds dKedc[2] = {-Ke*(-Ki*dKdc[0]*Ki + Gc[0])*Ke,-Ke*(-Ki*dKdc[1]*Ki + Gc[1])*Ke};
		
		// calculate all the matrices
		vec3d vtmp,gp,gc[2],qpu,qcu[2],wc[2],jc[2][2];
		mat3d wu,ju[2];
		double qcc[2][2];
		tmp = detJ*gw[n];
		for (i=0; i<neln; ++i)
		{
			for (j=0; j<neln; ++j)
			{
				// calculate the kpu matrix
				gp = gradp+((D[0]*gradc[0])*(kappa[0]/D0[0])+(D[1]*gradc[1])*(kappa[1]/D0[1]))*R*T;
				wu = vdotTdotv(-gp, dKedE, gradN[j])
				-(((Ke*(D[0]*gradc[0])) & gradN[j])*(J*dkdJ[0] - kappa[0])
				  +Ke*(2*kappa[0]*(gradN[j]*(D[0]*gradc[0]))))*R*T/D0[0]
				- Ke*vdotTdotv(gradc[0], dDdE[0], gradN[j])*(kappa[0]*R*T/D0[0])
				-(((Ke*(D[1]*gradc[1])) & gradN[j])*(J*dkdJ[1] - kappa[1])
				  +Ke*(2*kappa[1]*(gradN[j]*(D[1]*gradc[1]))))*R*T/D0[1]
				- Ke*vdotTdotv(gradc[1], dDdE[1], gradN[j])*(kappa[1]*R*T/D0[1]);
				qpu = -gradN[j]*(divv+1/dt)-gradv.transpose()*gradN[j];
				vtmp = (wu.transpose()*gradN[i] + qpu*H[i])*(tmp*dt);
				ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
				
				// calculate the kcu matrix for the cation
				gc[0] = -gradc[0]*phiw + w*c[0]/D0[0];
				ju[0] = ((D[0]*gc[0]) & gradN[j])*(J*dkdJ[0]) 
				+ vdotTdotv(gc[0], dDdE[0], gradN[j])*kappa[0]
				+ (((D[0]*gradc[0]) & gradN[j])*(-phis)
				   +(D[0]*((gradN[j]*w)*2) - ((D[0]*w) & gradN[j]))*c[0]/D0[0]
				   )*kappa[0]
				+D[0]*wu*(kappa[0]*c[0]/D0[0]);
				qcu[0] = -gradN[j]*(c[0]*dJdt*(2*(dpdJ*kappa[0]+phiw*dkdJ[0]+J*dpdJ*dkdJ[0])
											   +J*(dpdJJ*kappa[0]+phiw*dkdJJ[0]))
									+dcdt[0]*(phiw*kappa[0]+J*dpdJ*kappa[0]+J*phiw*dkdJ[0])
									+c[0]*(dcdt[0]*((phiw+J*dpdJ)*dkdc[0][0]+J*phiw*dkdJc[0][0])
										   +dcdt[1]*((phiw+J*dpdJ)*dkdc[0][1]+J*phiw*dkdJc[0][1])))
				+qpu*(c[0]*(phiw*kappa[0]+J*dpdJ*kappa[0]+J*phiw*dkdJ[0]));
				vtmp = (ju[0].transpose()*gradN[i] + qcu[0]*H[i])*(tmp*dt);
				ke[ndpn*i+4][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+4][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+4][ndpn*j+2] += vtmp.z;

				// calculate the kcu matrix for the anion
				gc[1] = -gradc[1]*phiw + w*c[1]/D0[1];
				ju[1] = ((D[1]*gc[1]) & gradN[j])*(J*dkdJ[1]) 
				+ vdotTdotv(gc[1], dDdE[1], gradN[j])*kappa[1]
				+ (((D[1]*gradc[1]) & gradN[j])*(-phis)
				   +(D[1]*((gradN[j]*w)*2) - ((D[1]*w) & gradN[j]))*c[1]/D0[1]
				   )*kappa[1]
				+D[1]*wu*(kappa[1]*c[1]/D0[1]);
				qcu[1] = -gradN[j]*(c[1]*dJdt*(2*(dpdJ*kappa[1]+phiw*dkdJ[1]+J*dpdJ*dkdJ[1])
											   +J*(dpdJJ*kappa[1]+phiw*dkdJJ[1]))
									+dcdt[1]*(phiw*kappa[1]+J*dpdJ*kappa[1]+J*phiw*dkdJ[1])
									+c[1]*(dcdt[0]*((phiw+J*dpdJ)*dkdc[1][0]+J*phiw*dkdJc[1][0])
										   +dcdt[1]*((phiw+J*dpdJ)*dkdc[1][1]+J*phiw*dkdJc[1][1])))
				+qpu*(c[1]*(phiw*kappa[1]+J*dpdJ*kappa[1]+J*phiw*dkdJ[1]));
				vtmp = (ju[1].transpose()*gradN[i] + qcu[1]*H[i])*(tmp*dt);
				ke[ndpn*i+5][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+5][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+5][ndpn*j+2] += vtmp.z;
				
				// calculate the kup matrix
				vtmp = -gradN[i]*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
				
				// calculate the kpp matrix
				ke[ndpn*i+3][ndpn*j+3] -= gradN[i]*(Ke*gradN[j])*(tmp*dt);
				
				// calculate the kcp matrix for the cation
				ke[ndpn*i+4][ndpn*j+3] -= (gradN[i]*((D[0]*Ke)*gradN[j]))*(kappa[0]*c[0]/D0[0])*(tmp*dt);
				
				// calculate the kcp matrix for the anion
				ke[ndpn*i+5][ndpn*j+3] -= (gradN[i]*((D[1]*Ke)*gradN[j]))*(kappa[1]*c[1]/D0[1])*(tmp*dt);
				
				// calculate the kuc matrix for the cation
				vtmp = (dTdc[0]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[0]
														  +c[0]*(dodc[0]*kappa[0]+osmc*dkdc[0][0])
														  +c[1]*(dodc[0]*kappa[1]+osmc*dkdc[1][0])
														  )))*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+4] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+4] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+4] += vtmp.z;
				
				// calculate the kuc matrix for the anion
				vtmp = (dTdc[1]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[1]
														  +c[0]*(dodc[1]*kappa[0]+osmc*dkdc[0][1])
														  +c[1]*(dodc[1]*kappa[1]+osmc*dkdc[1][1])
														  )))*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+5] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+5] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+5] += vtmp.z;
				
				// calculate the kpc matrix for the cation
				wc[0] = (dKedc[0]*gp)*(-H[j])
				-Ke*((D[0]*gradN[j])*(kappa[0]/D0[0])
					 +((D[0]*(dkdc[0][0]-kappa[0]/D0[0]*dD0dc[0][0])
						+dDdc[0][0]*kappa[0])*gradc[0]
					   +(D[1]*(dkdc[1][0]-kappa[1]/D0[1]*dD0dc[1][0])
						 +dDdc[1][0]*kappa[1])*gradc[1]
					   )*H[j]
					 )*(R*T);
				ke[ndpn*i+3][ndpn*j+4] += (gradN[i]*wc[0])*(tmp*dt);
				
				// calculate the kpc matrix for the anion
				wc[1] = (dKedc[1]*gp)*(-H[j])
				-Ke*((D[1]*gradN[j])*(kappa[1]/D0[1])
					 +((D[0]*(dkdc[0][1]-kappa[0]/D0[0]*dD0dc[0][1])
						+dDdc[0][1]*kappa[0])*gradc[0]
					   +(D[1]*(dkdc[1][1]-kappa[1]/D0[1]*dD0dc[1][1])
						 +dDdc[1][1]*kappa[1])*gradc[1]
					   )*H[j]
					 )*(R*T);
				ke[ndpn*i+3][ndpn*j+5] += (gradN[i]*wc[1])*(tmp*dt);
				
				// calculate the kcc matrix for the cation-cation
				jc[0][0] = (D[0]*(gradN[j]*(-phiw)+w*(H[j]/D0[0])))*kappa[0]
				+((D[0]*dkdc[0][0]+dDdc[0][0]*kappa[0])*(gradc[0]*(-phiw)+w*(c[0]/D0[0])))*H[j]
				+(D[0]*w*(-H[j]*dD0dc[0][0]/D0[0])+D[0]*wc[0])*(kappa[0]*c[0]/D0[0]);
				qcc[0][0] = -H[j]*(((phiw+J*dpdJ)*divv+phiw/dt)*(kappa[0]+c[0]*dkdc[0][0])
								   +phiw*(dkdt[0]+dkdc[0][0]*dcdt[0])
								   +phiw*c[0]*dJdt*dkdJc[0][0]
								   );
				ke[ndpn*i+4][ndpn*j+4] += (gradN[i]*jc[0][0] + H[i]*qcc[0][0])*(tmp*dt);
				
				// calculate the kcc matrix for the anion-anion
				jc[1][1] = (D[1]*(gradN[j]*(-phiw)+w*(H[j]/D0[1])))*kappa[1]
				+((D[1]*dkdc[1][1]+dDdc[1][1]*kappa[1])*(gradc[1]*(-phiw)+w*(c[1]/D0[1])))*H[j]
				+(D[1]*w*(-H[j]*dD0dc[1][1]/D0[1])+D[1]*wc[1])*(kappa[1]*c[1]/D0[1]);
				qcc[1][1] = -H[j]*(((phiw+J*dpdJ)*divv+phiw/dt)*(kappa[1]+c[1]*dkdc[1][1])
								   +phiw*(dkdt[1]+dkdc[1][1]*dcdt[1])
								   +phiw*c[1]*dJdt*dkdJc[1][1]
								   );
				ke[ndpn*i+5][ndpn*j+5] += (gradN[i]*jc[1][1] + H[i]*qcc[1][1])*(tmp*dt);
				
				// calculate the kcc matrix for the cation-anion
				jc[0][1] = 
				((D[0]*dkdc[0][1]+dDdc[0][1]*kappa[0])*(gradc[0]*(-phiw)+w*(c[0]/D0[0])))*H[j]
				+(D[0]*w*(-H[j]*dD0dc[0][1]/D0[0])+D[0]*wc[1])*(kappa[0]*c[0]/D0[0]);
				qcc[0][1] = -H[j]*(((phiw+J*dpdJ)*divv+phiw/dt)*(c[0]*dkdc[0][1])
								   +phiw*(dkdc[0][1]*dcdt[0])
								   +phiw*c[0]*dJdt*dkdJc[0][1]
								   );
				ke[ndpn*i+4][ndpn*j+5] += (gradN[i]*jc[0][1] + H[i]*qcc[0][1])*(tmp*dt);
				
				// calculate the kcc matrix for the anion-cation
				jc[1][0] = 
				((D[1]*dkdc[1][0]+dDdc[1][0]*kappa[1])*(gradc[1]*(-phiw)+w*(c[1]/D0[1])))*H[j]
				+(D[1]*w*(-H[j]*dD0dc[1][0]/D0[1])+D[1]*wc[0])*(kappa[1]*c[1]/D0[1]);
				qcc[1][0] = -H[j]*(((phiw+J*dpdJ)*divv+phiw/dt)*(c[1]*dkdc[1][0])
								   +phiw*(dkdc[1][0]*dcdt[1])
								   +phiw*c[1]*dJdt*dkdJc[1][0]
								   );
				ke[ndpn*i+5][ndpn*j+4] += (gradN[i]*jc[1][0] + H[i]*qcc[1][0])*(tmp*dt);
				
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
bool FETriphasicDomain::ElementTriphasicStiffnessSS(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	int ndpn = 6;
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz;
	
	// jacobian
	double Ji[3][3], detJ;
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	vector<vec3d> gradN(neln);
	double tmp;
	
	// gauss-weights
	double* gw = el.GaussWeights();
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d r0[8], rt[8], rp[8], v[8];
	double cp[2][8];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[0][i] = mesh.Node(el.m_node[i]).m_cp[0];
		cp[1][i] = mesh.Node(el.m_node[i]).m_cp[1];
		v[i]  = mesh.Node(el.m_node[i]).m_vt;
	}
	
	// zero stiffness matrix
	ke.zero();
	
	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(fem, el, ks);
	
	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[ndpn*i  ][ndpn*j] = ks[3*i  ][3*j  ]; ke[ndpn*i  ][ndpn*j+1] = ks[3*i  ][3*j+1]; ke[ndpn*i  ][ndpn*j+2] = ks[3*i  ][3*j+2];
			ke[ndpn*i+1][ndpn*j] = ks[3*i+1][3*j  ]; ke[ndpn*i+1][ndpn*j+1] = ks[3*i+1][3*j+1]; ke[ndpn*i+1][ndpn*j+2] = ks[3*i+1][3*j+2];
			ke[ndpn*i+2][ndpn*j] = ks[3*i+2][3*j  ]; ke[ndpn*i+2][ndpn*j+1] = ks[3*i+2][3*j+1]; ke[ndpn*i+2][ndpn*j+2] = ks[3*i+2][3*j+2];
		}
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	// check if we use the symmetric version of the poro-implementation
	bool bsymm = fem.m_bsym_poro;
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(el.m_State[n]->ExtractData<FEBiphasicMaterialPoint>());
		FESaltMaterialPoint& spt = *(el.m_State[n]->ExtractData<FESaltMaterialPoint>());
		
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
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
			gradN[i] = vec3d(Gx,Gy,Gz);
		}
		
		// next we get the determinant
		double J = ept.J;
		
		// get the fluid flux and pressure gradient
		vec3d w = ppt.m_w;
		vec3d gradp = ppt.m_gradp;
		
		// get the effective concentration and its gradient
		double c[2] = {spt.m_c[0],spt.m_c[1]};
		vec3d gradc[2] = {spt.m_gradc[0],spt.m_gradc[1]};

		// get the charge number
		double z[2] = {pm->m_pSolute[0]->ChargeNumber(),
			pm->m_pSolute[1]->ChargeNumber()};
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		double dcFdJJ = 2*cF/SQR(J-phi0);
		
		// evaluate the solubility and its derivatives w.r.t. J and c
		double khat[2] = {
			pm->m_pSolute[0]->m_pSolub->Solubility(mp),
			pm->m_pSolute[1]->m_pSolub->Solubility(mp)};
		double dkhdJ[2] = {
			pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain(mp),
			pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain(mp)};
		double dkhdJJ[2] = {
			pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Strain(mp),
			pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Strain(mp)};
		double dkhdc[2][2] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration(mp,1)}};
		double dkhdJc[2][2] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,1)}};
		double dkhdcc[2][3] = {
			{
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,0),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,1),
				pm->m_pSolute[0]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,1,1)},
			{
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,0),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,0,1),
				pm->m_pSolute[1]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,1,1)}};
		
		// evaluate electric potential (nondimensional exponential form) and its derivatives
		// also evaluate partition coefficients and their derivatives
		double zeta = pm->ElectricPotential(mp, true);
		double zz[2] = {pow(zeta, z[0]), pow(zeta, z[1])};
		double kappa[2] = {zz[0]*khat[0], zz[1]*khat[1]};
		double den = SQR(z[0])*kappa[0]*c[0]+SQR(z[1])*kappa[1]*c[1];
		double zidzdJ = 0;
		double zidzdJJ = 0;
		double zidzdc[2] = {0,0};
		double zidzdJc[2] = {0,0};
		if (den > 0) {
			zidzdJ = -(dcFdJ+z[0]*zz[0]*dkhdJ[0]*c[0]
					   +z[1]*zz[1]*dkhdJ[1]*c[1])/den;
			zidzdc[0] = -(z[0]*kappa[0]+z[0]*zz[0]*dkhdc[0][0]*c[0]+z[1]*zz[1]*dkhdc[1][0]*c[1])/den;
			zidzdc[1] = -(z[1]*kappa[1]+z[0]*zz[0]*dkhdc[0][1]*c[0]+z[1]*zz[1]*dkhdc[1][1]*c[1])/den;
			zidzdJJ = -(dcFdJJ+z[0]*zz[0]*c[0]*(z[0]*zidzdJ*dkhdJ[0]+dkhdJJ[0])
						+z[1]*zz[1]*c[1]*(z[1]*zidzdJ*dkhdJ[1]+dkhdJJ[1]))/den
			+zidzdJ*(SQR(z[0])*zz[0]*c[0]*(dkhdJ[0]+kappa[0]*(z[0]-1)*zidzdJ)
					 +SQR(z[1])*zz[1]*c[1]*(dkhdJ[1]+kappa[1]*(z[1]-1)*zidzdJ))/den;
			zidzdJc[0] = -(z[0]*zz[0]*dkhdJ[0]+z[0]*zz[0]*c[0]*(z[0]*zidzdc[0]*dkhdJ[0]+dkhdJc[0][0])
						   +z[1]*zz[1]*c[1]*(z[1]*zidzdc[0]*dkhdJ[1]+dkhdJc[1][0]))/den
			+zidzdJ*(SQR(z[0])*kappa[0]+SQR(z[0])*zz[0]*c[0]*(dkhdc[0][0]+kappa[0]*(z[0]-1)*zidzdc[0])
					 +SQR(z[1])*zz[1]*c[1]*(dkhdc[1][0]+kappa[1]*(z[1]-1)*zidzdc[0]))/den;
			zidzdJc[1] = -(z[1]*zz[1]*dkhdJ[1]+z[0]*zz[0]*c[0]*(z[0]*zidzdc[1]*dkhdJ[0]+dkhdJc[0][1])
						   +z[1]*zz[1]*c[1]*(z[1]*zidzdc[1]*dkhdJ[1]+dkhdJc[1][1]))/den
			+zidzdJ*(SQR(z[1])*kappa[1]+SQR(z[0])*zz[0]*c[0]*(dkhdc[0][1]+kappa[0]*(z[0]-1)*zidzdc[1])
					 +SQR(z[1])*zz[1]*c[1]*(dkhdc[1][1]+kappa[1]*(z[1]-1)*zidzdc[1]))/den;
		}
		double dkdJ[2] = {zz[0]*dkhdJ[0]+z[0]*kappa[0]*zidzdJ,
			zz[1]*dkhdJ[1]+z[1]*kappa[1]*zidzdJ};
		double dkdc[2][2] = {{zz[0]*dkhdc[0][0]+z[0]*kappa[0]*zidzdc[0],
			zz[0]*dkhdc[0][1]+z[0]*kappa[0]*zidzdc[1]},
			{zz[1]*dkhdc[1][0]+z[1]*kappa[1]*zidzdc[0],
				zz[1]*dkhdc[1][1]+z[1]*kappa[1]*zidzdc[1]}};
		double dkdJJ[2] = {zz[0]*dkhdJJ[0]+2*z[0]*zz[0]*dkhdJ[0]*zidzdJ
			+z[0]*kappa[0]*((z[0]-1)*SQR(zidzdJ)+zidzdJJ),
			zz[1]*dkhdJJ[1]+2*z[1]*zz[1]*dkhdJ[1]*zidzdJ
			+z[1]*kappa[1]*((z[1]-1)*SQR(zidzdJ)+zidzdJJ)};
		double dkdJc[2][2] = {
			{
				zz[0]*dkhdJc[0][0]+z[0]*zz[0]*(dkhdJ[0]*zidzdc[0]+dkhdc[0][0]*zidzdJ)
				+z[0]*kappa[0]*((z[0]-1)*zidzdc[0]*zidzdJ+zidzdJc[0]),
				zz[0]*dkhdJc[0][1]+z[0]*zz[0]*(dkhdJ[0]*zidzdc[1]+dkhdc[0][1]*zidzdJ)
				+z[0]*kappa[0]*((z[0]-1)*zidzdc[1]*zidzdJ+zidzdJc[1])},
			{
				zz[1]*dkhdJc[1][0]+z[1]*zz[1]*(dkhdJ[1]*zidzdc[0]+dkhdc[1][0]*zidzdJ)
				+z[1]*kappa[1]*((z[1]-1)*zidzdc[0]*zidzdJ+zidzdJc[0]),
				zz[1]*dkhdJc[1][1]+z[1]*zz[1]*(dkhdJ[1]*zidzdc[1]+dkhdc[1][1]*zidzdJ)
				+z[1]*kappa[1]*((z[1]-1)*zidzdc[1]*zidzdJ+zidzdJc[1])}};
		
		// evaluate the permeability and its derivatives
		mat3ds K = pm->m_pPerm->Permeability(mp);
		tens4ds dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
		mat3ds dKdc[2] = {pm->m_pPerm->Tangent_Permeability_Concentration(mp,0),
			pm->m_pPerm->Tangent_Permeability_Concentration(mp,1)};
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double phis = 1. - phiw;
		double dpdJ = phis/J;
		double dpdJJ = -2*phis/(J*J);
		
		// evaluate the diffusivity tensor and its derivatives
		mat3ds D[2] = {
			pm->m_pSolute[0]->m_pDiff->Diffusivity(mp),
			pm->m_pSolute[1]->m_pDiff->Diffusivity(mp)};
		tens4ds dDdE[2] = {
			pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Strain(mp),
			pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Strain(mp)};
		mat3ds dDdc[2][2] = {
			{
				pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[0]->m_pDiff->Tangent_Diffusivity_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[1]->m_pDiff->Tangent_Diffusivity_Concentration(mp,1)}};
		
		// evaluate the solute free diffusivity
		double D0[2] = {
			pm->m_pSolute[0]->m_pDiff->Free_Diffusivity(mp),
			pm->m_pSolute[1]->m_pDiff->Free_Diffusivity(mp)};
		double dD0dc[2][2] = {
			{
				pm->m_pSolute[0]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[0]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,1)},
			{
				pm->m_pSolute[1]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0),
				pm->m_pSolute[1]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,1)}};
		
		// evaluate the osmotic coefficient and its derivatives
		double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
		double dodc[2] = {
			pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,0),
			pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,1)};
		
		// evaluate the stress tangent with concentration
		mat3ds dTdc[2] = {
			pm->m_pSolid->Tangent_Concentration(mp,0),
			pm->m_pSolid->Tangent_Concentration(mp,1)};
		
		// Miscellaneous constants
		mat3dd I(1);
		double R = pm->m_Rgas;
		double T = pm->m_Tabs;
		
		// evaluate the effective permeability and its derivatives
		mat3ds Ki = K.inverse();
		mat3ds ImD[2] = {I-D[0]/D0[0],I-D[1]/D0[1]};
		mat3ds Ke = (Ki + (ImD[0]*(kappa[0]*c[0]/D0[0])+ImD[1]*(kappa[1]*c[1]/D0[1]))*(R*T/phiw)).inverse();
		tens4ds G = dyad1s(Ki,I) - dyad4s(Ki,I)*2 - ddots(dyad2s(Ki),dKdE)*0.5
		+dyad1s(ImD[0],I)*(R*T*c[0]*J/D0[0]/2/phiw*(dkdJ[0]-kappa[0]/phiw*dpdJ))
		+(dyad1s(I) - dyad4s(I)*2 - dDdE[0]/D0[0])*(R*T*kappa[0]*c[0]/phiw/D0[0])
		+dyad1s(ImD[1],I)*(R*T*c[1]*J/D0[1]/2/phiw*(dkdJ[1]-kappa[1]/phiw*dpdJ))
		+(dyad1s(I) - dyad4s(I)*2 - dDdE[1]/D0[1])*(R*T*kappa[1]*c[1]/phiw/D0[1]);
		tens4ds dKedE = dyad1s(Ke,I) - 2*dyad4s(Ke,I) - ddots(dyad2s(Ke),G)*0.5;
		mat3ds Gc[2] = {
			(ImD[0]*(kappa[0]/D0[0])
			 +ImD[0]*(c[0]/D0[0]*(dkdc[0][0]-kappa[0]/D0[0]*dD0dc[0][0]))
			 +ImD[1]*(c[1]/D0[1]*(dkdc[1][0]-kappa[1]/D0[1]*dD0dc[1][0]))
			 -(dDdc[0][0]-D[0]*(dD0dc[0][0]/D0[0])*(kappa[0]*c[0]/SQR(D0[0])))
			 -(dDdc[1][0]-D[1]*(dD0dc[1][0]/D0[1])*(kappa[1]*c[1]/SQR(D0[1])))
			 )*(R*T/phiw),
			(ImD[1]*(kappa[1]/D0[1])
			 +ImD[0]*(c[0]/D0[0]*(dkdc[0][1]-kappa[0]/D0[0]*dD0dc[0][1]))
			 +ImD[1]*(c[1]/D0[1]*(dkdc[1][1]-kappa[1]/D0[1]*dD0dc[1][1]))
			 -(dDdc[0][1]-D[0]*(dD0dc[0][1]/D0[0])*(kappa[0]*c[0]/SQR(D0[0])))
			 -(dDdc[1][1]-D[1]*(dD0dc[1][1]/D0[1])*(kappa[1]*c[1]/SQR(D0[1])))
			 )*(R*T/phiw)};
		mat3ds dKedc[2] = {-Ke*(-Ki*dKdc[0]*Ki + Gc[0])*Ke,-Ke*(-Ki*dKdc[1]*Ki + Gc[1])*Ke};
		
		// calculate all the matrices
		vec3d vtmp,gp,gc[2],wc[2],jc[2][2];
		mat3d wu,ju[2];
		tmp = detJ*gw[n];
		for (i=0; i<neln; ++i)
		{
			for (j=0; j<neln; ++j)
			{
				// calculate the kpu matrix
				gp = gradp+((D[0]*gradc[0])*(kappa[0]/D0[0])+(D[1]*gradc[1])*(kappa[1]/D0[1]))*R*T;
				wu = vdotTdotv(-gp, dKedE, gradN[j])
				-(((Ke*(D[0]*gradc[0])) & gradN[j])*(J*dkdJ[0] - kappa[0])
				  +Ke*(2*kappa[0]*(gradN[j]*(D[0]*gradc[0]))))*R*T/D0[0]
				- Ke*vdotTdotv(gradc[0], dDdE[0], gradN[j])*(kappa[0]*R*T/D0[0])
				-(((Ke*(D[1]*gradc[1])) & gradN[j])*(J*dkdJ[1] - kappa[1])
				  +Ke*(2*kappa[1]*(gradN[j]*(D[1]*gradc[1]))))*R*T/D0[1]
				- Ke*vdotTdotv(gradc[1], dDdE[1], gradN[j])*(kappa[1]*R*T/D0[1]);
				vtmp = (wu.transpose()*gradN[i])*(tmp*dt);
				ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
				
				// calculate the kcu matrix for the cation
				gc[0] = -gradc[0]*phiw + w*c[0]/D0[0];
				ju[0] = ((D[0]*gc[0]) & gradN[j])*(J*dkdJ[0]) 
				+ vdotTdotv(gc[0], dDdE[0], gradN[j])*kappa[0]
				+ (((D[0]*gradc[0]) & gradN[j])*(-phis)
				   +(D[0]*((gradN[j]*w)*2) - ((D[0]*w) & gradN[j]))*c[0]/D0[0]
				   )*kappa[0]
				+D[0]*wu*(kappa[0]*c[0]/D0[0]);
				vtmp = (ju[0].transpose()*gradN[i])*(tmp*dt);
				ke[ndpn*i+4][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+4][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+4][ndpn*j+2] += vtmp.z;
				
				// calculate the kcu matrix for the anion
				gc[1] = -gradc[1]*phiw + w*c[1]/D0[1];
				ju[1] = ((D[1]*gc[1]) & gradN[j])*(J*dkdJ[1]) 
				+ vdotTdotv(gc[1], dDdE[1], gradN[j])*kappa[1]
				+ (((D[1]*gradc[1]) & gradN[j])*(-phis)
				   +(D[1]*((gradN[j]*w)*2) - ((D[1]*w) & gradN[j]))*c[1]/D0[1]
				   )*kappa[1]
				+D[1]*wu*(kappa[1]*c[1]/D0[1]);
				vtmp = (ju[1].transpose()*gradN[i])*(tmp*dt);
				ke[ndpn*i+5][ndpn*j  ] += vtmp.x;
				ke[ndpn*i+5][ndpn*j+1] += vtmp.y;
				ke[ndpn*i+5][ndpn*j+2] += vtmp.z;
				
				// calculate the kup matrix
				vtmp = -gradN[i]*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
				
				// calculate the kpp matrix
				ke[ndpn*i+3][ndpn*j+3] -= gradN[i]*(Ke*gradN[j])*(tmp*dt);
				
				// calculate the kcp matrix for the cation
				ke[ndpn*i+4][ndpn*j+3] -= (gradN[i]*((D[0]*Ke)*gradN[j]))*(kappa[0]*c[0]/D0[0])*(tmp*dt);
				
				// calculate the kcp matrix for the anion
				ke[ndpn*i+5][ndpn*j+3] -= (gradN[i]*((D[1]*Ke)*gradN[j]))*(kappa[1]*c[1]/D0[1])*(tmp*dt);
				
				// calculate the kuc matrix for the cation
				vtmp = (dTdc[0]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[0]
														  +c[0]*(dodc[0]*kappa[0]+osmc*dkdc[0][0])
														  +c[1]*(dodc[0]*kappa[1]+osmc*dkdc[1][0])
														  )))*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+4] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+4] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+4] += vtmp.z;
				
				// calculate the kuc matrix for the anion
				vtmp = (dTdc[1]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[1]
														  +c[0]*(dodc[1]*kappa[0]+osmc*dkdc[0][1])
														  +c[1]*(dodc[1]*kappa[1]+osmc*dkdc[1][1])
														  )))*H[j]*tmp;
				ke[ndpn*i  ][ndpn*j+5] += vtmp.x;
				ke[ndpn*i+1][ndpn*j+5] += vtmp.y;
				ke[ndpn*i+2][ndpn*j+5] += vtmp.z;
				
				// calculate the kpc matrix for the cation
				wc[0] = (dKedc[0]*gp)*(-H[j])
				-Ke*((D[0]*gradN[j])*(kappa[0]/D0[0])
					 +((D[0]*(dkdc[0][0]-kappa[0]/D0[0]*dD0dc[0][0])
						+dDdc[0][0]*kappa[0])*gradc[0]
					   +(D[1]*(dkdc[1][0]-kappa[1]/D0[1]*dD0dc[1][0])
						 +dDdc[1][0]*kappa[1])*gradc[1]
					   )*H[j]
					 )*(R*T);
				ke[ndpn*i+3][ndpn*j+4] += (gradN[i]*wc[0])*(tmp*dt);
				
				// calculate the kpc matrix for the anion
				wc[1] = (dKedc[1]*gp)*(-H[j])
				-Ke*((D[1]*gradN[j])*(kappa[1]/D0[1])
					 +((D[0]*(dkdc[0][1]-kappa[0]/D0[0]*dD0dc[0][1])
						+dDdc[0][1]*kappa[0])*gradc[0]
					   +(D[1]*(dkdc[1][1]-kappa[1]/D0[1]*dD0dc[1][1])
						 +dDdc[1][1]*kappa[1])*gradc[1]
					   )*H[j]
					 )*(R*T);
				ke[ndpn*i+3][ndpn*j+5] += (gradN[i]*wc[1])*(tmp*dt);
				
				// calculate the kcc matrix for the cation-cation
				jc[0][0] = (D[0]*(gradN[j]*(-phiw)+w*(H[j]/D0[0])))*kappa[0]
				+((D[0]*dkdc[0][0]+dDdc[0][0]*kappa[0])*(gradc[0]*(-phiw)+w*(c[0]/D0[0])))*H[j]
				+(D[0]*w*(-H[j]*dD0dc[0][0]/D0[0])+D[0]*wc[0])*(kappa[0]*c[0]/D0[0]);
				ke[ndpn*i+4][ndpn*j+4] += (gradN[i]*jc[0][0])*(tmp*dt);
				
				// calculate the kcc matrix for the anion-anion
				jc[1][1] = (D[1]*(gradN[j]*(-phiw)+w*(H[j]/D0[1])))*kappa[1]
				+((D[1]*dkdc[1][1]+dDdc[1][1]*kappa[1])*(gradc[1]*(-phiw)+w*(c[1]/D0[1])))*H[j]
				+(D[1]*w*(-H[j]*dD0dc[1][1]/D0[1])+D[1]*wc[1])*(kappa[1]*c[1]/D0[1]);
				ke[ndpn*i+5][ndpn*j+5] += (gradN[i]*jc[1][1])*(tmp*dt);
				
				// calculate the kcc matrix for the cation-anion
				jc[0][1] = 
				((D[0]*dkdc[0][1]+dDdc[0][1]*kappa[0])*(gradc[0]*(-phiw)+w*(c[0]/D0[0])))*H[j]
				+(D[0]*w*(-H[j]*dD0dc[0][1]/D0[0])+D[0]*wc[1])*(kappa[0]*c[0]/D0[0]);
				ke[ndpn*i+4][ndpn*j+5] += (gradN[i]*jc[0][1])*(tmp*dt);
				
				// calculate the kcc matrix for the anion-cation
				jc[1][0] = 
				((D[1]*dkdc[1][0]+dDdc[1][0]*kappa[1])*(gradc[1]*(-phiw)+w*(c[1]/D0[1])))*H[j]
				+(D[1]*w*(-H[j]*dD0dc[1][0]/D0[1])+D[1]*wc[0])*(kappa[1]*c[1]/D0[1]);
				ke[ndpn*i+5][ndpn*j+4] += (gradN[i]*jc[1][0])*(tmp*dt);
				
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

void FETriphasicDomain::SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	TriphasicMaterialStiffness(fem, el, ke);
	
	// calculate geometrical stiffness
	GeometricalStiffness(el, ke);
	
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

void FETriphasicDomain::TriphasicMaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	assert(fem.GetCurrentStep()->m_nModule == FE_TRIPHASIC);
	
	int i, i3, j, j3, n;
	
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	
	// global derivatives of shape functions
	// NOTE: hard-coding of hex elements!
	// Gx = dH/dx
	double Gx[8], Gy[8], Gz[8];
	
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
	
	// nodal concentrations
	double ct[2][8];
	for (i=0; i<neln; ++i) {
		ct[0][i] = m_pMesh->Node(el.m_node[i]).m_ct[0];
		ct[1][i] = m_pMesh->Node(el.m_node[i]).m_ct[1];
	}
	
	// weights at gauss points
	const double *gw = el.GaussWeights();
	
	// see if this is a biphasic-solute material
	FETriphasic* pmat = dynamic_cast<FETriphasic*>(fem.GetMaterial(el.GetMatID()));
	assert(pmat);
	
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
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		
		// evaluate concentration at gauss-point
		FESaltMaterialPoint& spt = *(mp.ExtractData<FESaltMaterialPoint>());
		spt.m_c[0] = el.Evaluate(ct[0], n);
		spt.m_c[1] = el.Evaluate(ct[1], n);
		
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
void FETriphasicDomain::UpdateStresses(FEModel &fem)
{
	int i, n;
	int nint, neln;
	double* gw;
	vec3d r0[8];
	vec3d rt[8];
	double pn[8], ct[2][8];
	
	FEMesh& mesh = *m_pMesh;
	FEM& mdl = dynamic_cast<FEM&>(fem);
	double dt = mdl.GetCurrentStep()->m_dt;
	bool sstate = (mdl.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE);
	
	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		assert(!el.IsRigid());
		
		assert(el.Type() != FE_UDGHEX);
		
		// get the number of integration points
		nint = el.GaussPoints();
		
		// get the number of nodes
		neln = el.Nodes();
		assert(neln <= 8);
		
		// get the integration weights
		gw = el.GaussWeights();
		
		// get the nodal data
		for (int j=0; j<neln; ++j)
		{
			r0[j] = mesh.Node(el.m_node[j]).m_r0;
			rt[j] = mesh.Node(el.m_node[j]).m_rt;
			pn[j] = mesh.Node(el.m_node[j]).m_pt;
			ct[0][j] = mesh.Node(el.m_node[j]).m_ct[0];
			ct[1][j] = mesh.Node(el.m_node[j]).m_ct[1];
		}
		
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(fem.GetMaterial(el.GetMatID()));
		
		assert(dynamic_cast<FETriphasic*>(pm) != 0);
		
		// get the biphasic-solute material
		FETriphasic* pmb = dynamic_cast<FETriphasic*>(pm);
		
		// extract the elastic component
		FEElasticMaterial* pme = pmb->m_pSolid;
		
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
			FESaltMaterialPoint& spt = *(mp.ExtractData<FESaltMaterialPoint>());
			
			// evaluate fluid pressure at gauss-point
			ppt.m_p = el.Evaluate(pn, n);
			
			// calculate the gradient of p at gauss-point
			ppt.m_gradp = gradient(el, pn, n);
			
			// evaluate effective solute concentration at gauss-point
			spt.m_c[0] = el.Evaluate(ct[0], n);
			spt.m_c[1] = el.Evaluate(ct[1], n);
			
			// calculate the gradient of c at gauss-point
			spt.m_gradc[0] = gradient(el, ct[0], n);
			spt.m_gradc[1] = gradient(el, ct[1], n);
			
			// for biphasic-solute materials also update the porosity, fluid and solute fluxes
			// and evaluate the actual fluid pressure and solute concentration
			ppt.m_w = pmb->FluidFlux(mp);
			spt.m_psi = pmb->ElectricPotential(mp);
			spt.m_ca[0] = pmb->Concentration(mp,0);
			spt.m_ca[1] = pmb->Concentration(mp,1);
			ppt.m_pa = pmb->Pressure(mp);
			spt.m_j[0] = pmb->SoluteFlux(mp,0);
			spt.m_j[1] = pmb->SoluteFlux(mp,1);
			spt.m_cF = pmb->FixedChargeDensity(mp);
			spt.m_Ie = pmb->CurrentDensity(mp);
			
			pt.s = pmb->Stress(mp);
		}
	}
}
