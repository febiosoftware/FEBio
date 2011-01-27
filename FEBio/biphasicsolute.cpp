#include "stdafx.h"
#include "FEDomain.h"
#include "FESolidSolver.h"
#include "FEMaterial.h"
#include "FEMicroMaterial.h"
#include "log.h"

//-----------------------------------------------------------------------------
void FEBiphasicSoluteDomain::Residual(FESolidSolver* psolver, vector<double>& R)
{
	int i, j;
	
	FEM& fem = psolver->m_fem;
	
	// make sure we are in poro-solute mode
	assert(fem.m_pStep->m_nModule == FE_POROSOLUTE);
	
	// element force vector
	vector<double> fe;
	
	int NE = m_Elem.size();
	for (i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		
		// this element should not be rigid
		assert(!el.IsRigid());
		
		//! this element should not be UDG
		assert(el.Type() != FE_UDGHEX);
		
		// unpack the element
		UnpackElement(el);
		
		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
		
		// calculate internal force vector
		// (This function is inherited from FEElasticSolidDomain)
		InternalForces(el, fe);
		
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
		psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
		
		FEMaterial* pm = fem.GetMaterial(el.GetMatID());
		assert(dynamic_cast<FEBiphasicSolute*>(pm) != 0);
		
		// calculate fluid internal work
		InternalFluidWork(fem, el, fe);
		
		// add fluid work to global residual
		int neln = el.Nodes();
		vector<int>& lm = el.LM();
		int J;
		for (j=0; j<neln; ++j)
		{
			J = lm[3*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
		
		// calculate solute internal work
		InternalSoluteWork(fem, el, fe);
		
		// add solute work to global residual
		for (j=0; j<neln; ++j)
		{
			J = lm[11*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FEBiphasicSoluteDomain::InternalFluidWork(FEM& fem, FESolidElement& el, vector<double>& fe)
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
	
	// get the logfile
	Logfile& log = GetLogfile();
	
	// get the element's material
	FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.m_pStep->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutePoroElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FESolutePoroElasticMaterialPoint>());
		
		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		el.invjac0(J0i, n);
		
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
		vec3d& w = pt.m_w;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z - divv*H[i])*detJ*wg[n];
			//			fe[i] -= (B1[i]*w.x+B2[i]*w.y+B3[i]*w.z - divv*H[i])*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FEBiphasicSoluteDomain::InternalSoluteWork(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int i, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	double *Gr, *Gs, *Gt, *H;
	double *Grr, *Gsr, *Gtr, *Grs, *Gss, *Gts, *Grt, *Gst, *Gtt;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	vec3d* vt = el.vt();
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d r0[8], rt[8], rp[8];
	double cp[8];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[i] = mesh.Node(el.m_node[i]).m_cp;
	}
	
	// get the logfile
	Logfile& log = GetLogfile();
	
	// get the element's material
	FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	zero(fe);
	
	// get the time step value
	double dt = fem.m_pStep->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutePoroElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FESolutePoroElasticMaterialPoint>());
		
		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);
		vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
		vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
		vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		el.invjac0(J0i, n);
		vec3d G1(J0i[0][0],J0i[0][1],J0i[0][2]);
		vec3d G2(J0i[1][0],J0i[1][1],J0i[1][2]);
		vec3d G3(J0i[2][0],J0i[2][1],J0i[2][2]);
		
		// next we calculate the deformation gradient and the solid velocity
		mat3d Fp;
		Fp.zero();
		vec3d vs(0);
		vec3d gradJ(0);
		double cprev = 0;
		
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
			cprev += cp[i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.J;
		double dJdt = (J-Jp)/dt;
		gradJ *= J;
		
		// and then finally
		double divv = dJdt/J;
		
		// get the solute flux
		vec3d& j = pt.m_j;
		// get the effective concentration and its gradient
		double c = pt.m_c;
		vec3d gradc = pt.m_gradc;
		
		// evaluate the solubility and its derivatives w.r.t. J and c, and its gradient
		double kappa = pm->m_pSolub->Solubility(mp);
		double dkdJ = pm->m_pSolub->Tangent_Solubility_Strain(mp);
		double dkdc = pm->m_pSolub->Tangent_Solubility_Concentration(mp);
		vec3d gradk = gradJ*dkdJ;
		// evaluate the porosity, its derivative w.r.t. J, and its gradient
		double phiw = pm->Porosity(mp);
		double dpdJ = (1. - phiw)/J;
		vec3d gradp = gradJ*dpdJ;
		// evaluate time derivatives of concentration, solubility and porosity
		double dcdt = (c - cprev)/dt;
		double dkdt = dkdJ*dJdt + dkdc*dcdt;
		double dpdt = dpdJ*dJdt;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*j.x+B2[i]*j.y+B3[i]*j.z 
						 - H[i]*(dpdt*kappa*c+phiw*dkdt*c+phiw*kappa*dcdt
								 +vs*(gradp*kappa*c+gradk*phiw*c+gradc*phiw*kappa)
								 +phiw*kappa*c*divv)
						 )*detJ*wg[n];
		}
	}
	
	return true;
}


//-----------------------------------------------------------------------------

void FEBiphasicSoluteDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;
	
	// element stiffness matrix
	matrix ke;
	
	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		
		// this element should not be rigid
		assert(!el.IsRigid());
		
		UnpackElement(el);
		
		// get the elements material
		FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
		assert(dynamic_cast<FEBiphasicSolute*>(pmat) != 0);
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*5;
		ke.Create(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicSoluteStiffness(fem, el, ke);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackElement
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackElement function so that it returns
		// the LM vector in the right order for solute-solid elements.
		vector<int> lm(ndof);
		for (int i=0; i<neln; ++i)
		{
			lm[5*i  ] = el.LM()[3*i];
			lm[5*i+1] = el.LM()[3*i+1];
			lm[5*i+2] = el.LM()[3*i+2];
			lm[5*i+3] = el.LM()[3*neln+i];
			lm[5*i+4] = el.LM()[11*neln+i];
		}
		
		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEBiphasicSoluteDomain::ElementBiphasicSoluteStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;
	
	// jacobian
	double Ji[3][3], detJ, J0i[3][3];
	
	vec3d* vt = el.vt();
	
	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	vector<vec3d> gradN(neln);
	double tmp;
	
	// gauss-weights
	double* gw = el.GaussWeights();
	
	FEMesh& mesh = fem.m_mesh;
	
	vec3d r0[8], rt[8], rp[8], v[8];
	double cp[8];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[i] = mesh.Node(el.m_node[i]).m_cp;
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
			ke[5*i  ][5*j] = ks[3*i  ][3*j  ]; ke[5*i  ][5*j+1] = ks[3*i  ][3*j+1]; ke[5*i  ][5*j+2] = ks[3*i  ][3*j+2];
			ke[5*i+1][5*j] = ks[3*i+1][3*j  ]; ke[5*i+1][5*j+1] = ks[3*i+1][3*j+1]; ke[5*i+1][5*j+2] = ks[3*i+1][3*j+2];
			ke[5*i+2][5*j] = ks[3*i+2][3*j  ]; ke[5*i+2][5*j+1] = ks[3*i+2][3*j+1]; ke[5*i+2][5*j+2] = ks[3*i+2][3*j+2];
		}
	
	// get the logfile
	Logfile& log = GetLogfile();
	
	// get the element's material
	FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}
	
	// check if we use the symmetric version of the poro-implementation
	bool bsymm = fem.m_bsym_poro;
	double dt = fem.m_pStep->m_dt;
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutePoroElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FESolutePoroElasticMaterialPoint>());
		
		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);
		vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
		vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
		vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
		
		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		el.invjac0(J0i, n);
		vec3d G1(J0i[0][0],J0i[0][1],J0i[0][2]);
		vec3d G2(J0i[1][0],J0i[1][1],J0i[1][2]);
		vec3d G3(J0i[2][0],J0i[2][1],J0i[2][2]);
		
		// next we calculate the deformation gradient and the solid velocity
		mat3d Fp, gradv;
		Fp.zero();
		gradv.zero();
		vec3d vs(0);
		double cprev = 0;
		
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
			vs += vt[i]*H[i];
			gradv[0][0] += v[i].x*Gx; gradv[1][0] += v[i].y*Gx; gradv[2][0] += v[i].z*Gx;
			gradv[0][1] += v[i].x*Gy; gradv[1][1] += v[i].y*Gy; gradv[2][1] += v[i].z*Gy;
			gradv[0][2] += v[i].x*Gz; gradv[1][2] += v[i].y*Gz; gradv[2][2] += v[i].z*Gz;
			
			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
			gradN[i] = vec3d(Gx,Gy,Gz);
			
			// calculate effective concentration at previous time step
			cprev += cp[i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.J;
		
		// and then finally
		double dJdt = (J-Jp)/dt;
		double divv = dJdt/J;
		
		// get the fluid flux and pressure gradient
		vec3d w = pt.m_w;
		vec3d gradp = pt.m_gradp;
		
		// get the effective concentration, its gradient and its time derivative
		double c = pt.m_c;
		vec3d gradc = pt.m_gradc;
		double dcdt = (c - cprev)/dt;
		
		// evaluate the permeability and its derivatives
		mat3ds K = pm->m_pPerm->Permeability(mp);
		tens4ds dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
		mat3ds dKdc = pm->m_pPerm->Tangent_Permeability_Concentration(mp); 
		
		// evaluate the porosity and its derivative
		double phiw = pm->Porosity(mp);
		double dpdJ = (1. - phiw)/J;
		
		// evaluate the solubility and its derivatives
		double kappa = pm->m_pSolub->Solubility(mp);
		double dkdJ = pm->m_pSolub->Tangent_Solubility_Strain(mp);
		double dkdJJ = pm->m_pSolub->Tangent_Solubility_Strain_Strain(mp);
		double dkdc = pm->m_pSolub->Tangent_Solubility_Concentration(mp);
		double dkdJc = pm->m_pSolub->Tangent_Solubility_Strain_Concentration(mp);
		
		// evaluate the diffusivity tensor and its derivatives
		mat3ds D = pm->m_pDiff->Diffusivity(mp);
		mat3ds dDdc = pm->m_pDiff->Tangent_Diffusivity_Concentration(mp);
		tens4ds dDdE = pm->m_pDiff->Tangent_Diffusivity_Strain(mp);
		
		// evaluate the solute free diffusivity
		double D0 = pm->m_pDiff->Free_Diffusivity(mp);
		double dD0dc = 0;	// temporary, until I implement it in FESoluteDiffusivity
		
		// evaluate the osmotic coefficient and its derivatives
		double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
		double dodc = pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp);
		
		// Miscellaneous constants
		mat3dd I(1);
		double R = pm->m_Rgas;
		double T = pm->m_Tabs;
		
		// evaluate the effective permeability and its derivatives
		mat3ds Ki = K.inverse();
		mat3ds ImD = I-D/D0;
		mat3ds Ke = (Ki + ImD*(R*T*kappa*c/phiw/D0)).inverse();
		tens4ds G = dyad1s(Ki,I) - dyad4s(Ki,I)*2 - ddots(dyad2s(Ki),dKdE)*0.5
		+dyad1s(ImD,I)*(R*T*c*J/D0/2/phiw*(dkdJ-kappa/phiw*dpdJ))
		+(dyad1s(I) - dyad4s(I)*2 - dDdE/D0)*(R*T*kappa*c/phiw/D0);
		tens4ds dKedE = dyad1s(Ke,I) - 2*dyad4s(Ke,I) - ddots(dyad2s(Ke),G)*0.5;
		mat3ds Gc = -Ki*dKdc*Ki + ImD*(R*T/phiw/D0*(dkdc*c+kappa-kappa*c/D0*dD0dc))
		+R*T*kappa*c/phiw/D0/D0*(D*dD0dc/D0 - dDdc);
		mat3ds dKedc = -Ke*Gc*Ke;
		
		// calculate the kpp matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				ke[5*i+3][5*j+3] -= gradN[i]*(Ke*gradN[j])*tmp;
			}
		
		// calculate the kcc matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				ke[5*i+4][5*j+4] += (H[j]*(gradN[i]*((D*dkdc+dDdc*kappa)*(w*c/D0-gradc*phiw)))
									 -phiw*kappa*(gradN[i]*(D*gradN[j]))
									 -H[j]*kappa*c/D0*(gradN[i]*((D*dKdc)*(gradp+D*gradc*(R*T*kappa/D0))))
									 -R*T*kappa*c/D0*(gradN[i]*((D*Ke)*(H[j]*(D*(dkdc-kappa/D0*dD0dc)/D0+dDdc*kappa/D0)*gradc+D*gradN[j]*kappa/D0)))
									 +H[j]*phiw*(kappa+c*dkdc)*(vs*gradN[i])
									 -H[i]*H[j]*divv*(phiw*kappa+J*dpdJ*kappa+J*phiw*dkdJ+c*(phiw*dkdc+J*dpdJ*dkdc+J*phiw*dkdJc))
									 -phiw*H[i]*H[j]*(dkdc*dcdt+kappa/dt))*tmp;
			}
		
		// calculate the kup matrix
		for (i=0; i<neln; ++i) {
			for (j=0; j<neln; ++j)
			{
				tmp = detJ*gw[n]*H[j];
				ke[5*i  ][5*j+3] -= tmp*gradN[i].x;
				ke[5*i+1][5*j+3] -= tmp*gradN[i].y;
				ke[5*i+2][5*j+3] -= tmp*gradN[i].z;
			}
		}
		
		// calculate the kpu matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i) {
			for (j=0; j<neln; ++j)
			{
				vec3d vt = (-(vdotTdotv(gradN[i], dKedE, gradN[j])*(gradp+(D*gradc)*(R*T*kappa/D0)))
							-gradN[j]*(gradN[i]*((Ke*D)*gradc))*(R*T/D0*(J*dkdJ-kappa))
							-(Ke*gradN[i])*(gradN[j]*(D*gradc))*(2*R*T*kappa/D0)
							-vdotTdotv(gradc, dDdE, gradN[j])*(Ke*gradN[i])*(R*T*kappa/D0)
							-(I*(divv+1./dt) - gradv.transpose())*gradN[j]*H[i])*tmp;
				ke[5*i+3][5*j  ] += vt.x;
				ke[5*i+3][5*j+1] += vt.y;
				ke[5*i+3][5*j+2] += vt.z;
			}
		}
		
		// calculate the kuc matrix
		mat3ds dTdc = pm->m_pSolid->Tangent_Concentration(mp);
		for (i=0; i<neln; ++i) {
			vec3d vt = dTdc*gradN[i]-gradN[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa));
			for (j=0; j<neln; ++j)
			{
				tmp = detJ*gw[n]*H[j];
				ke[5*i  ][5*j+4] += tmp*vt.x;
				ke[5*i+1][5*j+4] += tmp*vt.y;
				ke[5*i+2][5*j+4] += tmp*vt.z;
			}
		}
		
		// calculate the kcu matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i) {
			for (j=0; j<neln; ++j)
			{
				vec3d vt = (D*gradN[i])*(gradc*gradN[j]) 
				+ ((vdotTdotv(gradN[i], dDdE, gradN[j])).transpose()*gradc)*kappa;
				vt += (gradN[j]*(kappa+J*dkdJ)*(vs*gradN[i])
					   +gradN[i]*kappa*(H[j]/dt-(vs*gradN[j])))*c;
				double s = (kappa+3*J*dkdJ+J*J*dkdJJ)*c*divv+(kappa+J*dkdJ)*(dcdt+c/dt);
				ke[5*i+4][5*j  ] += tmp*(vt.x-H[i]*s*gradN[j].x);
				ke[5*i+4][5*j+1] += tmp*(vt.y-H[i]*s*gradN[j].y);
				ke[5*i+4][5*j+2] += tmp*(vt.z-H[i]*s*gradN[j].z);
			}
		}
		
		// calculate the kcp matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i) {
			for (j=0; j<neln; ++j)
			{
				ke[5*i+4][5*j+3] -= (gradN[i]*((D*Ke)*gradN[j]))*(kappa*c/D0)*tmp;
			}
		}
		
		// calculate the kpc matrix
		tmp = dt*detJ*gw[n];
		for (i=0; i<neln; ++i) {
			for (j=0; j<neln; ++j)
			{
				ke[5*i+3][5*j+4] -= (gradN[i]*(dKedc*(gradp+(D*gradc)*(R*T*kappa/D0)))*H[j]
									 +(gradN[i]*((Ke*(D*(dkdc-kappa*dD0dc/D0)/D0+dDdc*(kappa/D0)))*gradc))*(R*T*H[j])
									 +(gradN[i]*((Ke*D)*gradN[j]))*R*T*kappa/D0)*tmp;
			}
		}
		
	}
	
	// Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
	if (bsymm) {
		for (i=0; i<5*neln; ++i)
			for (j=i+1; j<5*neln; ++j) {
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

void FEBiphasicSoluteDomain::SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	BiphasicSoluteMaterialStiffness(fem, el, ke);
	
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

void FEBiphasicSoluteDomain::BiphasicSoluteMaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	assert(fem.m_pStep->m_nModule == FE_POROSOLUTE);
	
	int i, i3, j, j3, n;
	
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;
	
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
	
	// weights at gauss points
	const double *gw = el.GaussWeights();
	
	// see if this is a biphasic-solute material
	FEBiphasicSolute* pmat = dynamic_cast<FEBiphasicSolute*>(fem.GetMaterial(el.GetMatID()));
	assert(pmat);
	
	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];
		
		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);
		
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		
		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		//		el.defgrad(pt.F, n);
		//		pt.J = el.detF(n);
		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;
		
		// evaluate concentration at gauss-point
		FESolutePoroElasticMaterialPoint& ppt = *(mp.ExtractData<FESolutePoroElasticMaterialPoint>());
		ppt.m_c = el.Evaluate(el.ct(), n);
		
		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
		C.extract(D);
		
		if (dynamic_cast<FEMicroMaterial*>(pmat))
		{
			// the micro-material screws up the currently unpacked elements
			// so I have to unpack the element data again
			UnpackElement(el);
		}
		
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
void FEBiphasicSoluteDomain::UpdateStresses(FEM &fem)
{
	int i, n;
	int nint;
	double* gw;
	
	assert(fem.m_pStep->m_nModule == FE_POROSOLUTE);
	
	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		assert(!el.IsRigid());
		
		assert(el.Type() != FE_UDGHEX);
		
		// unpack the element data
		UnpackElement(el);
		
		// get the number of integration points
		nint = el.GaussPoints();
		
		// get the integration weights
		gw = el.GaussWeights();
		
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(fem.GetMaterial(el.GetMatID()));
		
		assert(dynamic_cast<FEBiphasicSolute*>(pm) != 0);
		
		// get the biphasic-solute material
		FEBiphasicSolute* pmb = dynamic_cast<FEBiphasicSolute*>(pm);
		
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
			pt.r0 = el.Evaluate(el.r0(), n);
			pt.rt = el.Evaluate(el.rt(), n);
			
			// get the deformation gradient and determinant
			pt.J = el.defgrad(pt.F, n);
			
			// three-field element variables
			pt.avgJ = el.m_eJ;
			pt.avgp = el.m_ep;
			
			// solute-poroelastic data
			FESolutePoroElasticMaterialPoint& ppt = *(mp.ExtractData<FESolutePoroElasticMaterialPoint>());
			
			// evaluate fluid pressure at gauss-point
			ppt.m_p = el.Evaluate(el.pt(), n);
			
			// calculate the gradient of p at gauss-point
			ppt.m_gradp = el.gradient(el.pt(), n);
			
			// evaluate effective solute concentration at gauss-point
			ppt.m_c = el.Evaluate(el.ct(), n);
			
			// calculate the gradient of c at gauss-point
			ppt.m_gradc = el.gradient(el.ct(), n);
			
			if (dynamic_cast<FEMicroMaterial*>(pme))
			{
				// the micro-material screws up the currently unpacked elements
				// so I have to unpack the element data again
				UnpackElement(el);
			}
			
			// for biphasic-solute materials also update the fluid and solute fluxes
			// and evaluate the actual fluid pressure and solute concentration
			ppt.m_w = pmb->FluidFlux(mp);
			ppt.m_pa = pmb->Pressure(mp);
			ppt.m_j = pmb->SoluteFlux(mp);
			ppt.m_ca = pmb->Concentration(mp);
			
			// calculate the stress at this material point (must be done after evaluating m_pa)
			pt.s = pmb->Stress(mp);
		}
	}
}
