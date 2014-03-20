#include "FETriphasicDomain.h"
#include "FECore/FEMaterial.h"
#include "FETriphasic.h"
#include "FECore/log.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FETriphasicDomain::FETriphasicDomain(FEMesh* pm, FEMaterial* pmat) : FEBiphasicSoluteDomain(pm, pmat)
{
	m_ntype = FE_TRIPHASIC_DOMAIN; 
}

//-----------------------------------------------------------------------------
bool FETriphasicDomain::Initialize(FEModel &mdl)
{
	// initialize base class
	FEElasticSolidDomain::Initialize(mdl);
    
	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());
    
	// get the triphasic material
	FETriphasic* pmb = dynamic_cast<FETriphasic*>(pm); assert(pmb);
	const int nsol = 2;
	const int nsbm = 0;
    
	const int NE = FEElement::MAX_NODES;
    double p0[NE];
    vector< vector<double> > c0(nsol, vector<double>(NE));
	FEMesh& m = *GetMesh();
    
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
		for (int i=0; i<neln; ++i)
		{
			p0[i] = m.Node(el.m_node[i]).m_p0;
            for (int isol=0; isol<nsol; ++isol)
                c0[isol][i] = m.Node(el.m_node[i]).m_c0[isol];
		}

		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
            FEElasticMaterialPoint& pm = *(mp.ExtractData<FEElasticMaterialPoint>());
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
            
            // initialize effective fluid pressure, its gradient, and fluid flux
            pt.m_p = el.Evaluate(p0, n);
            pt.m_gradp = gradient(el, p0, n);
            pt.m_w = pmb->FluidFlux(mp);

            
			// initialize multiphasic solutes
			ps.m_nsol = nsol;
			ps.m_nsbm = nsbm;
            
            // initialize effective solute concentrations
            for (int isol=0; isol<nsol; ++isol) {
                ps.m_c[isol] = el.Evaluate(c0[isol].data(), n);
                ps.m_gradc[isol] = gradient(el, c0[isol].data(), n);
            }
			
            ps.m_psi = pmb->ElectricPotential(mp);
            for (int isol=0; isol<nsol; ++isol) {
                ps.m_ca[isol] = pmb->Concentration(mp,isol);
                ps.m_j[isol] = pmb->SoluteFlux(mp,isol);
            }
            pt.m_pa = pmb->Pressure(mp);
            ps.m_cF = pmb->FixedChargeDensity(mp);
            ps.m_Ie = pmb->CurrentDensity(mp);
			
            pm.m_s = pmb->Stress(mp);

		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Reset()
{
	// reset base class
	FEElasticSolidDomain::Reset();
	
	// get the material
	FEMaterial* pm = dynamic_cast<FEMaterial*>(GetMaterial());
    
	// get the multiphasic material
	FETriphasic* pmb = dynamic_cast<FETriphasic*>(pm); assert(pmb);
	const int nsol = 2;
	const int nsbm = 0;
	
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
		}
	}
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

/*
//-----------------------------------------------------------------------------
void FETriphasicDomain::Residual(FESolver* psolver, vector<double>& R)
{
	int i, j;
	
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());
	double dt = fem.GetCurrentStep()->m_dt;
	
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
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::ElementInternalForce(el, fe);
			
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			FEMaterial* pm = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pm) != 0);
			
			// calculate fluid internal work
			ElementInternalFluidWorkSS(el, fe, dt);
			
			// add fluid work to global residual
			int neln = el.Nodes();
			int J;
			for (j=0; j<neln; ++j)
			{
				J = elm[3*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate cation internal work
			ElementInternalSoluteWorkSS(el, fe, dt, 0);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[11*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate anion internal work
			ElementInternalSoluteWorkSS(el, fe, dt, 1);
			
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
			
			// unpack the element
			UnpackLM(el, elm);
			
			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);
			
			// calculate internal force vector
			// (This function is inherited from FEElasticSolidDomain)
			FEElasticSolidDomain::ElementInternalForce(el, fe);
			
			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, elm, fe, R);
			
			FEMaterial* pm = fem.GetMaterial(el.GetMatID());
			assert(dynamic_cast<FETriphasic*>(pm) != 0);
			
			// calculate fluid internal work
			ElementInternalFluidWork(el, fe, dt);
			
			// add fluid work to global residual
			int neln = el.Nodes();
			int J;
			for (j=0; j<neln; ++j)
			{
				J = elm[3*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate cation internal work
			ElementInternalSoluteWork(el, fe, dt, 0);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[11*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
			
			// calculate anion internal work
			ElementInternalSoluteWork(el, fe, dt, 1);
			
			// add solute work to global residual
			for (j=0; j<neln; ++j)
			{
				J = elm[12*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}
	}
}
*/

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalSoluteWorkSS(vector<double>& R, double dt)
{
	FETriphasic* pm = dynamic_cast<FETriphasic*>(GetMaterial()); assert(pm);
	
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE, pm)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> elm;
		
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
			
		// calculate cation internal work
		ElementInternalSoluteWorkSS(el, fe, dt, 0);
			
		// add solute work to global residual
		int dofc = DOF_C + pm->m_pSolute[0]->GetSoluteID();
		for (int j=0; j<neln; ++j)
		{
			int J = elm[dofc*neln+j];
			if (J >= 0) R[J] += fe[j];
		}

		// calculate anion internal work
		ElementInternalSoluteWorkSS(el, fe, dt, 1);
			
		// add solute work to global residual
        #pragma omp critical
        {
            dofc = DOF_C + pm->m_pSolute[1]->GetSoluteID();
            for (int j=0; j<neln; ++j)
            {
                int J = elm[dofc*neln+j];
                if (J >= 0) R[J] += fe[j];
            }
        }
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalSoluteWork(vector<double>& R, double dt)
{
	// element force vector
	vector<double> fe;
	vector<int> elm;

	FETriphasic* pm = dynamic_cast<FETriphasic*>(GetMaterial()); assert(pm);
	
	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// unpack the element
		UnpackLM(el, elm);
			
		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);
			
		// calculate cation internal work
		ElementInternalSoluteWork(el, fe, dt, 0);
			
		// add solute work to global residual
		int dofc = DOF_C + pm->m_pSolute[0]->GetSoluteID();
		for (int j=0; j<neln; ++j)
		{
			int J = elm[dofc*neln+j];
			if (J >= 0) R[J] += fe[j];
		}

		// calculate anion internal work
		ElementInternalSoluteWork(el, fe, dt, 1);
			
		// add solute work to global residual
		dofc = DOF_C + pm->m_pSolute[1]->GetSoluteID();
		for (int j=0; j<neln; ++j)
		{
			int J = elm[dofc*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalFluidWorkSS(vector<double>& R, double dt)
{
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> elm;
		
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
        #pragma omp critical
        {
            int neln = el.Nodes();
            for (int j=0; j<neln; ++j)
            {
                int J = elm[3*neln+j];
                if (J >= 0) R[J] += fe[j];
            }
        }
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalFluidWork(vector<double>& R, double dt)
{
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> elm;
		
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
        #pragma omp critical
        {
            int neln = el.Nodes();
            for (int j=0; j<neln; ++j)
            {
                int J = elm[3*neln+j];
                if (J >= 0) R[J] += fe[j];
            }
        }
	}
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FETriphasicDomain::ElementInternalFluidWork(FESolidElement& el, vector<double>& fe, double dt)
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

bool FETriphasicDomain::ElementInternalFluidWorkSS(FESolidElement& el, vector<double>& fe, double dt)
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
	
	zero(fe);
	
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

bool FETriphasicDomain::ElementInternalSoluteWork(FESolidElement& el, vector<double>& fe, double dt, const int ion)
{
	int i, n;
	
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
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial());
	assert(pm);
	
	int id0 = pm->m_pSolute[0]->GetSoluteID();
	int id1 = pm->m_pSolute[1]->GetSoluteID();
	
	const int NE = FEElement::MAX_NODES;
	vec3d r0[NE], rt[NE], rp[NE], vt[NE];
	double cp[2][NE];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		cp[0][i] = mesh.Node(el.m_node[i]).m_cp[id0];
		cp[1][i] = mesh.Node(el.m_node[i]).m_cp[id1];
		vt[i] = mesh.Node(el.m_node[i]).m_vt;
	}

	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(el.m_State[n]->ExtractData<FESolutesMaterialPoint>());
		
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
			
			// calculate solid velocity
			vs += vt[i]*H[i];
			
			// save spatial gradient of shape functions
			gradN[i] = vec3d(Gx,Gy,Gz);
			
			// calculate effective concentration at previous time step
			cprev[0] += cp[0][i]*H[i];
			cprev[1] += cp[1][i]*H[i];
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.m_J;
		double dJdt = (J-Jp)/dt;

		// and then finally
		double divv = dJdt/J;
		
		// get the solute flux
		vec3d j[2] = {spt.m_j[0],spt.m_j[1]};
		// get the effective concentration
		double c[2] = {spt.m_c[0],spt.m_c[1]};

		// get the charge number
		int z[2] = {pm->m_pSolute[0]->ChargeNumber(),
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
		// evaluate time derivatives of concentration, solubility and porosity
		double dcdt[2] = {(c[0] - cprev[0])/dt, (c[1] - cprev[1])/dt};
		double dkdt[2] = {dkdJ[0]*dJdt + dkdc[0][0]*dcdt[0] + dkdc[0][1]*dcdt[1],
			dkdJ[1]*dJdt + dkdc[1][0]*dcdt[0] + dkdc[1][1]*dcdt[1]};
		double dpdt = dpdJ*dJdt;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[ion]+(j[0]*z[0]+j[1]*z[1])*pm->m_penalty)
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

bool FETriphasicDomain::ElementInternalSoluteWorkSS(FESolidElement& el, vector<double>& fe, double dt, const int ion)
{
	int i, n;
	
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
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial());
	assert(pm);
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
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
		
		// get the solute flux
		vec3d j[2] = {spt.m_j[0],spt.m_j[1]};

		// get the charge number
		int z[2] = {pm->m_pSolute[0]->ChargeNumber(),
			pm->m_pSolute[1]->ChargeNumber()};

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*(j[ion]+(j[0]*z[0]+j[1]*z[1])*pm->m_penalty)
						 )*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
/*
void FETriphasicDomain::StiffnessMatrix(FESolver* psolver)
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
*/

//-----------------------------------------------------------------------------
void FETriphasicDomain::StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp)
{
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial()); assert(pm);
	int dofc0 = DOF_C + pm->m_pSolute[0]->GetSoluteID();
	int dofc1 = DOF_C + pm->m_pSolute[1]->GetSoluteID();
	
	// repeat over all solid elements
	int NE = m_Elem.size();
    
	#pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> elm;
		
		FESolidElement& el = m_Elem[iel];
		UnpackLM(el, elm);
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 6;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementTriphasicStiffness(el, ke, bsymm, tp.dt);
		
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
			lm[ndpn*i+4] = elm[dofc0*neln+i];
			lm[ndpn*i+5] = elm[dofc1*neln+i];
		}
		
		// assemble element matrix in global stiffness matrix
        #pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::StiffnessMatrixSS(FESolver* psolver, bool bsymm, const FETimePoint& tp)
{
	// get the elements material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial()); assert(pm);
	int dofc0 = DOF_C + pm->m_pSolute[0]->GetSoluteID();
	int dofc1 = DOF_C + pm->m_pSolute[1]->GetSoluteID();
	
	// repeat over all solid elements
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> elm;
		
		FESolidElement& el = m_Elem[iel];
		UnpackLM(el, elm);
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 6;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementTriphasicStiffnessSS(el, ke, bsymm, tp.dt);
		
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
			lm[ndpn*i+4] = elm[dofc0*neln+i];
			lm[ndpn*i+5] = elm[dofc1*neln+i];
		}
		
		// assemble element matrix in global stiffness matrix
        #pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FETriphasicDomain::ElementTriphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt)
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
	
	FEMesh& mesh = *GetMesh();
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial());
	assert(pm);
    
	const int nsol = 2;
	int ndpn = 4+nsol;
    int sid[nsol];
	for (isol=0; isol<nsol; ++isol)
		sid[isol] = pm->m_pSolute[isol]->GetSoluteID();
	
	vec3d rp[FEElement::MAX_NODES], v[FEElement::MAX_NODES];
    double cp[nsol][FEElement::MAX_NODES];
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
		
		// get the charge density and its derivatives
		double phi0 = ppt.m_phi0;
		double cF = pm->FixedChargeDensity(mp);
		double dcFdJ = -cF/(J - phi0);
		double dcFdJJ = 2*cF/SQR(J-phi0);
		
		double c[nsol];
		vec3d gradc[nsol];
		double dcdt[nsol];
		int z[nsol];
		double khat[nsol];
		double dkhdJ[nsol];
		double dkhdJJ[nsol];
		double dkhdc[nsol][nsol];
		double dkhdJc[nsol][nsol];
		double dkhdcc[nsol][nsol][nsol];
		
		double zeta = pm->ElectricPotential(mp, true);
		double zz[nsol];
		double kappa[nsol];
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
		double zidzdc[nsol] = {0};
		double zidzdJc[nsol] = {0}, zidzdJc1[nsol] = {0}, zidzdJc2[nsol] = {0};
		double zidzdcc[nsol][nsol] = {0};
		double zidzdcc1[nsol][nsol] = {0};
		double zidzdcc2[nsol] = {0};
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
		
		double dkdJ[nsol];
		double dkdc[nsol][nsol];
		double dkdJJ[nsol];
		double dkdJc[nsol][nsol];
		double dkdcc[nsol][nsol][nsol];
		
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
		
		mat3ds dKdc[nsol];
		mat3ds D[nsol];
		tens4ds dDdE[nsol];
		mat3ds dDdc[nsol][nsol];
		vector<double> D0(nsol);
		double dD0dc[nsol][nsol];
		double dodc[nsol];
		mat3ds dTdc[nsol];
		mat3ds ImD[nsol];
		mat3dd I(1);
		
		// evaluate the solvent supply and its derivatives
		double phiwhat = 0;
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		double Phic[nsol] = {0};
		
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
//			dTdc[isol] = pm->m_pSolid->Tangent_Concentration(mp,isol);
			dTdc[isol] = mat3ds(0,0,0,0,0,0);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
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
		mat3ds Gc[nsol];
		mat3ds dKedc[nsol];
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
		vec3d gc[nsol],qcu[nsol],wc[nsol],jce[nsol];
		vec3d jc[nsol][nsol];
		mat3d wu, jue;
		mat3d ju[nsol];
		double qcc[nsol][nsol];
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
				qpu = -gradN[j]*(divv+1.0/dt)+(Phie + gradv.transpose())*gradN[j];
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
					ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol]+H[i]*H[j]*Phic[isol])*(tmp*dt);
					
				}
				
				// calculate data for the kcc matrix
                jce[0] = jce[1] = vec3d(0,0,0);
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
bool FETriphasicDomain::ElementTriphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt)
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
	
	FEMesh& mesh = *GetMesh();
	
	// get the element's material
	FETriphasic* pm = dynamic_cast<FETriphasic*> (GetMaterial());
	assert(pm);
    
	const int nsol = 2;
	int ndpn = 4+nsol;
	int sid[nsol];
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
		double J = ept.m_J;
		
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
		double zz[nsol];
		double kappa[nsol];
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
		double zidzdc[nsol] = {0};
		double zidzdJc[nsol] = {0}, zidzdJc1[nsol] = {0}, zidzdJc2[nsol] = {0};
		double zidzdcc[nsol][nsol] = {0}, zidzdcc1[nsol][nsol] = {0};
		double zidzdcc2[nsol] = {0};
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
		
		double dkdJ[nsol];
		double dkdc[nsol][nsol];
		double dkdJJ[nsol];
		double dkdJc[nsol][nsol];
		double dkdcc[nsol][nsol][nsol];
		
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
		
		mat3ds dKdc[nsol];
		mat3ds D[nsol];
		tens4ds dDdE[nsol];
		mat3ds dDdc[nsol][nsol];
		double D0[nsol];
		double dD0dc[nsol][nsol];
		double dodc[nsol];
		mat3ds dTdc[nsol];
		mat3ds ImD[nsol];
		mat3dd I(1);
		
		// evaluate the solvent supply and its derivatives
		double phiwhat = 0;
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		double Phic[nsol] = {0};
		
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
//			dTdc[isol] = pm->m_pSolid->Tangent_Concentration(mp,isol);
			dTdc[isol] = mat3ds(0,0,0,0,0,0);
			
			ImD[isol] = I-D[isol]/D0[isol];
			
			for (jsol=0; jsol<nsol; ++jsol) {
				dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
				dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
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
		mat3ds Gc[nsol];
		mat3ds dKedc[nsol];
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
		vec3d gc[nsol],wc[nsol],jce[nsol];
		vec3d jc[nsol][nsol];
		mat3d wu, jue;
		mat3d ju[nsol];
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
				jce[0] = jce[1] = vec3d(0,0,0);
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

void FETriphasicDomain::SolidElementStiffness(FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	ElementTriphasicMaterialStiffness(el, ke);
	
	// calculate geometrical stiffness (inherited from FEElasticSolidDomain)
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

void FETriphasicDomain::ElementTriphasicMaterialStiffness(FESolidElement &el, matrix &ke)
{
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

	// see if this is a biphasic-solute material
	FETriphasic* pmat = dynamic_cast<FETriphasic*>(GetMaterial());	assert(pmat);
	int id0 = pmat->m_pSolute[0]->GetSoluteID();
	int id1 = pmat->m_pSolute[1]->GetSoluteID();
	
	// nodal concentrations
	double ct[2][FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) {
		ct[0][i] = m_pMesh->Node(el.m_node[i]).m_ct[id0];
		ct[1][i] = m_pMesh->Node(el.m_node[i]).m_ct[id1];
	}
	
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
		
		// evaluate concentration at gauss-point
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
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
	bool berr = false;
	int NE = (int) m_Elem.size();
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			UpdateElementStress(i);
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
void FETriphasicDomain::UpdateElementStress(int iel)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	// get the number of nodes
	int neln = el.Nodes();
		
	// get the integration weights
	double* gw = el.GaussWeights();

	// get the biphasic-solute material
	FETriphasic* pmb = dynamic_cast<FETriphasic*>(GetMaterial()); assert(pmb);
	int id0 = pmb->m_pSolute[0]->GetSoluteID();
	int id1 = pmb->m_pSolute[1]->GetSoluteID();
		
	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES], ct[2][FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
		pn[j] = mesh.Node(el.m_node[j]).m_pt;
		ct[0][j] = mesh.Node(el.m_node[j]).m_ct[id0];
		ct[1][j] = mesh.Node(el.m_node[j]).m_ct[id1];
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
			
		pt.m_s = pmb->Stress(mp);
	}
}
