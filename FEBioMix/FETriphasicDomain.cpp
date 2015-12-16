#include "FETriphasicDomain.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FETriphasicDomain::FETriphasicDomain(FEModel* pfem) : FESolidDomain(&pfem->GetMesh()), FEElasticDomain(pfem)
{
	m_pMat = 0;
	m_dofP = pfem->GetDOFIndex("p");
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FETriphasic*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FETriphasicDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int dofc0 = DOF_C + m_pMat->m_pSolute[0]->GetSoluteID();
	int dofc1 = DOF_C + m_pMat->m_pSolute[1]->GetSoluteID();

	int N = el.Nodes();
	lm.resize(N*9);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];

		// now the pressure dofs
		lm[3*N+i] = id[m_dofP];

		// rigid rotational dofs
		// TODO: Do I really need these?
		lm[4*N + 3*i  ] = id[DOF_RU];
		lm[4*N + 3*i+1] = id[DOF_RV];
		lm[4*N + 3*i+2] = id[DOF_RW];

		// concentration dofs
		lm[7*N + i] = id[dofc0];
		lm[8*N + i] = id[dofc1];
	}
}

//-----------------------------------------------------------------------------
bool FETriphasicDomain::Initialize(FEModel &fem)
{
	// initialize base class
	FESolidDomain::Initialize(fem);

	// initialize local coordinate systems (can I do this elsewhere?)
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
	}
    
	return true;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Activate()
{
	int dofc0 = DOF_C + m_pMat->m_pSolute[0]->GetSoluteID();
	int dofc1 = DOF_C + m_pMat->m_pSolute[1]->GetSoluteID();

	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[m_dofX] = DOF_ACTIVE;
				node.m_ID[m_dofY] = DOF_ACTIVE;
				node.m_ID[m_dofZ] = DOF_ACTIVE;
			}

			node.m_ID[m_dofP] = DOF_ACTIVE;
			node.m_ID[dofc0] = DOF_ACTIVE;
			node.m_ID[dofc1] = DOF_ACTIVE;
		}
	}

	// get the triphasic material
	FETriphasic* pmb = m_pMat;
	const int nsol = 2;
	const int nsbm = 0;

	const int NE = FEElement::MAX_NODES;
	double p0[NE];
	vector< vector<double> > c0(nsol, vector<double>(NE));
	FEMesh& m = *GetMesh();

	int id[2] = { m_pMat->m_pSolute[0]->GetSoluteID(), m_pMat->m_pSolute[1]->GetSoluteID() };

	for (int i = 0; i<(int)m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the number of nodes
		int neln = el.Nodes();
		// get initial values of fluid pressure and solute concentrations
		for (int i = 0; i<neln; ++i)
		{
//			p0[i] = m.Node(el.m_node[i]).m_p0;
			p0[i] = m.Node(el.m_node[i]).get(m_dofP);
			for (int isol = 0; isol<nsol; ++isol)
				//				c0[isol][i] = m.Node(el.m_node[i]).m_c0[id[isol]];
				c0[isol][i] = m.Node(el.m_node[i]).get(DOF_C + id[isol]);
		}

		// get the number of integration points
		int nint = el.GaussPoints();

		// loop over the integration points
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
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
			for (int isol = 0; isol<nsol; ++isol) {
				ps.m_c[isol] = el.Evaluate(c0[isol], n);
				ps.m_gradc[isol] = gradient(el, c0[isol], n);
			}

			ps.m_psi = pmb->ElectricPotential(mp);
			for (int isol = 0; isol<nsol; ++isol) {
				ps.m_ca[isol] = pmb->Concentration(mp, isol);
				ps.m_j[isol] = pmb->SoluteFlux(mp, isol);
				ps.m_crp[isol] = pm.m_J*m_pMat->Porosity(mp)*ps.m_ca[isol];
			}
			pt.m_pa = pmb->Pressure(mp);
			ps.m_cF = pmb->FixedChargeDensity(mp);
			ps.m_Ie = pmb->CurrentDensity(mp);

			pm.m_s = pmb->Stress(mp);

		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Reset()
{
	// reset base class
	FESolidDomain::Reset();
	
	// get the multiphasic material
	FETriphasic* pmb = m_pMat;
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
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
			
			// initialize multiphasic solutes
			ps.m_nsol = nsol;
			ps.m_c.assign(nsol,0);
			ps.m_ca.assign(nsol,0);
            ps.m_crp.assign(nsol, 0);
			ps.m_gradc.assign(nsol,vec3d(0,0,0));
			ps.m_k.assign(nsol, 0);
			ps.m_dkdJ.assign(nsol, 0);
			ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
			ps.m_j.assign(nsol,vec3d(0,0,0));
			ps.m_nsbm = nsbm;
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InitElements()
{
	FESolidDomain::InitElements();

	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
	FEMesh& m = *GetMesh();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int i=0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
			xt[i] = m.Node(el.m_node[i]).m_rt;
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			pt.m_r0 = r0;
			pt.m_rt = rt;

			pt.m_J = defgrad(el, pt.m_F, j);

			mp.Init(false);
		}
	}
	
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
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
			// reset referential solid volume fraction at previous time
			pt.m_phi0p = pt.m_phi0;
            
            // reset determinant of solid deformation gradient at previous time
            pt.m_Jp = pe.m_J;
            
            // reset referential actual solute concentration at previous time
            for (int j=0; j<2; ++j) {
                ps.m_crp[j] = pe.m_J*m_pMat->Porosity(mp)*ps.m_ca[j];
            }
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalForces(FEGlobalVector& R)
{
	size_t NE = m_Elem.size();
	#pragma omp parallel for shared (NE)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> lm;
		
		// get the element
		FESolidElement& el = m_Elem[i];

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		ElementInternalForce(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		#pragma omp critical
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FETriphasicDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;
	mat3ds s;

	const double* Gr, *Gs, *Gt;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		detJt = invjact(el, Ji, n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		s = pt.m_s;

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

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalSoluteWorkSS(vector<double>& R, double dt)
{
	size_t NE = m_Elem.size();
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();

		// element force vector
		vector<double> fe(neln);
		
		// unpack the element
		vector<int> elm;
		UnpackLM(el, elm);
			
		// calculate cation internal work
		ElementInternalSoluteWorkSS(el, fe, dt, 0);
			
		// add solute work to global residual
        #pragma omp critical
        {
			for (int j=0; j<neln; ++j)
			{
				int J = elm[7*neln+j];
				if (J >= 0) R[J] += fe[j];
			}
		}

		// calculate anion internal work
		ElementInternalSoluteWorkSS(el, fe, dt, 1);
			
		// add solute work to global residual
        #pragma omp critical
        {
            for (int j=0; j<neln; ++j)
            {
                int J = elm[8*neln+j];
                if (J >= 0) R[J] += fe[j];
            }
        }
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalSoluteWork(vector<double>& R, double dt)
{
	size_t NE = m_Elem.size();
	for (size_t i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// unpack the element
		vector<int> elm;
		UnpackLM(el, elm);
			
		// the element force vector
		vector<double> fe(neln);
			
		// calculate cation internal work
		ElementInternalSoluteWork(el, fe, dt, 0);
			
		// add solute work to global residual
		for (int j=0; j<neln; ++j)
		{
			int J = elm[7*neln+j];
			if (J >= 0) R[J] += fe[j];
		}

		// calculate anion internal work
		ElementInternalSoluteWork(el, fe, dt, 1);
			
		// add solute work to global residual
		for (int j=0; j<neln; ++j)
		{
			int J = elm[8*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalFluidWorkSS(vector<double>& R, double dt)
{
	size_t NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// unpack the element
		vector<int> elm;
		UnpackLM(el, elm);
			
		// calculate fluid internal work
		vector<double> fe(neln);
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
	size_t NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// unpack the element
		vector<int> elm;
		UnpackLM(el, elm);
			
		// calculate fluid internal work
		vector<double> fe(neln);
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
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt, *H;
	
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
	// gauss-weights
	double* wg = el.GaussWeights();
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
		H = el.H(n);
		
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
		
		// next we get the determinant
		double Jp = ppt.m_Jp;
		double J = ept.m_J;
		
		// and then finally
		double divv = ((J-Jp)/dt)/J;
		
		// get the flux
		vec3d& w = ppt.m_w;
		
		// update force vector
		for (i=0; i<neln; ++i)
			fe[i] -= dt*(gradN[i]*w - divv*H[i])*detJ*wg[n];
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
	
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
	// gauss-weights
	double* wg = el.GaussWeights();
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEBiphasicMaterialPoint& ppt = *(el.GetMaterialPoint(n)->ExtractData<FEBiphasicMaterialPoint>());
		
		// calculate jacobian
		detJ = invjact(el, Ji, n);
		
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);
		
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
		// get the flux
		vec3d& w = ppt.m_w;
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(gradN[i]*w)*detJ*wg[n];
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
	double Ji[3][3], detJ;
	
	double *Gr, *Gs, *Gt, *H;
	
	// gradient of shape functions
	vector<vec3d> gradN(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FETriphasic* pm = m_pMat;
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FESolutesMaterialPoint& spt = *(el.GetMaterialPoint(n)->ExtractData<FESolutesMaterialPoint>());
		
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
            // save spatial gradient of shape functions
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
		
		// next we get the determinant
		double J = ept.m_J;
		
		// get the solute flux
		vec3d j[2] = {spt.m_j[0],spt.m_j[1]};
        // get the charge number
        int z[2] = {pm->m_pSolute[0]->ChargeNumber(),
            pm->m_pSolute[1]->ChargeNumber()};
        
        vec3d je = j[0]*z[0] + j[1]*z[1];

		// evaluate the porosity, its derivative w.r.t. J, and its gradient
		double phiw = pm->Porosity(mp);
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
            fe[i] -= dt*(gradN[i]*(j[ion]+je*pm->m_penalty)
                         - H[i]*((phiw*spt.m_ca[ion] - spt.m_crp[ion]/J)/dt)
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
	
	// gradient of shape functions
	vector<vec3d> gradN(neln);
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// get the element's material
	FETriphasic* pm = m_pMat;
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
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

        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
		// get the solute flux
		vec3d j[2] = {spt.m_j[0],spt.m_j[1]};

		// get the charge number
		int z[2] = {pm->m_pSolute[0]->ChargeNumber(),
			pm->m_pSolute[1]->ChargeNumber()};

        vec3d je = j[0]*z[0] + j[1]*z[1];
        
		// update force vector
		for (i=0; i<neln; ++i)
			fe[i] -= dt*(gradN[i]*(j[ion]+je*pm->m_penalty))*detJ*wg[n];
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp)
{
	FETriphasic* pm = m_pMat;

	// repeat over all solid elements
	size_t NE = m_Elem.size();
    
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
			lm[ndpn*i+4] = elm[7*neln+i];
			lm[ndpn*i+5] = elm[8*neln+i];
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
	FETriphasic* pm = m_pMat;
	
	// repeat over all solid elements
	size_t NE = m_Elem.size();
    
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
			lm[ndpn*i+4] = elm[7*neln+i];
			lm[ndpn*i+5] = elm[8*neln+i];
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
    int i, j, isol, jsol, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    FETriphasic* pm = m_pMat;
    const int nsol = 2;
    int ndpn = 4+nsol;
    
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
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        vec3d w = ppt.m_w;
        vec3d gradp = ppt.m_gradp;
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (isol=0; isol<nsol; ++isol)
            z[isol] = pm->m_pSolute[isol]->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        vector< vector<double> > dkdr(spt.m_dkdr);
        vector< vector<double> > dkdJr(spt.m_dkdJr);
        vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
        
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
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
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
                qpu = -gradN[j]*(1.0/dt);
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
                ke[ndpn*i+3][ndpn*j+3] += (- gradN[i]*(Ke*gradN[j]))*(tmp*dt);
                
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
                    qcu[isol] = qpu*(c[isol]*(kappa[isol]+J*phiw*dkdJ[isol]));
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
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*c[isol]*dkdc[isol][jsol];
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*(c[isol]*dkdc[isol][jsol] + kappa[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                             + H[i]*(qcc[isol][jsol]))*(tmp*dt);
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
    int i, j, isol, jsol, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    FETriphasic* pm = m_pMat;
    const int nsol = 2;
    int ndpn = 4+nsol;
    
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
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        vec3d w = ppt.m_w;
        vec3d gradp = ppt.m_gradp;
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (isol=0; isol<nsol; ++isol)
            z[isol] = pm->m_pSolute[isol]->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        vector< vector<double> > dkdr(spt.m_dkdr);
        vector< vector<double> > dkdJr(spt.m_dkdJr);
        vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
        
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
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
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
        vector<vec3d> gc(nsol),qcu(nsol),wc(nsol),jce(nsol);
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
                ke[ndpn*i+3][ndpn*j+3] += (- gradN[i]*(Ke*gradN[j]))*(tmp*dt);
                
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
	FETriphasic* pmat = m_pMat;
	int id0 = pmat->m_pSolute[0]->GetSoluteID();
	int id1 = pmat->m_pSolute[1]->GetSoluteID();
	
	// nodal concentrations
	double ct[2][FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) {
		ct[0][i] = m_pMesh->Node(el.m_node[i]).get(DOF_C + id0);
		ct[1][i] = m_pMesh->Node(el.m_node[i]).get(DOF_C + id1);
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
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		
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
//! calculates element's geometrical stiffness component for integration point n
void FETriphasicDomain::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
{
	int n, i, j;

	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];
	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// nr of nodes
	int neln = el.Nodes();

	// nr of integration points
	int nint = el.GaussPoints();

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// calculate geometrical element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

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

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds& s = pt.m_s;

		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
					   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
					   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
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
			#pragma omp critical
			{
				berr = true;
				if (NegativeJacobian::m_boutput) e.print();
			}
		}
	}

	// if we encountered an error, we request a running restart
	if (berr)
	{
		if (NegativeJacobian::m_boutput == false) felog.printbox("ERROR", "Negative jacobian was detected.");
		throw DoRunningRestart();
	}
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
		pn[j] = mesh.Node(el.m_node[j]).get(m_dofP);
		ct[0][j] = mesh.Node(el.m_node[j]).get(DOF_C + id0);
		ct[1][j] = mesh.Node(el.m_node[j]).get(DOF_C + id1);
	}
		
	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
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
        pmb->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc);
			
		pt.m_s = pmb->Stress(mp);
	}
}
