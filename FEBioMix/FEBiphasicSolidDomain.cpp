#include "FEBiphasicSolidDomain.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEBiphasicSolidDomain::FEBiphasicSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_BIPHASIC_DOMAIN, pm)
{
	m_pMat = dynamic_cast<FEBiphasic*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEBiphasicSolidDomain::InitElements()
{
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
}

//-----------------------------------------------------------------------------
bool FEBiphasicSolidDomain::Initialize(FEModel &fem)
{
	// initialize base class
	FESolidDomain::Initialize(fem);
    
    // initialize body forces
    m_pMat->m_bf.clear();
    for (int j=0; j<fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
        if (pbf) m_pMat->m_bf.push_back(pbf);
    }

	// initialize local coordinate systems (can I do this elsewhere?)
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
	}

	// initialize all element data
	const int NE = FEElement::MAX_NODES;
    double p0[NE];
	FEMesh& m = *GetMesh();
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
		for (int i=0; i<neln; ++i)
			p0[i] = m.Node(el.m_node[i]).m_p0;
        
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEElasticMaterialPoint&  pm = *(mp.ExtractData<FEElasticMaterialPoint >());
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
            // initialize effective fluid pressure, its gradient, and fluid flux
            pt.m_p = el.Evaluate(p0, n);
            pt.m_gradp = gradient(el, p0, n);
            pt.m_w = m_pMat->Flux(mp);
            pt.m_pa = m_pMat->Pressure(mp);
            pm.m_s = m_pMat->Stress(mp);
           
			// initialize referential solid volume fraction
			pt.m_phi0 = m_pMat->m_phi0;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEBiphasicSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	int N = el.Nodes();
	lm.resize(N*MAX_NDOFS);
	
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[0];
		lm[3*i+1] = id[1];
		lm[3*i+2] = id[2];

		// now the pressure dofs
		lm[3*N+i] = id[6];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[7];
		lm[4*N + 3*i+1] = id[8];
		lm[4*N + 3*i+2] = id[9];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;

		lm[10*N + i] = id[10];
		
		// concentration dofs
		for (int k=0; k<MAX_CDOFS; ++k)
			lm[(11+k)*N + i] = id[11+k];
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Reset()
{
	// reset base class data
	FESolidDomain::Reset();

	// get the biphasic material
	FEBiphasic* pmb = m_pMat;

	// initialize all element data
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
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
		}
	}
}


//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::InternalForces(FEGlobalVector& R)
{
	int NE = m_Elem.size();
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
		//#pragma omp critical
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
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
// Calculate the work due to the internal fluid pressure (for steady-state analysis)
void FEBiphasicSolidDomain::InternalFluidWorkSS(vector<double>& R, double dt)
{
	int NE = m_Elem.size();

	#pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		
		// calculate fluid internal work
		vector<double> fe(neln);
		ElementInternalFluidWorkSS(el, fe, dt);
			
		// assemble element 'fe'-vector into global R vector
		vector<int> elm;
		UnpackLM(el, elm);
		
		//#pragma omp critical
		{
			// add fluid work to global residual
			int neln = el.Nodes();
			for (int j=0; j<neln; ++j)
			{
				int J = elm[3*neln+j];
				if (J >= 0)
					{
#pragma omp atomic
						R[J] += fe[j];
				}
			}
		}
	}
}


//-----------------------------------------------------------------------------
// Calculate the work due to the internal fluid pressure
void FEBiphasicSolidDomain::InternalFluidWork(vector<double>& R, double dt)
{
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
			
		// calculate fluid internal work
		vector<double> fe(neln);
		ElementInternalFluidWork(el, fe, dt);

		// unpack the element
		vector<int> elm;
		UnpackLM(el, elm);
			
        {
            // add fluid work to global residual
            int neln = el.Nodes();
            for (int j=0; j<neln; ++j)
            {
                int J = elm[3*neln+j];
                if (J >= 0)
				{
		#pragma omp atomic
					R[J] += fe[j];
				}
            }
        }
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
bool FEBiphasicSolidDomain::ElementInternalFluidWork(FESolidElement& el, vector<double>& fe, double dt)
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
	
	FEMesh& mesh = *GetMesh();
	
	vec3d rp[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
	}
	
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		
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
		vec3d& w = pt.m_w;

		// get the solvent supply
		double phiwhat = m_pMat->SolventSupply(mp);
		
		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z + (phiwhat - divv)*H[i])*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! for a steady-state analysis (solid velocity = 0)

bool FEBiphasicSolidDomain::ElementInternalFluidWorkSS(FESolidElement& el, vector<double>& fe, double dt)
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
		
	zero(fe);
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		
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
		vec3d& w = pt.m_w;

		// get the solvent supply
		double phiwhat = m_pMat->SolventSupply(mp);

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z + H[i]*phiwhat)*detJ*wg[n];
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrix(FESolver* psolver, bool bsymm, double dt)
{
	// repeat over all solid elements
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> elm;
		
		FESolidElement& el = m_Elem[iel];
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffness(el, ke, bsymm, dt);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
		UnpackLM(el, elm);

        vector<int> lm(ndof);
        for (int i=0; i<neln; ++i)
        {
            lm[4*i  ] = elm[3*i];
            lm[4*i+1] = elm[3*i+1];
            lm[4*i+2] = elm[3*i+2];
            lm[4*i+3] = elm[3*neln+i];
        }
        
        // assemble element matrix in global stiffness matrix
        #pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrixSS(FESolver* psolver, bool bsymm, double dt)
{
	// repeat over all solid elements
	int NE = m_Elem.size();

	#pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> elm;
		
		FESolidElement& el = m_Elem[iel];
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffnessSS(el, ke, bsymm, dt);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
		UnpackLM(el, elm);

		vector<int> lm(ndof);
		for (int i=0; i<neln; ++i)
		{
			lm[4*i  ] = elm[3*i];
			lm[4*i+1] = elm[3*i+1];
			lm[4*i+2] = elm[3*i+2];
			lm[4*i+3] = elm[3*neln+i];
		}
		
		// assemble element matrix in global stiffness matrix
		#pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEBiphasicSolidDomain::ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
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
	
	FEMesh& mesh = *GetMesh();
	
	const int NE = FEElement::MAX_NODES;
	vec3d r0[NE], rt[NE], rp[NE], v[NE];
	for (i=0; i<neln; ++i) 
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		v[i]  = mesh.Node(el.m_node[i]).m_vt;
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
			ke[4*i  ][4*j] = ks[3*i  ][3*j  ]; ke[4*i  ][4*j+1] = ks[3*i  ][3*j+1]; ke[4*i  ][4*j+2] = ks[3*i  ][3*j+2];
			ke[4*i+1][4*j] = ks[3*i+1][3*j  ]; ke[4*i+1][4*j+1] = ks[3*i+1][3*j+1]; ke[4*i+1][4*j+2] = ks[3*i+1][3*j+2];
			ke[4*i+2][4*j] = ks[3*i+2][3*j  ]; ke[4*i+2][4*j+1] = ks[3*i+2][3*j+1]; ke[4*i+2][4*j+2] = ks[3*i+2][3*j+2];
		}
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		
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
		vec3d vs(0,0,0);
		
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
			
		}
		
		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.m_J;
		
		// and then finally
		double dJdt = (J-Jp)/dt;
		double divv = dJdt/J;
		
		// get the fluid flux and pressure gradient
		vec3d w = pt.m_w;
		vec3d gradp = pt.m_gradp;
		
		// evaluate the permeability and its derivatives
		mat3ds K = m_pMat->Permeability(mp);
		tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);

		// evaluate the solvent supply and its derivatives
		double phiwhat = 0;
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		if (m_pMat->GetSolventSupply()) {
			phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
			Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
			Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
		}
		
		// Miscellaneous constants
		mat3dd I(1);
		
		// calculate the kpp matrix
		//		tmp = detJ*gw[n];
		tmp = detJ*gw[n]*dt;
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				ke[4*i+3][4*j+3] += (H[i]*H[j]*Phip - gradN[i]*(K*gradN[j]))*tmp;
			}
		
		if (!bsymm) {
			// calculate the kup matrix
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					tmp = detJ*gw[n]*H[j];
					ke[4*i  ][4*j+3] -= tmp*gradN[i].x;
					ke[4*i+1][4*j+3] -= tmp*gradN[i].y;
					ke[4*i+2][4*j+3] -= tmp*gradN[i].z;
				}
			}
			
			// calculate the kpu matrix
			//			tmp = detJ*gw[n];
			tmp = detJ*gw[n]*dt;
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					vec3d vt = (-(vdotTdotv(gradN[i], dKdE, gradN[j])*(gradp))
								-(I*(divv+1./dt) - gradv.transpose() - Phie)*gradN[j]*H[i])*tmp;
					ke[4*i+3][4*j  ] += vt.x;
					ke[4*i+3][4*j+1] += vt.y;
					ke[4*i+3][4*j+2] += vt.z;
				}
			}
			
		} else {
			// calculate the kup matrix and let kpu be its symmetric part
			tmp = detJ*gw[n];
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					ke[4*i  ][4*j+3] -= tmp*H[j]*gradN[i].x;
					ke[4*i+1][4*j+3] -= tmp*H[j]*gradN[i].y;
					ke[4*i+2][4*j+3] -= tmp*H[j]*gradN[i].z;
					
					ke[4*i+3][4*j  ] -= tmp*H[i]*gradN[j].x;
					ke[4*i+3][4*j+1] -= tmp*H[i]*gradN[j].y;
					ke[4*i+3][4*j+2] -= tmp*H[i]*gradN[j].z;
				}
			}
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for the steady-state response (zero solid velocity)
//!
bool FEBiphasicSolidDomain::ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
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
			ke[4*i  ][4*j] = ks[3*i  ][3*j  ]; ke[4*i  ][4*j+1] = ks[3*i  ][3*j+1]; ke[4*i  ][4*j+2] = ks[3*i  ][3*j+2];
			ke[4*i+1][4*j] = ks[3*i+1][3*j  ]; ke[4*i+1][4*j+1] = ks[3*i+1][3*j+1]; ke[4*i+1][4*j+2] = ks[3*i+1][3*j+2];
			ke[4*i+2][4*j] = ks[3*i+2][3*j  ]; ke[4*i+2][4*j+1] = ks[3*i+2][3*j+1]; ke[4*i+2][4*j+2] = ks[3*i+2][3*j+2];
		}
	
	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		
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
		
		// get the fluid flux and pressure gradient
		vec3d w = pt.m_w;
		vec3d gradp = pt.m_gradp;
		
		// evaluate the permeability and its derivatives
		mat3ds K = m_pMat->Permeability(mp);
		tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);

		// evaluate the solvent supply and its derivatives
		double phiwhat = 0;
		mat3ds Phie; Phie.zero();
		double Phip = 0;
		if (m_pMat->GetSolventSupply()) {
			phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
			Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
			Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
		}

		// Miscellaneous constants
		mat3dd I(1);
		
		// calculate the kpp matrix
		//		tmp = detJ*gw[n];
		tmp = detJ*gw[n]*dt;
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				ke[4*i+3][4*j+3] += (H[i]*H[j]*Phip - gradN[i]*(K*gradN[j]))*tmp;
			}
		
		if (!bsymm) {
			// calculate the kup matrix
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					tmp = detJ*gw[n]*H[j];
					ke[4*i  ][4*j+3] -= tmp*gradN[i].x;
					ke[4*i+1][4*j+3] -= tmp*gradN[i].y;
					ke[4*i+2][4*j+3] -= tmp*gradN[i].z;
				}
			}
			
			// calculate the kpu matrix
			//			tmp = detJ*gw[n];
			tmp = detJ*gw[n]*dt;
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					vec3d vt = (-(vdotTdotv(gradN[i], dKdE, gradN[j])*(gradp))
								+Phie*gradN[j]*H[i])*tmp;
					ke[4*i+3][4*j  ] += vt.x;
					ke[4*i+3][4*j+1] += vt.y;
					ke[4*i+3][4*j+2] += vt.z;
				}
			}
			
		} else {
			// calculate the kup matrix and let kpu be its symmetric part
			tmp = detJ*gw[n];
			for (i=0; i<neln; ++i) {
				for (j=0; j<neln; ++j)
				{
					ke[4*i  ][4*j+3] -= tmp*H[j]*gradN[i].x;
					ke[4*i+1][4*j+3] -= tmp*H[j]*gradN[i].y;
					ke[4*i+2][4*j+3] -= tmp*H[j]*gradN[i].z;
					
					ke[4*i+3][4*j  ] -= tmp*H[i]*gradN[j].x;
					ke[4*i+3][4*j+1] -= tmp*H[i]*gradN[j].y;
					ke[4*i+3][4*j+2] -= tmp*H[i]*gradN[j].z;
				}
			}
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

void FEBiphasicSolidDomain::SolidElementStiffness(FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	ElementBiphasicMaterialStiffness(el, ke);
	
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

void FEBiphasicSolidDomain::ElementBiphasicMaterialStiffness(FESolidElement &el, matrix &ke)
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
	
	// weights at gauss points
	const double *gw = el.GaussWeights();
	
	// nodal pressures
	double pn[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) pn[i] = m_pMesh->Node(el.m_node[i]).m_pt;
	
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
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// evaluate fluid pressure at gauss-point
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		ppt.m_p = el.Evaluate(pn, n);
		
		// get the 'D' matrix
		tens4ds C = m_pMat->Tangent(mp);
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
void FEBiphasicSolidDomain::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
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
void FEBiphasicSolidDomain::UpdateStresses(FEModel &fem)
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
//			felog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
			#pragma omp critical
			berr = true;
		}
	}
	// if we encountered an error, we request a running restart
	if (berr)
	{
		felog.printbox("ERROR","Negative jacobian was detected.");
		throw DoRunningRestart();
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::UpdateElementStress(int iel)
{
   // extract the elastic component
    FEElasticMaterial* pme = m_pMat->GetElasticMaterial();

	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	// get the integration weights
	double* gw = el.GaussWeights();

	// get the number of nodes
	int neln = el.Nodes();

	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
		pn[j] = mesh.Node(el.m_node[j]).m_pt;
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
			
		// poroelasticity data
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
		// evaluate fluid pressure at gauss-point
		ppt.m_p = el.Evaluate(pn, n);
			
		// calculate the gradient of p at gauss-point
		ppt.m_gradp = gradient(el, pn, n);
			
		// for biphasic materials also update the fluid flux
		ppt.m_w = m_pMat->Flux(mp);
		ppt.m_pa = m_pMat->Pressure(mp);
			
		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);

		int a = 0;
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 3*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEBiphasicSolidDomain::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
    double rhoTw = m_pMat->FluidDensity();
    
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d b;
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
        rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    }
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
        pt.m_r0 = el.Evaluate(r0, n);
        pt.m_rt = el.Evaluate(rt, n);
        
        detJ = detJt(el, n)*gw[n];
        
        // get the force
        b = BF.force(mp);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[3*i  ] -= H[i]*rho*b.x*detJ;
            fe[3*i+1] -= H[i]*rho*b.y*detJ;
            fe[3*i+2] -= H[i]*rho*b.z*detJ;
        }						
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
    FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(GetMaterial()); assert(pmb);
    
    // element stiffness matrix
    matrix ke;
    vector<int> elm;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int neln = el.Nodes();
        int ndof = 4*neln;
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // TODO: the problem here is that the LM array that is returned by the UnpackLM
        // function does not give the equation numbers in the right order. For this reason we
        // have to create a new lm array and place the equation numbers in the right order.
        // What we really ought to do is fix the UnpackLM function so that it returns
        // the LM vector in the right order for poroelastic elements.
        UnpackLM(el, elm);
        
        vector<int> lm(ndof);
        for (int i=0; i<neln; ++i)
        {
            lm[4*i  ] = elm[3*i];
            lm[4*i+1] = elm[3*i+1];
            lm[4*i+2] = elm[3*i+2];
            lm[4*i+3] = elm[3*neln+i];
        }
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEBiphasicSolidDomain::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    int neln = el.Nodes();
    int ndof = ke.columns()/neln;
    
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
    double rhoTw = m_pMat->FluidDensity();
    
    // jacobian
    double detJ, detJt, Ji[3][3];
    double *H;
    double* gw = el.GaussWeights();
    vec3d G[FEElement::MAX_NODES];
    double *Grn, *Gsn, *Gtn;
    double Gr, Gs, Gt;
    
    vec3d b, kpu;
    mat3ds Kb;
    mat3d Kw, Kuu;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        detJ = detJ0(el, n)*gw[n];
        
        // get the body force
        b = BF.force(mp);
        
        // get the body force stiffness
        Kb = BF.stiffness(mp);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        H = el.H(n);
        
        // calculate jacobian
        detJt = invjact(el, Ji, n)*gw[n];
        
        Grn = el.Gr(n);
        Gsn = el.Gs(n);
        Gtn = el.Gt(n);
        
        for (int i=0; i<neln; ++i)
        {
            Gr = Grn[i];
            Gs = Gsn[i];
            Gt = Gtn[i];
            
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            G[i].x = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
            G[i].y = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
            G[i].z = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
        }
        
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                Kw = b & G[j];
                Kuu = (Kb*(H[j]*rho) + Kw*rhoTw)*(H[i]*detJ);
                ke[ndof*i  ][ndof*j  ] -= Kuu(0,0);
                ke[ndof*i  ][ndof*j+1] -= Kuu(0,1);
                ke[ndof*i  ][ndof*j+2] -= Kuu(0,2);
                
                ke[ndof*i+1][ndof*j  ] -= Kuu(1,0);
                ke[ndof*i+1][ndof*j+1] -= Kuu(1,1);
                ke[ndof*i+1][ndof*j+2] -= Kuu(1,2);
                
                ke[ndof*i+2][ndof*j  ] -= Kuu(2,0);
                ke[ndof*i+2][ndof*j+1] -= Kuu(2,1);
                ke[ndof*i+2][ndof*j+2] -= Kuu(2,2);
                
                kpu = (vdotTdotv(G[i], dKdE, G[j])*b
                       + (Kw + Kb*H[j])*K*G[i])*(rhoTw*detJ);
                ke[ndof*i+3][ndof*j  ] += kpu.x;
                ke[ndof*i+3][ndof*j+1] += kpu.y;
                ke[ndof*i+3][ndof*j+2] += kpu.z;
            }
    }	
}
