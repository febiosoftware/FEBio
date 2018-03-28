#include "FEThermoElasticSolidDomain.h"
#include "FECore/FEMesh.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEHeatTransferMaterial.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEThermoElasticSolidDomain::FEThermoElasticSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticDomain(pfem)
{
	m_pMat = 0;
	m_dofT = pfem->GetDOFIndex("T");
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FEThermoElasticMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEThermoElasticSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
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

			mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
bool FEThermoElasticSolidDomain::Init()
{
	// initialize base class
	FESolidDomain::Init();
	FEModel& fem = *GetFEModel();
	const int dof_T = fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return false; }
    
	// initialize local coordinate systems (can I do this elsewhere?)
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
	}

	// initialize all element data
	const int NE = FEElement::MAX_NODES;
    double T0[NE];
	FEMesh& m = *GetMesh();
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of temperature
		for (int i=0; i<neln; ++i)
			T0[i] = m.Node(el.m_node[i]).get(dof_T);
        
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEElasticMaterialPoint&  pm = *(mp.ExtractData<FEElasticMaterialPoint >());
			
            // initialize stress
            pm.m_s.zero(); // m_pMat->Stress(mp);
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[m_dofX] = DOF_ACTIVE;
				node.m_ID[m_dofY] = DOF_ACTIVE;
				node.m_ID[m_dofZ] = DOF_ACTIVE;
			}

			node.m_ID[m_dofT] = DOF_ACTIVE;
		}
	}
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEThermoElasticSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.assign(N*4, -1);
	
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];

		// now the temperature dofs
		lm[3*N+i] = id[m_dofT];
	}
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::Reset()
{
	// reset base class data
	FESolidDomain::Reset();
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::InternalForces(FEGlobalVector& R)
{
	int NE = (int)m_Elem.size();
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

void FEThermoElasticSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
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
// Calculate the work due to the heat flux
void FEThermoElasticSolidDomain::InternalThermalWork(vector<double>& R)
{
	int NE = (int)m_Elem.size();

	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		const int neln = el.Nodes();
		
		// calculate thermal internal work
		vector<double> fe(neln);
		ElementInternalThermalWork(el, fe);
			
		// assemble element 'fe'-vector into global R vector
		vector<int> elm;
		UnpackLM(el, elm);
		
		// add forces to global residual
		for (int j=0; j<neln; ++j)
		{
			int J = elm[3*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
bool FEThermoElasticSolidDomain::ElementInternalThermalWork(FESolidElement& el, vector<double>& fe)
{
	// jacobian
	double Ji[3][3];
	
	// gauss-weights
	double* wg = el.GaussWeights();
	
	// zero force vector
	zero(fe);
	
	// loop over gauss-points
	const int neln = el.Nodes();
	const int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEHeatMaterialPoint& pt = *(mp.ExtractData<FEHeatMaterialPoint >());

		// get the heat flux vector
		vec3d q = pt.m_q;
		
		// calculate jacobian
		double detJ = invjact(el, Ji, n);

		// loop over all nodes
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);
		double* Gt = el.Gt(n);
		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// update thermal "force" vector
			fe[i] += detJ*wg[n]*(q.x*Gx + q.y*Gy + q.z*Gz);
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::StiffnessMatrix(FESolver* psolver)
{
	// repeat over all solid elements
	int NE = (int)m_Elem.size();
    
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
		ElementStiffness(el, ke);
		
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
void FEThermoElasticSolidDomain::ElementStiffness(FESolidElement& el, matrix& ke)
{
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// get the solid stiffness
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(el, ks);

	// calculate the thermal tangent stiffness
	matrix kt(ndof, neln); kt.zero();
	ElementThermalStiffness(el, kt);

	// calculate the thermal conductivity stiffness
	matrix kc(neln, neln); kc.zero();
	ElementConductionStiffness(el, kc);

	// calculate the conducitivity gradient stiffness
	matrix kg(neln, ndof); kg.zero();
	ElementGradientStiffness(el, kg);
	
	// Compose all the sub-matrices into a single element stiffness matrix
	for (int i=0; i<neln; ++i)
		for (int j=0; j<neln; ++j)
		{
			// expand solid stiffess
			ke[4*i  ][4*j] = ks[3*i  ][3*j  ]; ke[4*i  ][4*j+1] = ks[3*i  ][3*j+1]; ke[4*i  ][4*j+2] = ks[3*i  ][3*j+2];
			ke[4*i+1][4*j] = ks[3*i+1][3*j  ]; ke[4*i+1][4*j+1] = ks[3*i+1][3*j+1]; ke[4*i+1][4*j+2] = ks[3*i+1][3*j+2];
			ke[4*i+2][4*j] = ks[3*i+2][3*j  ]; ke[4*i+2][4*j+1] = ks[3*i+2][3*j+1]; ke[4*i+2][4*j+2] = ks[3*i+2][3*j+2];

			// expand thermal stiffness
			ke[4*i  ][4*j+3] = kt[3*i  ][j];
			ke[4*i+1][4*j+3] = kt[3*i+1][j];
			ke[4*i+2][4*j+3] = kt[3*i+2][j];

			// expand flux gradient stiffness
			ke[4*i+3][4*j  ] = kg[i][3*j  ];
			ke[4*i+3][4*j+1] = kg[i][3*j+1];
			ke[4*i+3][4*j+2] = kg[i][3*j+2];

			// expand conductivity stiffness
			ke[4*i+3][4*j+3] = kc[i][j];
		}
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function. Note that these functions
//! only calculate the upper diagonal matrix due to the symmetry of the element
//! stiffness matrix. The last section of this function fills the rest of the element
//! stiffness matrix.
void FEThermoElasticSolidDomain::SolidElementStiffness(FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	ElementMaterialStiffness(el, ke);
	
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
//
void FEThermoElasticSolidDomain::ElementMaterialStiffness(FESolidElement &el, matrix &ke)
{
	int i, i3, j, j3;

	// global derivatives of shape functions
	// Gx = dH/dx
	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// jacobian
	double Ji[3][3];
	
	// calculate element stiffness matrix
	const double *gw = el.GaussWeights();
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian
		double detJt = invjact(el, Ji, n)*gw[n];

		double* Grn = el.Gr(n);
		double* Gsn = el.Gs(n);
		double* Gtn = el.Gt(n);

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// get the 'D' matrix
		tens4ds C = m_pMat->Tangent(mp);
		C.extract(D);

		for (i=0; i<neln; ++i)
		{
			double Gr = Grn[i];
			double Gs = Gsn[i];
			double Gt = Gtn[i];

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
			double Gxi = Gx[i];
			double Gyi = Gy[i];
			double Gzi = Gz[i];

			for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				double Gxj = Gx[j];
				double Gyj = Gy[j];
				double Gzj = Gz[j];

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
void FEThermoElasticSolidDomain::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
{
	// spatial shape function gradients
	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];

	// jacobian
	double Ji[3][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate geometrical element stiffness matrix
	const int neln = el.Nodes();
	const int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian
		double detJt = invjact(el, Ji, n)*gw[n];

		double* Grn = el.Gr(n);
		double* Gsn = el.Gs(n);
		double* Gtn = el.Gt(n);

		for (int i=0; i<neln; ++i)
		{
			double Gr = Grn[i];
			double Gs = Gsn[i];
			double Gt = Gtn[i];

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
		mat3ds& s = pt.m_s;

		for (int i=0; i<neln; ++i)
			for (int j=i; j<neln; ++j)
			{
				double kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
					          Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
					          Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::ElementConductionStiffness(FESolidElement &el, matrix& ke)
{
	// global derivatives of shape functions
	// Gx = dH/dx
	const int EN = FEElement::MAX_NODES;
	double Gx[EN], Gy[EN], Gz[EN];

	double Gi[3], Gj[3];
	double DB[3];

	// jacobian
	double Ji[3][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// loop over all integration points
	const int ne = el.Nodes();
	const int ni = el.GaussPoints();
	for (int n=0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = invjact(el, Ji, n);

		// evaluate the conductivity
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		mat3ds K = m_pMat->Conductivity(mp);

		for (int i=0; i<ne; ++i)
		{
			double Gr = el.Gr(n)[i];
			double Gs = el.Gs(n)[i];
			double Gt = el.Gt(n)[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}		

		for (int i=0; i<ne; ++i)
		{
			Gi[0] = Gx[i];
			Gi[1] = Gy[i];
			Gi[2] = Gz[i];

			for (int j=0; j<ne; ++j)
			{
				Gj[0] = Gx[j];
				Gj[1] = Gy[j];
				Gj[2] = Gz[j];

				DB[0] = K(0,0)*Gj[0] + K(0,1)*Gj[1] + K(0,2)*Gj[2];
				DB[1] = K(1,0)*Gj[0] + K(1,1)*Gj[1] + K(1,2)*Gj[2];
				DB[2] = K(2,0)*Gj[0] + K(2,1)*Gj[1] + K(2,2)*Gj[2];

				ke[i][j] += (Gi[0]*DB[0] + Gi[1]*DB[1] + Gi[2]*DB[2] )*detJt*gw[n];
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Thermal stiffness (contribution of the dependency of stress on temperature)
void FEThermoElasticSolidDomain::ElementThermalStiffness(FESolidElement &el, matrix& ke)
{
	// global derivatives of shape functions
	// Gx = dH/dx
	const int EN = FEElement::MAX_NODES;
	double Gx[EN], Gy[EN], Gz[EN];

	double Gi[3];
	double BT[3];

	// jacobian
	double Ji[3][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// loop over all integration points
	const int ne = el.Nodes();
	const int ni = el.GaussPoints();
	for (int n=0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = invjact(el, Ji, n);

		// evaluate the thermal tangent
		// (i.e. the derivative of the stress with respect with the temperature)
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		mat3ds T = m_pMat->ThermalTangent(mp);

		for (int i=0; i<ne; ++i)
		{
			double Gr = el.Gr(n)[i];
			double Gs = el.Gs(n)[i];
			double Gt = el.Gt(n)[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}		

		// loop over displacement dofs
		for (int i=0; i<ne; ++i)
		{
			Gi[0] = Gx[i];
			Gi[1] = Gy[i];
			Gi[2] = Gz[i];

			// loop over temperature dofs
			for (int j=0; j<ne; ++j)
			{
				double Hj = el.H(n)[j];

				BT[0] = T(0,0)*Gi[0] + T(0,1)*Gi[1] + T(0,2)*Gi[2];
				BT[1] = T(1,0)*Gi[0] + T(1,1)*Gi[1] + T(1,2)*Gi[2];
				BT[2] = T(2,0)*Gi[0] + T(2,1)*Gi[1] + T(2,2)*Gi[2];

				ke[3*i  ][j] += (BT[0]*Hj)*detJt*gw[n];
				ke[3*i+1][j] += (BT[1]*Hj)*detJt*gw[n];
				ke[3*i+2][j] += (BT[2]*Hj)*detJt*gw[n];
			}
		}
	}
}

//-----------------------------------------------------------------------------
// implements the weird product.
vec3d weird_product(double Ga[3], double Gb[3], double GT[3], const tens4ds& t)
{
	vec3d r(0,0,0);

	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
			{
				double G3 = Ga[i]*GT[j]*Gb[k];
				r.x += G3*t(i,j,k,0);
				r.y += G3*t(i,j,k,1);
				r.z += G3*t(i,j,k,2);
			}

	return r;
}

//-----------------------------------------------------------------------------
//! Conductivity gradient stiffness (i.e. derivative of conductivity wrt strain
void FEThermoElasticSolidDomain::ElementGradientStiffness(FESolidElement &el, matrix& ke)
{
	const int dof_T = GetFEModel()->GetDOFS().GetDOF("T");

	// global derivatives of shape functions
	// Gx = dH/dx
	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// jacobian
	double Ji[3][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// current nodal temperatures
	FEMesh& mesh = *GetMesh();
	double T[EN];
	const int ne = el.Nodes();
	for (int i=0; i<ne; ++i) T[i] = mesh.Node(el.m_node[i]).get(dof_T);

	// loop over all integration points
	const int ni = el.GaussPoints();
	for (int n=0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = invjact(el, Ji, n);

		// evaluate the conductivity gradient
		// (i.e. the derivative of the conductivity with respect with to strain)
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		tens4ds D = m_pMat->ConductivityGradient(mp);

		// calculate the spatial gradients of the shape functions
		for (int i=0; i<ne; ++i)
		{
			double Gr = el.Gr(n)[i];
			double Gs = el.Gs(n)[i];
			double Gt = el.Gt(n)[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			G[i].x = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			G[i].y = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			G[i].z = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// evaluate the gradient at the integration point
		vec3d Tgrad = gradient(el, T, n);
		double GT[3] = {Tgrad.x, Tgrad.y, Tgrad.z};

		// loop over temperature dofs
		for (int i=0; i<ne; ++i)
		{
			double Ga[3] = {G[i].x, G[i].y, G[i].z};
			// loop over displacement dofs
			for (int j=0; j<ne; ++j)
			{
				double Gb[3] = {G[j].x, G[j].y, G[j].z};

				vec3d r = weird_product(Ga, Gb, GT, D);
				ke[i][3*j  ] += detJt*gw[n]*r.x;
				ke[i][3*j+1] += detJt*gw[n]*r.y;
				ke[i][3*j+2] += detJt*gw[n]*r.z;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEThermoElasticSolidDomain::Update(const FETimeInfo& tp)
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
// This function evaluates the state variables at the integration points of element iel.
// It evaluates the Cauchy stress tensor, as well as the spatial heat flux vector.
void FEThermoElasticSolidDomain::UpdateElementStress(int iel)
{
	const int dof_T = GetFEModel()->GetDOFS().GetDOF("T");

	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the nodal data
	FEMesh& mesh = *GetMesh();
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double u0[FEElement::MAX_NODES];
	double ut[FEElement::MAX_NODES];
	int neln = el.Nodes();
	for (int j=0; j<neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
		rt[j] = mesh.Node(el.m_node[j]).m_rt;

		// TODO: After I make the transition to domain specific data I need to fix this
//		u0[j] = mesh.Node(el.m_node[j]).m_T0;
//		ut[j] = mesh.Node(el.m_node[j]).get(dof_T);
		assert(false);
	}

	// loop over the integration points and calculate
	// the state data at the integration point
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			
		// material point coordinates
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);
			
		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
			
		// evaluate temperature
		FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());
		ht.m_T0 = el.Evaluate(u0, n);
		ht.m_T  = el.Evaluate(ut, n);

		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);

		// evaluate the conductivity
		mat3ds K = m_pMat->Conductivity(mp);

		// update the heat flux vector
		ht.m_q = -(K*gradient(el, ut, n));
	}
}
