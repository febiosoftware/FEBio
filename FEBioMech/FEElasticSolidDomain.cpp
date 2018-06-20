#include "stdafx.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticMaterial.h"
#include "FEBodyForce.h"
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEElasticSolidDomain::FEElasticSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticDomain(pfem)
{
	m_pMat = 0;
    m_alphaf = m_beta = 1;
    m_alpham = 2;
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEElasticSolidDomain& FEElasticSolidDomain::operator = (FEElasticSolidDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
//! Assign material
void FEElasticSolidDomain::SetMaterial(FEMaterial* pmat)
{
	if (pmat)
	{
		m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
		assert(m_pMat);
	}
	else m_pMat = 0;
}

//-----------------------------------------------------------------------------
//! \todo The material point initialization needs to move to the base class.
bool FEElasticSolidDomain::Init()
{
	// initialize base class
	if (FESolidDomain::Init() == false) return false;

	// get the elements material
	if (m_pMat)
	{
		FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
		if (pme)
		{
			// assign local coordinate system to each integration point
			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FESolidElement& el = m_Elem[i];
				for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
			}
		}
	}

	// check for initially inverted elements
	int ninverted = 0;
	for (int i=0; i<Elements(); ++i)
	{
		FESolidElement& el = Element(i);

		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			double J0 = detJ0(el, n);
			if (J0 <= 0)
			{
				felog.printf("**************************** E R R O R ****************************\n");
				felog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.GetID());
				felog.printf("Jacobian = %lg\n", J0);
				felog.printf("Did you use the right node numbering?\n");
				felog.printf("Nodes:");
				for (int l=0; l<el.Nodes(); ++l)
				{
					felog.printf("%d", el.m_node[l]+1);
					if (l+1 != el.Nodes()) felog.printf(","); else felog.printf("\n");
				}
				felog.printf("*******************************************************************\n\n");
				++ninverted;
			}
		}
	}

	return (ninverted == 0);
}


//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Activate()
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
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEElasticSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    m_alphaf = timeInfo.alphaf;
    m_alpham = timeInfo.alpham;
    m_beta = timeInfo.beta;

	vec3d r0, rt;
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            pt.m_Wp = pt.m_Wt;

			mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::InternalForces(FEGlobalVector& R)
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

void FEElasticSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3];

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		double detJt = invjact(el, Ji, n, m_alphaf);

		detJt *= gw[n];

		// get the stress vector for this integration point
        mat3ds& s = pt.m_s;

		const double* Gr = el.Gr(n);
		const double* Gs = el.Gs(n);
		const double* Gt = el.Gt(n);

		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

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
void FEElasticSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();

	// TODO: Evaluate of body forces is not thread-safe due to the use of FEMathDouble in FENonConstBodyForce.
	//       I need to turn of parallelization until this issue is resolved.
//#pragma omp parallel for 
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

void FEElasticSolidDomain::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
	// don't forget to multiply with the density
	double dens = m_pMat->Density();

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	vec3d f;

	// number of nodes
	int neln = el.Nodes();

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		detJ = detJ0(el, n)*gw[n];

		// get the force
		f = BF.force(mp);

		H = el.H(n);

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] -= H[i]*dens*f.x*detJ;
			fe[3*i+1] -= H[i]*dens*f.y*detJ;
			fe[3*i+2] -= H[i]*dens*f.z*detJ;
		}						
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEElasticSolidDomain::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
	int neln = el.Nodes();
	int ndof = ke.columns()/neln;

	// don't forget to multiply with the density
	double dens = m_pMat->Density();

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	mat3ds K;

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		detJ = detJ0(el, n)*gw[n]*m_alphaf;

		// get the stiffness
		K = BF.stiffness(mp);

		H = el.H(n);

		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				ke[ndof*i  ][ndof*j  ] -= H[i]*H[j]*dens*K(0,0)*detJ;
				ke[ndof*i  ][ndof*j+1] -= H[i]*H[j]*dens*K(0,1)*detJ;
				ke[ndof*i  ][ndof*j+2] -= H[i]*H[j]*dens*K(0,2)*detJ;

				ke[ndof*i+1][ndof*j  ] -= H[i]*H[j]*dens*K(1,0)*detJ;
				ke[ndof*i+1][ndof*j+1] -= H[i]*H[j]*dens*K(1,1)*detJ;
				ke[ndof*i+1][ndof*j+2] -= H[i]*H[j]*dens*K(1,2)*detJ;

				ke[ndof*i+2][ndof*j  ] -= H[i]*H[j]*dens*K(2,0)*detJ;
				ke[ndof*i+2][ndof*j+1] -= H[i]*H[j]*dens*K(2,1)*detJ;
				ke[ndof*i+2][ndof*j+2] -= H[i]*H[j]*dens*K(2,2)*detJ;
			}
	}	
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n
void FEElasticSolidDomain::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
{
	// spatial derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate geometrical element stiffness matrix
	int neln = el.Nodes();
	int nint = el.GaussPoints();
	for (int n = 0; n<nint; ++n)
	{
		// calculate shape function gradients and jacobian
		double w = ShapeGradient(el, n, G, m_alphaf)*gw[n]*m_alphaf;

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		mat3ds& s = pt.m_s;

		for (int i = 0; i<neln; ++i)
			for (int j = i; j<neln; ++j)
			{
				double kab = (G[i]*(s * G[j]))*w;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEElasticSolidDomain::ElementMaterialStiffness(FESolidElement &el, matrix &ke)
{
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();

	// global derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// jacobian
	double detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian and shape function gradients
		detJt = ShapeGradient(el, n, G, m_alphaf)*gw[n]*m_alphaf;

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		// get the 'D' matrix
		tens4ds C = m_pMat->Tangent(mp);
		C.extract(D);

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (int i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = G[i].x;
			Gyi = G[i].y;
			Gzi = G[i].z;

			for (int j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = G[j].x;
				Gyj = G[j].y;
				Gzj = G[j].z;

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
void FEElasticSolidDomain::StiffnessMatrix(FESolver* psolver)
{
	// repeat over all solid elements
	int NE = m_Elem.size();
	
	#pragma omp parallel for shared (NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> lm;
		
		FESolidElement& el = m_Elem[iel];

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate geometrical stiffness
		ElementGeometricalStiffness(el, ke);

		// calculate material stiffness
		ElementMaterialStiffness(el, ke);

		// assign symmetic parts
		// TODO: Can this be omitted by changing the Assemble routine so that it only
		// grabs elements from the upper diagonal matrix?
		for (int i=0; i<ndof; ++i)
			for (int j=i+1; j<ndof; ++j)
				ke[j][i] = ke[i][j];

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		#pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::MassMatrix(FESolver* psolver, double scale)
{
	// element stiffness matrix
    matrix ke;
    vector<int> lm;

    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
            
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
            
        // calculate inertial stiffness
        ElementMassMatrix(el, ke, scale);
            
        // get the element's LM vector
        UnpackLM(el, lm);
            
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial()); assert(pme);

	// element stiffness matrix
    matrix ke;
    vector<int> lm;
        
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
            
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
            
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
            
        // get the element's LM vector
        UnpackLM(el, lm);
            
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FEElasticSolidDomain::ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke)
{
	FESolidElement& el = Element(iel);

	// calculate material stiffness (i.e. constitutive component)
	ElementMaterialStiffness(el, ke);

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
//! calculates element inertial stiffness matrix
void FEElasticSolidDomain::ElementMassMatrix(FESolidElement& el, matrix& ke, double a)
{
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// density
	double D = m_pMat->Density();
    
	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// shape functions
		double* H = el.H(n);

		// Jacobian
		double J0 = detJ0(el, n)*gw[n];

		for (int i=0; i<neln; ++i)
			for (int j=i; j<neln; ++j)
			{
				double kab = a*D*H[i]*H[j]*J0;
				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
    
	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	for (int i=0; i<ndof; ++i)
		for (int j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
    
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Update(const FETimeInfo& tp)
{
	// TODO: This is temporary hack for running micro-materials in parallel.
	//	     Evaluating the stress for a micro-material will make FEBio solve
	//       a new FE problem. We don't want to see the output of that problem.
	//       The logfile is a shared resource between the master FEM and the RVE
	//       in order not to corrupt the logfile we don't print anything for
	//       the RVE problem.
	// TODO: Maybe I need to create a new domain class for micro-material.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::LOG_NEVER);

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
				// reset the logfile mode
				felog.SetMode(nmode);
				berr = true;
				if (NegativeJacobian::m_boutput) e.print();
			}
		}
	}

	// reset the logfile mode
	felog.SetMode(nmode);

	// if we encountered an error, we request a running restart
	if (berr)
	{
		if (NegativeJacobian::m_boutput == false) felog.printbox("ERROR", "Negative jacobian was detected.");
		throw DoRunningRestart();
	}
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FEElasticSolidDomain::UpdateElementStress(int iel)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d r0[NELN], r[NELN], v[NELN], a[NELN];
	for (int j=0; j<neln; ++j)
	{
        FENode& node = m_pMesh->Node(el.m_node[j]);
		r0[j] = node.m_r0;
		r[j] = node.m_rt*m_alphaf + node.m_rp*(1-m_alphaf);
        v[j] = node.get_vec3d(m_dofVX, m_dofVY, m_dofVZ)*m_alphaf + node.m_vp*(1-m_alphaf);
        a[j] = node.m_at*m_alpham + node.m_ap*(1-m_alpham);
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
		pt.m_rt = el.Evaluate(r, n);

		// get the deformation gradient and determinant at intermediate time
        double Jt;
        mat3d Ft, Fp;
        Jt = defgrad(el, Ft, n);

		if (m_alphaf == 1.0)
		{
			pt.m_F = Ft;
		}
		else
		{
			defgradp(el, Fp, n);
			pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
		}

		pt.m_J = pt.m_F.det();
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (Ft - Fp)*Fi/dt;
        pt.m_v = el.Evaluate(v, n);
        pt.m_a = el.Evaluate(a, n);

        // evaluate strain energy at current time
        FEElasticMaterialPoint et = pt;
        et.m_F = Ft;
        et.m_J = Jt;

		// calculate the stress at this material point
        pt.m_s = m_pMat->Stress(mp);
        
        // adjust stress for strain energy conservation
        if (m_alphaf == 0.5) 
		{
			// evaluate strain-energy density
			pt.m_Wt = m_pMat->GetElasticMaterial()->StrainEnergyDensity(et);

            mat3ds D = pt.RateOfDeformation();
            double D2 = D.dotdot(D);
            if (D2 > 0)
                pt.m_s += D*(((pt.m_Wt-pt.m_Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEElasticSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*6);
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];

		// rigid rotational dofs
		lm[3*N + 3*i  ] = id[m_dofRU];
		lm[3*N + 3*i+1] = id[m_dofRV];
		lm[3*N + 3*i+2] = id[m_dofRW];
	}
    
    // substitute interface dofs for solid-shell interfaces
	FESolidElement& sel = static_cast<FESolidElement&>(el);
	for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofSX];
            lm[3*i+1] = id[m_dofSY];
            lm[3*i+2] = id[m_dofSZ];
        }
    }
}

//-----------------------------------------------------------------------------
// Calculate inertial forces \todo Why is F no longer needed?
void FEElasticSolidDomain::InertialForces(FEGlobalVector& R, vector<double>& F)
{
    int NE = (int)m_Elem.size();
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
        ElementInertialForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        double dens = m_pMat->Density();
        double J0 = detJ0(el, n)*gw[n];
        
        double* H = el.H(n);
        for (int i=0; i<neln; ++i)
        {
            double tmp = H[i]*J0*dens;
            fe[3*i  ] -= tmp*pt.m_a.x;
            fe[3*i+1] -= tmp*pt.m_a.y;
            fe[3*i+2] -= tmp*pt.m_a.z;
        }
    }
}

