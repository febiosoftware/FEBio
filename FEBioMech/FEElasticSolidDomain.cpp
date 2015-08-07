#include "stdafx.h"
#include "FEElasticSolidDomain.h"
#include "FETransverselyIsotropic.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"

#ifdef WIN32
extern "C" int __cdecl omp_get_num_threads(void);
extern "C" int __cdecl omp_get_thread_num(void);
#else
extern "C" int omp_get_num_threads(void);
extern "C" int omp_get_thread_num(void);
#endif


//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEElasticSolidDomain::FEElasticSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_ELASTIC_SOLID_DOMAIN, pm)
{
	SetMaterial(pmat);
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
//! create a copy (overridden from FEDomain).
//! Node that this creates a copy without a material assignment
FEDomain* FEElasticSolidDomain::Copy()
{
	FEElasticSolidDomain* pd = new FEElasticSolidDomain(0, 0);
	pd->m_Elem = m_Elem;
	pd->m_Node = m_Node;
	return pd;
}

//-----------------------------------------------------------------------------
//! \todo The material point initialization needs to move to the base class.
bool FEElasticSolidDomain::Initialize(FEModel &fem)
{
	// initialize base class
	FESolidDomain::Initialize(fem);

	// get the elements material
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();

	// assign local coordinate system to each integration point
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
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
				felog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.m_nID);
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
//! Initialize element data
void FEElasticSolidDomain::InitElements()
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
void FEElasticSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
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
		detJ = detJ0(el, n)*gw[n];

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
//! Calculates element material stiffness element matrix

void FEElasticSolidDomain::ElementMaterialStiffness(FESolidElement &el, matrix &ke)
{
	int i, i3, j, j3, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
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
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

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

void FEElasticSolidDomain::ElementStiffness(FEModel& fem, int iel, matrix& ke)
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
void FEElasticSolidDomain::UpdateStresses(FEModel &fem)
{
	double dt = fem.GetCurrentStep()->m_dt;

	// TODO: This is temporary hack for running micro-materials in parallel. 
	//	     Evaluating the stress for a micro-material will make FEBio solve
	//       a new FE problem. We don't want to see the output of that problem.
	//       The logfile is a shared resource between the master FEM and the RVE
	//       in order not to corrupt the logfile we don't print anything for
	//       the RVE problem.
	// TODO: Maybe I need to create a new domain class for micro-material.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	bool berr = false;
	int NE = (int) m_Elem.size();
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			UpdateElementStress(i, dt);
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
void FEElasticSolidDomain::UpdateElementStress(int iel, double dt)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
		rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
	}

	// get the integration weights
	double* gw = el.GaussWeights();

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

		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);
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
		lm[3*i  ] = id[DOF_X];
		lm[3*i+1] = id[DOF_Y];
		lm[3*i+2] = id[DOF_Z];

		// rigid rotational dofs
		lm[3*N + 3*i  ] = id[DOF_RU];
		lm[3*N + 3*i+1] = id[DOF_RV];
		lm[3*N + 3*i+2] = id[DOF_RW];
	}
}

//-----------------------------------------------------------------------------
// Calculate inertial forces
void FEElasticSolidDomain::InertialForces(FEGlobalVector& R, vector<double>& F)
{
	double d = m_pMat->Density();

	// element mass matrix
	matrix ke;

	// element force vector
	vector<double> fe;

	// element's LM vector
	vector<int> lm;

	// loop over all elements
	int NE = Elements();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);

		int nint = el.GaussPoints();
		int neln = el.Nodes();

		ke.resize(3*neln, 3*neln);
		ke.zero();

		fe.resize(3*neln);
		
		// create the element mass matrix
		for (int n=0; n<nint; ++n)
		{
			double J0 = detJ0(el, n)*el.GaussWeights()[n];

			double* H = el.H(n);
			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					double kab = H[i]*H[j]*J0*d;
					ke[3*i  ][3*j  ] += kab;
					ke[3*i+1][3*j+1] += kab;
					ke[3*i+2][3*j+2] += kab;
				}	
		}

		// now, multiply M with F and add to R
		int* en = &el.m_node[0];
		for (int i=0; i<3*neln; ++i)
		{
			fe[i] = 0;
			for (int j=0; j<3*neln; ++j)
			{
				fe[i] -= ke[i][j]*F[3*(en[j/3]) + j%3];
			}
		}

		// get the element degrees of freedom
		UnpackLM(el, lm);

		// assemble fe into R
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
// Calculate inertial forces \todo Why is F no longer needed?
void FEElasticSolidDomain::InertialForces2(FEGlobalVector& R, vector<double>& F)
{
	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial()); assert(pme);
	double d = pme->Density();
    
    const int MN = FEElement::MAX_NODES;
    vec3d at[MN];
        
    // acceleration at integration point
    vec3d a;
        
    // loop over all elements
    int NE = Elements();
       
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
        vector<double> fe;
        vector<int> lm;
            
        FESolidElement& el = Element(iel);
            
        int nint = el.GaussPoints();
        int neln = el.Nodes();
            
        fe.assign(3*neln,0);
            
        // get the nodal accelerations
        for (int i=0; i<neln; ++i)
        {
            at[i] = m_pMesh->Node(el.m_node[i]).m_at;
        }
            
        // evaluate the element inertial force vector
        for (int n=0; n<nint; ++n)
        {
            double J0 = detJ0(el, n)*el.GaussWeights()[n];
                
            // get the acceleration for this integration point
            a = el.Evaluate(at, n);
                
            double* H = el.H(n);
            for (int i=0; i<neln; ++i)
            {
                double tmp = H[i]*J0*d;
                fe[3*i  ] -= tmp*a.x;
                fe[3*i+1] -= tmp*a.y;
                fe[3*i+2] -= tmp*a.z;
            }
        }
            
        // get the element degrees of freedom
        UnpackLM(el, lm);
            
        // assemble fe into R
        R.Assemble(el.m_node, lm, fe);
    }
}
