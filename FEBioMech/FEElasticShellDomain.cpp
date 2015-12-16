#include "stdafx.h"
#include "FEElasticShellDomain.h"
#include "FEElasticMaterial.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <math.h>

//-----------------------------------------------------------------------------
FEElasticShellDomain::FEElasticShellDomain(FEModel* pfem) : FEShellDomain(&pfem->GetMesh()), FEElasticDomain(pfem)
{
	m_pMat = 0;
	m_dofU = pfem->GetDOFIndex("u");
	m_dofV = pfem->GetDOFIndex("v");
	m_dofW = pfem->GetDOFIndex("w");
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
bool FEElasticShellDomain::Initialize(FEModel& mdl)
{
	// initialize base class
	FEShellDomain::Initialize(mdl);

	// error flag (set true on error)
	bool bmerr = false;

	// set the local coordinate system for each integration point
	FECoordSysMap* pmap = m_pMat->GetElasticMaterial()->GetCoordinateSystemMap();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		// unpack element data
		FEShellElement& el = m_Elem[i];

		// set the local element coordinates
		if (pmap)
		{
			for (int n=0; n<el.GaussPoints(); ++n)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>();
				pt.m_Q = pmap->LocalElementCoord(el, n);
			}
		}
	}

	// check for initially inverted shells
	for (int i=0; i<Elements(); ++i)
	{
		FEShellElement& el = Element(i);
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
				bmerr = true;
			}
		}
	}

	return (bmerr == false);
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::Activate()
{
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

			if (node.m_bshell)
			{
				node.m_ID[m_dofU] = DOF_ACTIVE;
				node.m_ID[m_dofV] = DOF_ACTIVE;
				node.m_ID[m_dofW] = DOF_ACTIVE;
			}
		}
	}
}

//-----------------------------------------------------------------------------
double FEElasticShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
{
	int i;

	int neln = el.Nodes();

	double* Hrn = el.Hr(n);
	double* Hsn = el.Hs(n);
	double* Hn  = el.H(n);
	double NX, NY, NZ, MX, MY, MZ;
	double za;

	// current nodal coordinates and directors
	vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		FENode& ni = m_pMesh->Node(el.m_node[i]);
		r[i] = ni.m_rt;
		D[i] = el.m_D0[i] + ni.get_vec3d(m_dofU, m_dofV, m_dofW);
	}

	double g = el.gt(n);

	double Ji[3][3];
	invjac0(el, Ji, n);

	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	for (i=0; i<neln; ++i)
	{
		const double& Hri = Hrn[i];
		const double& Hsi = Hsn[i];
		const double& Hi  = Hn[i];

		const double& x = r[i].x;
		const double& y = r[i].y;
		const double& z = r[i].z;

		const double& dx = D[i].x;
		const double& dy = D[i].y;
		const double& dz = D[i].z;

		za = 0.5*g*el.m_h0[i];

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		NX = Ji[0][0]*Hri+Ji[1][0]*Hsi;
		NY = Ji[0][1]*Hri+Ji[1][1]*Hsi;
		NZ = Ji[0][2]*Hri+Ji[1][2]*Hsi;

		MX = za*Ji[0][0]*Hri + za*Ji[1][0]*Hsi + Ji[2][0]*0.5*el.m_h0[i]*Hi;
		MY = za*Ji[0][1]*Hri + za*Ji[1][1]*Hsi + Ji[2][1]*0.5*el.m_h0[i]*Hi;
		MZ = za*Ji[0][2]*Hri + za*Ji[1][2]*Hsi + Ji[2][2]*0.5*el.m_h0[i]*Hi;

		// calculate deformation gradient F
		F[0][0] += NX*x + MX*dx; F[0][1] += NY*x + MY*dx; F[0][2] += NZ*x + MZ*dx;
		F[1][0] += NX*y + MX*dy; F[1][1] += NY*y + MY*dy; F[1][2] += NZ*y + MZ*dy;
		F[2][0] += NX*z + MX*dz; F[2][1] += NY*z + MY*dz; F[2][2] += NZ*z + MZ*dz;
	}

	double V = F.det();
	if (V <= 0) throw NegativeJacobian(el.m_nID, n, V, &el);

	return V;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the current frame at 
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEElasticShellDomain::invjact(FEShellElement& el, double Ji[3][3], int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d rt[FEElement::MAX_NODES], Dt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		FENode& ni = m_pMesh->Node(el.m_node[i]);
		rt[i] = ni.m_rt;
		Dt[i] = el.m_D0[i] + ni.get_vec3d(m_dofU, m_dofV, m_dofW);
	}

	// calculate jacobian
	double* h0 = &el.m_h0[0];
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Hri = el.Hr(n)[i];
		const double& Hsi = el.Hs(n)[i];
		const double& Hi = el.H(n)[i];
		
		const double& x = rt[i].x;
		const double& y = rt[i].y;
		const double& z = rt[i].z;
		
		const double& dx = Dt[i].x;
		const double& dy = Dt[i].y;
		const double& dz = Dt[i].z;
			
		double za = 0.5*el.gt(n)*h0[i];
			
		J[0][0] += Hri*x + Hri*za*dx; J[0][1] += Hsi*x + Hsi*za*dx; J[0][2] += 0.5*h0[i]*Hi*dx;
		J[1][0] += Hri*y + Hri*za*dy; J[1][1] += Hsi*y + Hsi*za*dy; J[1][2] += 0.5*h0[i]*Hi*dy;
		J[2][0] += Hri*z + Hri*za*dz; J[2][1] += Hsi*z + Hsi*za*dz; J[2][2] += 0.5*h0[i]*Hi*dz;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);
		
	// calculate the inverse of the jacobian
	double deti = 1.0 / det;
			
	Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
	Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
	Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
	
	Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
	Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
	Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
	
	Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
	Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
	Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);

	return det;
}

//-----------------------------------------------------------------------------
// Calculates the forces due to the stress
void FEElasticShellDomain::InternalForces(FEGlobalVector& R)
{
	// element force vector
	vector<double> fe;

	vector<int> lm;

	int NS = m_Elem.size();
	for (int i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = m_Elem[i];

		// create the element force vector and initialize to zero
		int ndof = 6*el.Nodes();
		fe.assign(ndof, 0);

		// calculate element's internal force
		ElementInternalForce(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble the residual
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements
//! Note that we use a one-point gauss integration rule for the thickness
//! integration. This will integrate linear functions exactly.

void FEElasticShellDomain::ElementInternalForce(FEShellElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	const double* Gr, *Gs, *H;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();
	double gt;

	double* h0 = &el.m_h0[0];
	double za;

	double Nx, Ny, Nz;
	double Mx, My, Mz;

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		detJt = invjact(el, Ji, n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.m_s;

		gt = el.gt(n);

		Gr = el.Hr(n);
		Gs = el.Hs(n);
		H  = el.H(n);

		for (i=0; i<neln; ++i)
		{
			za = 0.5*gt*h0[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Nx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i];
			Ny = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i];
			Nz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i];

			Mx = za*Ji[0][0]*Gr[i] + za*Ji[1][0]*Gs[i] + Ji[2][0]*0.5*h0[i]*H[i];
			My = za*Ji[0][1]*Gr[i] + za*Ji[1][1]*Gs[i] + Ji[2][1]*0.5*h0[i]*H[i];
			Mz = za*Ji[0][2]*Gr[i] + za*Ji[1][2]*Gs[i] + Ji[2][2]*0.5*h0[i]*H[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[6*i  ] -= ( Nx*s.xx()  +
				           Ny*s.xy() +
					       Nz*s.xz() )*detJt;

			fe[6*i+1] -= ( Ny*s.yy()  +
				           Nx*s.xy() +
					       Nz*s.yz() )*detJt;

			fe[6*i+2] -= ( Nz*s.zz()  +
				           Ny*s.yz() +
					       Nx*s.xz() )*detJt;

			fe[6*i+3] -= ( Mx*s.xx()  +
				           My*s.xy() +
					       Mz*s.xz() )*detJt;

			fe[6*i+4] -= ( My*s.yy()  +
				           Mx*s.xy() +
					       Mz*s.yz() )*detJt;

			fe[6*i+5] -= ( Mz*s.zz()  +
				           My*s.yz() +
					       Mx*s.xz() )*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
	// element force vector
	vector<double> fe;

	vector<int> lm;

	int NS = m_Elem.size();
	for (int i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = m_Elem[i];

		// create the element force vector and initialize to zero
		int ndof = 6*el.Nodes();
		fe.assign(ndof, 0);

		// apply body forces to shells
		ElementBodyForce(BF, el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble the residual
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! Calculates element body forces for shells

void FEElasticShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
	double dens = m_pMat->Density();

	// calculate the average thickness
	double* h0 = &el.m_h0[0], gt, za;

	// integration weights
	double* gw = el.GaussWeights();
	double *Hn, detJ;

	// loop over integration points
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[4], rt[4];
	for (int i=0; i<neln; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);

		detJ = detJ0(el, n)*gw[n];
		Hn  = el.H(n);
		gt = el.gt(n);

		// get the force
		vec3d f = BF.force(mp);

		for (int i=0; i<neln; ++i)
		{
			za = 0.5*gt*h0[i];

			fe[6*i  ] -= Hn[i]*f.x*dens*detJ;
			fe[6*i+1] -= Hn[i]*f.y*dens*detJ;
			fe[6*i+2] -= Hn[i]*f.z*dens*detJ;

			fe[6*i+3] -= za*Hn[i]*dens*f.x*detJ;
			fe[6*i+4] -= za*Hn[i]*dens*f.y*detJ;
			fe[6*i+5] -= za*Hn[i]*dens*f.z*detJ;
		}
	}
}

//-----------------------------------------------------------------------------

void FEElasticShellDomain::StiffnessMatrix(FESolver* psolver)
{
	FEModel& fem = psolver->GetFEModel();

	matrix ke;

	vector<int> lm;

	int NS = m_Elem.size();
	for (int iel=0; iel<NS; ++iel)
	{
		FEShellElement& el = m_Elem[iel];

		// get the elements material
		FEMaterial* pmat = m_pMat;

		// create the element's stiffness matrix
		int ndof = 6*el.Nodes();
		ke.resize(ndof, ndof);

		// calculate the element stiffness matrix
		ElementStiffness(iel, ke);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);

		if (fem.GetCurrentStep()->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
		{
			fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*(iel)/ NS);
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the shell element stiffness matrix

void FEElasticShellDomain::ElementStiffness(int iel, matrix& ke)
{
	FEShellElement& el = Element(iel);

	int i, i6, j, j6, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 6*neln;

	// stiffness components for the initial stress component of stiffness matrix
	double kab;

	// global derivatives of shape functions
	// NOTE: hard-coding of quad elements!
	// Gx = dH/dx
	double Nx[4], Ny[4], Nz[4];
	double Mx[4], My[4], Mz[4];

	double Nxi, Nyi, Nzi;
	double Nxj, Nyj, Nzj;
	double Mxi, Myi, Mzi;
	double Mxj, Myj, Mzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][6];

	// element stress
	mat3ds s;

	double *Grn, *Gsn, *Hn;
	double Gr, Gs, H;

	// jacobian
	double Ji[3][3], detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate the average thickness
	double* h0 = &el.m_h0[0], gt, za;

	// calculate element stiffness matrix
	ke.zero();
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate jacobian
		detJt = invjact(el, Ji, n)*gw[n];

		Grn = el.Hr(n);
		Gsn = el.Hs(n);
		Hn  = el.H(n);

		gt = el.gt(n);

		// ------------ constitutive component --------------

		// setup the material point
		// NOTE: deformation gradient has already been calculated in stress routine

		// get the 'D' matrix
		tens4ds C = m_pMat->Tangent(mp);
		C.extract(D);

		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			H  = Hn[i];

			za = 0.5*gt*h0[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Nx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs;
			Ny[i] = Ji[0][1]*Gr+Ji[1][1]*Gs;
			Nz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs;

			Mx[i] = za*Ji[0][0]*Gr + za*Ji[1][0]*Gs + Ji[2][0]*0.5*h0[i]*H;
			My[i] = za*Ji[0][1]*Gr + za*Ji[1][1]*Gs + Ji[2][1]*0.5*h0[i]*H;
			Mz[i] = za*Ji[0][2]*Gr + za*Ji[1][2]*Gs + Ji[2][2]*0.5*h0[i]*H;
		}

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (i=0, i6=0; i<neln; ++i, i6 += 6)
		{
			Nxi = Nx[i];
			Nyi = Ny[i];
			Nzi = Nz[i];

			Mxi = Mx[i];
			Myi = My[i];
			Mzi = Mz[i];

			for (j=i, j6 = i6; j<neln; ++j, j6 += 6)
			{
				Nxj = Nx[j];
				Nyj = Ny[j];
				Nzj = Nz[j];

				Mxj = Mx[j];
				Myj = My[j];
				Mzj = Mz[j];


				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Nxj+D[0][3]*Nyj+D[0][5]*Nzj);
				DBL[0][1] = (D[0][1]*Nyj+D[0][3]*Nxj+D[0][4]*Nzj);
				DBL[0][2] = (D[0][2]*Nzj+D[0][4]*Nyj+D[0][5]*Nxj);
				DBL[0][3] = (D[0][0]*Mxj+D[0][3]*Myj+D[0][5]*Mzj);
				DBL[0][4] = (D[0][1]*Myj+D[0][3]*Mxj+D[0][4]*Mzj);
				DBL[0][5] = (D[0][2]*Mzj+D[0][4]*Myj+D[0][5]*Mxj);

				DBL[1][0] = (D[1][0]*Nxj+D[1][3]*Nyj+D[1][5]*Nzj);
				DBL[1][1] = (D[1][1]*Nyj+D[1][3]*Nxj+D[1][4]*Nzj);
				DBL[1][2] = (D[1][2]*Nzj+D[1][4]*Nyj+D[1][5]*Nxj);
				DBL[1][3] = (D[1][0]*Mxj+D[1][3]*Myj+D[1][5]*Mzj);
				DBL[1][4] = (D[1][1]*Myj+D[1][3]*Mxj+D[1][4]*Mzj);
				DBL[1][5] = (D[1][2]*Mzj+D[1][4]*Myj+D[1][5]*Mxj);

				DBL[2][0] = (D[2][0]*Nxj+D[2][3]*Nyj+D[2][5]*Nzj);
				DBL[2][1] = (D[2][1]*Nyj+D[2][3]*Nxj+D[2][4]*Nzj);
				DBL[2][2] = (D[2][2]*Nzj+D[2][4]*Nyj+D[2][5]*Nxj);
				DBL[2][3] = (D[2][0]*Mxj+D[2][3]*Myj+D[2][5]*Mzj);
				DBL[2][4] = (D[2][1]*Myj+D[2][3]*Mxj+D[2][4]*Mzj);
				DBL[2][5] = (D[2][2]*Mzj+D[2][4]*Myj+D[2][5]*Mxj);

				DBL[3][0] = (D[3][0]*Nxj+D[3][3]*Nyj+D[3][5]*Nzj);
				DBL[3][1] = (D[3][1]*Nyj+D[3][3]*Nxj+D[3][4]*Nzj);
				DBL[3][2] = (D[3][2]*Nzj+D[3][4]*Nyj+D[3][5]*Nxj);
				DBL[3][3] = (D[3][0]*Mxj+D[3][3]*Myj+D[3][5]*Mzj);
				DBL[3][4] = (D[3][1]*Myj+D[3][3]*Mxj+D[3][4]*Mzj);
				DBL[3][5] = (D[3][2]*Mzj+D[3][4]*Myj+D[3][5]*Mxj);

				DBL[4][0] = (D[4][0]*Nxj+D[4][3]*Nyj+D[4][5]*Nzj);
				DBL[4][1] = (D[4][1]*Nyj+D[4][3]*Nxj+D[4][4]*Nzj);
				DBL[4][2] = (D[4][2]*Nzj+D[4][4]*Nyj+D[4][5]*Nxj);
				DBL[4][3] = (D[4][0]*Mxj+D[4][3]*Myj+D[4][5]*Mzj);
				DBL[4][4] = (D[4][1]*Myj+D[4][3]*Mxj+D[4][4]*Mzj);
				DBL[4][5] = (D[4][2]*Mzj+D[4][4]*Myj+D[4][5]*Mxj);

				DBL[5][0] = (D[5][0]*Nxj+D[5][3]*Nyj+D[5][5]*Nzj);
				DBL[5][1] = (D[5][1]*Nyj+D[5][3]*Nxj+D[5][4]*Nzj);
				DBL[5][2] = (D[5][2]*Nzj+D[5][4]*Nyj+D[5][5]*Nxj);
				DBL[5][3] = (D[5][0]*Mxj+D[5][3]*Myj+D[5][5]*Mzj);
				DBL[5][4] = (D[5][1]*Myj+D[5][3]*Mxj+D[5][4]*Mzj);
				DBL[5][5] = (D[5][2]*Mzj+D[5][4]*Myj+D[5][5]*Mxj);

				ke[i6  ][j6  ] += (Nxi*DBL[0][0] + Nyi*DBL[3][0] + Nzi*DBL[5][0] )*detJt;
				ke[i6  ][j6+1] += (Nxi*DBL[0][1] + Nyi*DBL[3][1] + Nzi*DBL[5][1] )*detJt;
				ke[i6  ][j6+2] += (Nxi*DBL[0][2] + Nyi*DBL[3][2] + Nzi*DBL[5][2] )*detJt;
				ke[i6  ][j6+3] += (Nxi*DBL[0][3] + Nyi*DBL[3][3] + Nzi*DBL[5][3] )*detJt;
				ke[i6  ][j6+4] += (Nxi*DBL[0][4] + Nyi*DBL[3][4] + Nzi*DBL[5][4] )*detJt;
				ke[i6  ][j6+5] += (Nxi*DBL[0][5] + Nyi*DBL[3][5] + Nzi*DBL[5][5] )*detJt;

				ke[i6+1][j6  ] += (Nyi*DBL[1][0] + Nxi*DBL[3][0] + Nzi*DBL[4][0] )*detJt;
				ke[i6+1][j6+1] += (Nyi*DBL[1][1] + Nxi*DBL[3][1] + Nzi*DBL[4][1] )*detJt;
				ke[i6+1][j6+2] += (Nyi*DBL[1][2] + Nxi*DBL[3][2] + Nzi*DBL[4][2] )*detJt;
				ke[i6+1][j6+3] += (Nyi*DBL[1][3] + Nxi*DBL[3][3] + Nzi*DBL[4][3] )*detJt;
				ke[i6+1][j6+4] += (Nyi*DBL[1][4] + Nxi*DBL[3][4] + Nzi*DBL[4][4] )*detJt;
				ke[i6+1][j6+5] += (Nyi*DBL[1][5] + Nxi*DBL[3][5] + Nzi*DBL[4][5] )*detJt;

				ke[i6+2][j6  ] += (Nzi*DBL[2][0] + Nyi*DBL[4][0] + Nxi*DBL[5][0] )*detJt;
				ke[i6+2][j6+1] += (Nzi*DBL[2][1] + Nyi*DBL[4][1] + Nxi*DBL[5][1] )*detJt;
				ke[i6+2][j6+2] += (Nzi*DBL[2][2] + Nyi*DBL[4][2] + Nxi*DBL[5][2] )*detJt;
				ke[i6+2][j6+3] += (Nzi*DBL[2][3] + Nyi*DBL[4][3] + Nxi*DBL[5][3] )*detJt;
				ke[i6+2][j6+4] += (Nzi*DBL[2][4] + Nyi*DBL[4][4] + Nxi*DBL[5][4] )*detJt;
				ke[i6+2][j6+5] += (Nzi*DBL[2][5] + Nyi*DBL[4][5] + Nxi*DBL[5][5] )*detJt;

				ke[i6+3][j6  ] += (Mxi*DBL[0][0] + Myi*DBL[3][0] + Mzi*DBL[5][0] )*detJt;
				ke[i6+3][j6+1] += (Mxi*DBL[0][1] + Myi*DBL[3][1] + Mzi*DBL[5][1] )*detJt;
				ke[i6+3][j6+2] += (Mxi*DBL[0][2] + Myi*DBL[3][2] + Mzi*DBL[5][2] )*detJt;
				ke[i6+3][j6+3] += (Mxi*DBL[0][3] + Myi*DBL[3][3] + Mzi*DBL[5][3] )*detJt;
				ke[i6+3][j6+4] += (Mxi*DBL[0][4] + Myi*DBL[3][4] + Mzi*DBL[5][4] )*detJt;
				ke[i6+3][j6+5] += (Mxi*DBL[0][5] + Myi*DBL[3][5] + Mzi*DBL[5][5] )*detJt;

				ke[i6+4][j6  ] += (Myi*DBL[1][0] + Mxi*DBL[3][0] + Mzi*DBL[4][0] )*detJt;
				ke[i6+4][j6+1] += (Myi*DBL[1][1] + Mxi*DBL[3][1] + Mzi*DBL[4][1] )*detJt;
				ke[i6+4][j6+2] += (Myi*DBL[1][2] + Mxi*DBL[3][2] + Mzi*DBL[4][2] )*detJt;
				ke[i6+4][j6+3] += (Myi*DBL[1][3] + Mxi*DBL[3][3] + Mzi*DBL[4][3] )*detJt;
				ke[i6+4][j6+4] += (Myi*DBL[1][4] + Mxi*DBL[3][4] + Mzi*DBL[4][4] )*detJt;
				ke[i6+4][j6+5] += (Myi*DBL[1][5] + Mxi*DBL[3][5] + Mzi*DBL[4][5] )*detJt;

				ke[i6+5][j6  ] += (Mzi*DBL[2][0] + Myi*DBL[4][0] + Mxi*DBL[5][0] )*detJt;
				ke[i6+5][j6+1] += (Mzi*DBL[2][1] + Myi*DBL[4][1] + Mxi*DBL[5][1] )*detJt;
				ke[i6+5][j6+2] += (Mzi*DBL[2][2] + Myi*DBL[4][2] + Mxi*DBL[5][2] )*detJt;
				ke[i6+5][j6+3] += (Mzi*DBL[2][3] + Myi*DBL[4][3] + Mxi*DBL[5][3] )*detJt;
				ke[i6+5][j6+4] += (Mzi*DBL[2][4] + Myi*DBL[4][4] + Mxi*DBL[5][4] )*detJt;
				ke[i6+5][j6+5] += (Mzi*DBL[2][5] + Myi*DBL[4][5] + Mxi*DBL[5][5] )*detJt;
			}
		}

		// ------------ initial stress component --------------
	
		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		s = pt.m_s;

		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				// the v-v component
				kab = (Nx[i]*(s.xx()*Nx[j]+s.xy()*Ny[j]+s.xz()*Nz[j]) +
					   Ny[i]*(s.xy()*Nx[j]+s.yy()*Ny[j]+s.yz()*Nz[j]) + 
					   Nz[i]*(s.xz()*Nx[j]+s.yz()*Ny[j]+s.zz()*Nz[j]))*detJt;

				ke[6*i  ][6*j  ] += kab;
				ke[6*i+1][6*j+1] += kab;
				ke[6*i+2][6*j+2] += kab;

				// the v-t component
				kab = (Nx[i]*(s.xx()*Mx[j]+s.xy()*My[j]+s.xz()*Mz[j]) +
					   Ny[i]*(s.xy()*Mx[j]+s.yy()*My[j]+s.yz()*Mz[j]) + 
					   Nz[i]*(s.xz()*Mx[j]+s.yz()*My[j]+s.zz()*Mz[j]))*detJt;

				ke[6*i  ][6*j+3] += kab;
				ke[6*i+1][6*j+4] += kab;
				ke[6*i+2][6*j+5] += kab;

				// the t-v component
				kab = (Mx[i]*(s.xx()*Nx[j]+s.xy()*Ny[j]+s.xz()*Nz[j]) +
					   My[i]*(s.xy()*Nx[j]+s.yy()*Ny[j]+s.yz()*Nz[j]) + 
					   Mz[i]*(s.xz()*Nx[j]+s.yz()*Ny[j]+s.zz()*Nz[j]))*detJt;

				ke[6*i+3][6*j  ] += kab;
				ke[6*i+4][6*j+1] += kab;
				ke[6*i+5][6*j+2] += kab;

				// the t-t component
				kab = (Mx[i]*(s.xx()*Mx[j]+s.xy()*My[j]+s.xz()*Mz[j]) +
					   My[i]*(s.xy()*Mx[j]+s.yy()*My[j]+s.yz()*Mz[j]) + 
					   Mz[i]*(s.xz()*Mx[j]+s.yz()*My[j]+s.zz()*Mz[j]))*detJt;

				ke[6*i+3][6*j+3] += kab;
				ke[6*i+4][6*j+4] += kab;
				ke[6*i+5][6*j+5] += kab;
			}

	} // end loop over gauss-points

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}


//-----------------------------------------------------------------------------
//! Calculates body forces for shells

void FEElasticShellDomain::ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe)
{
	int NF = fem.BodyLoads();
	for (int nf = 0; nf < NF; ++nf)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(nf));
		if (pbf)
		{
			double dens = m_pMat->Density();

			// calculate the average thickness
			double* h0 = &el.m_h0[0], gt, za;

			// integration weights
			double* gw = el.GaussWeights();
			double *Hn, detJ;

			// loop over integration points
			int nint = el.GaussPoints();
			int neln = el.Nodes();

			// nodal coordinates
			vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
			for (int i=0; i<neln; ++i)
			{
				r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
				rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
			}

			for (int n=0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				pt.m_r0 = el.Evaluate(r0, n);
				pt.m_rt = el.Evaluate(rt, n);

				detJ = detJ0(el, n)*gw[n];
				Hn  = el.H(n);
				gt = el.gt(n);

				// get the force
				vec3d f = pbf->force(mp);

				for (int i=0; i<neln; ++i)
				{
					za = 0.5*gt*h0[i];

					fe[6*i  ] -= Hn[i]*f.x*dens*detJ;
					fe[6*i+1] -= Hn[i]*f.y*dens*detJ;
					fe[6*i+2] -= Hn[i]*f.z*dens*detJ;

					fe[6*i+3] -= za*Hn[i]*dens*f.x*detJ;
					fe[6*i+4] -= za*Hn[i]*dens*f.y*detJ;
					fe[6*i+5] -= za*Hn[i]*dens*f.z*detJ;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::UpdateStresses(FEModel &fem)
{
	FEMesh& mesh = *GetMesh();
	vec3d r0[9], rt[9];

	int n;
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FEShellElement& el = m_Elem[i];

		// get the number of integration points
		int nint = el.GaussPoints();

		// number of nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (int j=0; j<neln; ++j)
		{
			r0[j] = mesh.Node(el.m_node[j]).m_r0;
			rt[j] = mesh.Node(el.m_node[j]).m_rt;
		}

		// get the integration weights
		double* gw = el.GaussWeights();

		// loop over the integration points and calculate
		// the stress at the integration point
		for (n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
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
}


//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure
//! Note that for the shell elements the lm order is different compared
//! to the solid element ordering. This is because for shell elements the
//! nodes have six degrees of freedom each, where for solids they only
//! have 3 dofs.
void FEElasticShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*9);
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[6*i  ] = id[m_dofX];
		lm[6*i+1] = id[m_dofY];
		lm[6*i+2] = id[m_dofZ];

		// next the rotational dofs
		lm[6*i+3] = id[m_dofU];
		lm[6*i+4] = id[m_dofV];
		lm[6*i+5] = id[m_dofW];

		// rigid rotational dofs
		lm[6*N + 3*i  ] = id[m_dofRU];
		lm[6*N + 3*i+1] = id[m_dofRV];
		lm[6*N + 3*i+2] = id[m_dofRW];
	}
}
