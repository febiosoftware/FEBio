#include "stdafx.h"
#include "FESolidDomain.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void FESolidDomain::CopyFrom(FEDomain* pd)
{
	FESolidDomain* psd = dynamic_cast<FESolidDomain*>(pd);
	m_Elem = psd->m_Elem;
}

//-----------------------------------------------------------------------------
// Reset data
void FESolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
void solve_3x3(double A[3][3], double b[3], double x[3])
{
	double D = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[1][0]*A[2][1]*A[0][2] \
			 - A[1][1]*A[2][0]*A[0][2] - A[2][2]*A[1][0]*A[0][1] - A[0][0]*A[2][1]*A[1][2];

	assert(D != 0);

	double Ai[3][3];
	Ai[0][0] = A[1][1]*A[2][2] - A[2][1]*A[1][2];
	Ai[0][1] = A[2][1]*A[0][2] - A[0][1]*A[2][2];
	Ai[0][2] = A[0][1]*A[1][2] - A[1][1]*A[0][2];

	Ai[1][0] = A[2][0]*A[1][2] - A[1][0]*A[2][2];
	Ai[1][1] = A[0][0]*A[2][2] - A[2][0]*A[0][2];
	Ai[1][2] = A[1][0]*A[0][2] - A[0][0]*A[1][2];

	Ai[2][0] = A[1][0]*A[2][1] - A[2][0]*A[1][1];
	Ai[2][1] = A[2][0]*A[0][1] - A[0][0]*A[2][1];
	Ai[2][2] = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	x[0] = (Ai[0][0]*b[0] + Ai[0][1]*b[1] + Ai[0][2]*b[2])/D;
	x[1] = (Ai[1][0]*b[0] + Ai[1][1]*b[1] + Ai[1][2]*b[2])/D;
	x[2] = (Ai[2][0]*b[0] + Ai[2][1]*b[1] + Ai[2][2]*b[2])/D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2]); 
	r[1] = b[1] - (A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2]); 
	r[2] = b[2] - (A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2]);

	double nr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
}

//-----------------------------------------------------------------------------
//! This function finds the element in which point y lies and returns
//! the isoparametric coordinates in r if an element is found
//! (This has only been implemeneted for hexes!)
FESolidElement* FESolidDomain::FindElement(vec3d y, double r[3])
{
	int i, j;
	int NE = Elements();
	vec3d x[FEElement::MAX_NODES];
	for (i=0; i<NE; ++i)
	{
		// get the next element
		FESolidElement& e = Element(i);
		assert(e.Type() == FE_HEX8G8);

		// get the element nodal coordinates
		int neln = e.Nodes();
		for (j=0; j<neln; ++j) x[j] = m_pMesh->Node(e.m_node[j]).m_rt;

		// first, as a quick check, we see if y lies in the bounding box defined by x
		FE_BOUNDING_BOX box;
		box.r0 = box.r1 = x[0];
		for (j=1; j<neln; ++j) box += x[j];

		if (box.IsInside(y))
		{
			// If the point y lies inside the box, we apply a Newton method to find
			// the isoparametric coordinates r
			r[0] = r[1] = r[2] = 0;
			const double tol = 1e-5;
			double dr[3], norm;
			double H[8], G[8][3];
			do
			{
				H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
				H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
				H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
				H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
				H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
				H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
				H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
				H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

				G[0][0] = -0.125*(1 - r[1])*(1 - r[2]); G[0][1] = -0.125*(1 - r[0])*(1 - r[2]); G[0][2] = -0.125*(1 - r[0])*(1 - r[1]);
				G[1][0] =  0.125*(1 - r[1])*(1 - r[2]); G[1][1] = -0.125*(1 + r[0])*(1 - r[2]); G[1][2] = -0.125*(1 + r[0])*(1 - r[1]);
				G[2][0] =  0.125*(1 + r[1])*(1 - r[2]); G[2][1] =  0.125*(1 + r[0])*(1 - r[2]); G[2][2] = -0.125*(1 + r[0])*(1 + r[1]);
				G[3][0] = -0.125*(1 + r[1])*(1 - r[2]); G[3][1] =  0.125*(1 - r[0])*(1 - r[2]); G[3][2] = -0.125*(1 - r[0])*(1 + r[1]);
				G[4][0] = -0.125*(1 - r[1])*(1 + r[2]); G[4][1] = -0.125*(1 - r[0])*(1 + r[2]); G[4][2] =  0.125*(1 - r[0])*(1 - r[1]);
				G[5][0] =  0.125*(1 - r[1])*(1 + r[2]); G[5][1] = -0.125*(1 + r[0])*(1 + r[2]); G[5][2] =  0.125*(1 + r[0])*(1 - r[1]);
				G[6][0] =  0.125*(1 + r[1])*(1 + r[2]); G[6][1] =  0.125*(1 + r[0])*(1 + r[2]); G[6][2] =  0.125*(1 + r[0])*(1 + r[1]);
				G[7][0] = -0.125*(1 + r[1])*(1 + r[2]); G[7][1] =  0.125*(1 - r[0])*(1 + r[2]); G[7][2] =  0.125*(1 - r[0])*(1 + r[1]);

				double R[3] = {0}, A[3][3] = {0};
				for (j=0; j<8; ++j)
				{
					R[0] += x[j].x*H[j];
					R[1] += x[j].y*H[j];
					R[2] += x[j].z*H[j];

					A[0][0] -= x[j].x*G[j][0]; A[0][1] -= x[j].x*G[j][1]; A[0][2] -= x[j].x*G[j][2];
					A[1][0] -= x[j].y*G[j][0]; A[1][1] -= x[j].y*G[j][1]; A[1][2] -= x[j].y*G[j][2];
					A[2][0] -= x[j].z*G[j][0]; A[2][1] -= x[j].z*G[j][1]; A[2][2] -= x[j].z*G[j][2];
				}
				R[0] = y.x - R[0];
				R[1] = y.y - R[1];
				R[2] = y.z - R[2];

				solve_3x3(A, R, dr);
				r[0] -= dr[0];
				r[1] -= dr[1];
				r[2] -= dr[2];

				norm = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			}
			while (norm > tol);

			// see if the point r lies inside the element
			const double eps = 1.0001;
			if ((r[0] >= -eps) && (r[0] <= eps) &&
				(r[1] >= -eps) && (r[1] <= eps) && 
				(r[2] >= -eps) && (r[2] <= eps)) return &e;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element el at integration point n.
//! The deformation gradient is returned in F and its determinant is the return
//! value of the function
double FESolidDomain::defgrad(FESolidElement &el, mat3d &F, int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// shape function derivatives
	double *Grn = el.Gr(n);
	double *Gsn = el.Gs(n);
	double *Gtn = el.Gt(n);

	// nodal points
	vec3d r[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// calculate inverse jacobian
	double Ji[3][3];
	invjac0(el, Ji, n);

	// calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	for (i=0; i<neln; ++i)
	{
		double Gri = Grn[i];
		double Gsi = Gsn[i];
		double Gti = Gtn[i];

		double x = r[i].x;
		double y = r[i].y;
		double z = r[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		double GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
		double GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
		double GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
	
		// calculate deformation gradient F
		F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
		F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
		F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
	}

	double D = F.det();
	if (D <= 0) throw NegativeJacobian(el.GetID(), n, D, &el);

	return D;
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element at point r,s,t
double FESolidDomain::defgrad(FESolidElement &el, mat3d &F, double r, double s, double t)
{
	// number of nodes
	int neln = el.Nodes();

	// shape function derivatives
	const int NME = FEElement::MAX_NODES;
	double Gr[NME], Gs[NME], Gt[NME];
	el.shape_deriv(Gr, Gs, Gt, r, s, t);

	// nodal points
	vec3d rt[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// calculate inverse jacobian
	double Ji[3][3];
	invjac0(el, Ji, r, s, t);

	// calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	for (int i=0; i<neln; ++i)
	{
		double Gri = Gr[i];
		double Gsi = Gs[i];
		double Gti = Gt[i];

		double x = rt[i].x;
		double y = rt[i].y;
		double z = rt[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		double GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
		double GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
		double GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
	
		// calculate deformation gradient F
		F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
		F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
		F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
	}

	double D = F.det();
	if (D <= 0) throw NegativeJacobian(el.GetID(), -1, D, &el);

	return D;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the reference frame at  
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjac0(FESolidElement& el, double Ji[3][3], int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;

	// calculate Jacobian
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Gri = el.Gr(n)[i];
		const double& Gsi = el.Gs(n)[i];
		const double& Gti = el.Gt(n)[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		
	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);

	// calculate the inverse jacobian
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
//! Calculate the inverse jacobian with respect to the reference frame at  
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjac0(FESolidElement& el, double Ji[3][3], double r, double s, double t)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	const int NMAX = FEElement::MAX_NODES;
	vec3d r0[NMAX];
	for (i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;

	// evaluate shape function derivatives
	double Gr[NMAX], Gs[NMAX], Gt[NMAX];
	el.shape_deriv(Gr, Gs, Gt, r, s, t);

	// calculate Jacobian
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Gri = Gr[i];
		const double& Gsi = Gs[i];
		const double& Gti = Gt[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		
	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), -1, det);

	// calculate the inverse jacobian
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
//! Calculate the inverse jacobian with respect to the current frame at  
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjact(FESolidElement& el, double Ji[3][3], int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// calculate jacobian
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Gri = el.Gr(n)[i];
		const double& Gsi = el.Gs(n)[i];
		const double& Gti = el.Gt(n)[i];
		
		const double& x = rt[i].x;
		const double& y = rt[i].y;
		const double& z = rt[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
			
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);

	// calculate inverse jacobian
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
//! Calculate the inverse jacobian with respect to the reference frame at  
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjact(FESolidElement& el, double Ji[3][3], double r, double s, double t)
{
	// number of nodes
	const int neln = el.Nodes();

	// nodal coordinates
	const int NMAX = FEElement::MAX_NODES;
	vec3d r0[NMAX];
	for (int i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// evaluate shape function derivatives
	double Gr[NMAX], Gs[NMAX], Gt[NMAX];
	el.shape_deriv(Gr, Gs, Gt, r, s, t);

	// calculate Jacobian
	double J[3][3] = {0};
	for (int i=0; i<neln; ++i)
	{
		const double& Gri = Gr[i];
		const double& Gsi = Gs[i];
		const double& Gti = Gt[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		
	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), -1, det);

	// calculate the inverse jacobian
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
//! calculate gradient of function at integration points
vec3d FESolidDomain::gradient(FESolidElement& el, double* fn, int n)
{
	double Ji[3][3];
	invjact(el, Ji, n);
				
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);

	double Gx, Gy, Gz;

	vec3d gradf;
	int N = el.Nodes();
	for (int i=0; i<N; ++i)
	{
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i]+Ji[2][0]*Gtn[i];
		Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i]+Ji[2][1]*Gtn[i];
		Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i]+Ji[2][2]*Gtn[i];

		// calculate pressure gradient
		gradf.x += Gx*fn[i];
		gradf.y += Gy*fn[i];
		gradf.z += Gz*fn[i];
	}

	return gradf;
}

//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
vec3d FESolidDomain::gradient(FESolidElement& el, vector<double>& fn, int n)
{
	double Ji[3][3];
	invjact(el, Ji, n);
				
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);

	double Gx, Gy, Gz;

	vec3d gradf;
	int N = el.Nodes();
	for (int i=0; i<N; ++i)
	{
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i]+Ji[2][0]*Gtn[i];
		Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i]+Ji[2][1]*Gtn[i];
		Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i]+Ji[2][2]*Gtn[i];

		// calculate pressure gradient
		gradf.x += Gx*fn[i];
		gradf.y += Gy*fn[i];
		gradf.z += Gz*fn[i];
	}

	return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
mat3d FESolidDomain::gradient(FESolidElement& el, vec3d* fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
    vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
    vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
    
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    mat3d gradf;
    gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gradf += fn[i] & (g1*Gr[i] + g2*Gs[i] + g3*Gt[i]);
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate material gradient of function at integration points
mat3d FESolidDomain::Gradient(FESolidElement& el, vec3d* fn, int n)
{
    double Ji[3][3];
    invjac0(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
    vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
    vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
    
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    mat3d Gradf;
    Gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        Gradf += fn[i] & (g1*Gr[i] + g2*Gs[i] + g3*Gt[i]);
    
    return Gradf;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to current frame
double FESolidDomain::detJt(FESolidElement &el, int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// shape function derivatives
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);
	
	// jacobian matrix
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Gri = Grn[i];
		const double& Gsi = Gsn[i];
		const double& Gti = Gtn[i];
		
		const double& x = rt[i].x;
		const double& y = rt[i].y;
		const double& z = rt[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

	return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FESolidDomain::detJ0(FESolidElement &el, int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;

	// shape function derivatives
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);
	
	// jacobian matrix
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Gri = Grn[i];
		const double& Gsi = Gsn[i];
		const double& Gti = Gtn[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
		J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
		J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
	}
		
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

	return det;
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant basis vectors of a solid element
//! at an integration point

void FESolidDomain::CoBaseVectors(FESolidElement& el, int j, vec3d g[3])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    double* Ht = el.Gt(j);
    
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    for (int i=0; i<n; ++i)
    {
        vec3d rt = m.Node(el.m_node[i]).m_rt;
        g[0] += rt*Hr[i];
        g[1] += rt*Hs[i];
        g[2] += rt*Ht[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the contravariant basis vectors of a solid element
//! at an integration point

void FESolidDomain::ContraBaseVectors(FESolidElement& el, int j, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors(el, j, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of covariant basis
//! vectors of a solid element at an integration point

void FESolidDomain::CoBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hrr = el.Grr(j); double* Hrs = el.Grs(j); double* Hrt = el.Grt(j);
    double* Hsr = el.Gsr(j); double* Hss = el.Gss(j); double* Hst = el.Gst(j);
    double* Htr = el.Gtr(j); double* Hts = el.Gts(j); double* Htt = el.Gtt(j);
    
    dg[0][0] = dg[0][1] = dg[0][2] = vec3d(0,0,0);  // derivatives of g[0]
    dg[1][0] = dg[1][1] = dg[1][2] = vec3d(0,0,0);  // derivatives of g[1]
    dg[2][0] = dg[2][1] = dg[2][2] = vec3d(0,0,0);  // derivatives of g[2]
    
    for (int i=0; i<n; ++i)
    {
        vec3d rt = m.Node(el.m_node[i]).m_rt;
        dg[0][0] += rt*Hrr[i]; dg[0][1] += rt*Hsr[i]; dg[0][2] += rt*Htr[i];
        dg[1][0] += rt*Hrs[i]; dg[1][1] += rt*Hss[i]; dg[1][2] += rt*Hts[i];
        dg[2][0] += rt*Hrt[i]; dg[2][1] += rt*Hst[i]; dg[2][2] += rt*Htt[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of contravariant basis
//! vectors of a solid element at an integration point

void FESolidDomain::ContraBaseVectorDerivatives(FESolidElement& el, int j, vec3d dgcnt[3][3])
{
    vec3d gcnt[3];
    vec3d dgcov[3][3];
    ContraBaseVectors(el, j, gcnt);
    CoBaseVectorDerivatives(el, j, dgcov);
    
    // derivatives of gcnt[0]
    dgcnt[0][0] = -gcnt[0]*(gcnt[0]*dgcov[0][0])-gcnt[1]*(gcnt[0]*dgcov[1][0])-gcnt[2]*(gcnt[0]*dgcov[2][0]);
    dgcnt[0][1] = -gcnt[0]*(gcnt[0]*dgcov[0][1])-gcnt[1]*(gcnt[0]*dgcov[1][1])-gcnt[2]*(gcnt[0]*dgcov[2][1]);
    dgcnt[0][2] = -gcnt[0]*(gcnt[0]*dgcov[0][2])-gcnt[1]*(gcnt[0]*dgcov[1][2])-gcnt[2]*(gcnt[0]*dgcov[2][2]);

    // derivatives of gcnt[1]
    dgcnt[1][0] = -gcnt[0]*(gcnt[1]*dgcov[0][0])-gcnt[1]*(gcnt[1]*dgcov[1][0])-gcnt[2]*(gcnt[1]*dgcov[2][0]);
    dgcnt[1][1] = -gcnt[0]*(gcnt[1]*dgcov[0][1])-gcnt[1]*(gcnt[1]*dgcov[1][1])-gcnt[2]*(gcnt[1]*dgcov[2][1]);
    dgcnt[1][2] = -gcnt[0]*(gcnt[1]*dgcov[0][2])-gcnt[1]*(gcnt[1]*dgcov[1][2])-gcnt[2]*(gcnt[1]*dgcov[2][2]);

    // derivatives of gcnt[2]
    dgcnt[2][0] = -gcnt[0]*(gcnt[2]*dgcov[0][0])-gcnt[1]*(gcnt[2]*dgcov[1][0])-gcnt[2]*(gcnt[2]*dgcov[2][0]);
    dgcnt[2][1] = -gcnt[0]*(gcnt[2]*dgcov[0][1])-gcnt[1]*(gcnt[2]*dgcov[1][1])-gcnt[2]*(gcnt[2]*dgcov[2][1]);
    dgcnt[2][2] = -gcnt[0]*(gcnt[2]*dgcov[0][2])-gcnt[1]*(gcnt[2]*dgcov[1][2])-gcnt[2]*(gcnt[2]*dgcov[2][2]);
}

//-----------------------------------------------------------------------------
//! calculate the laplacian of a vector function at an integration point
vec3d FESolidDomain::lapvec(FESolidElement& el, vec3d* fn, int n)
{
    vec3d gcnt[3];
    vec3d dgcnt[3][3];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(n); double* Hrs = el.Grs(n); double* Hrt = el.Grt(n);
    double* Hsr = el.Gsr(n); double* Hss = el.Gss(n); double* Hst = el.Gst(n);
    double* Htr = el.Gtr(n); double* Hts = el.Gts(n); double* Htt = el.Gtt(n);
    
    vec3d lapv(0,0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        lapv += fn[i]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*gcnt[0]+
                       (gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*gcnt[1]+
                       (gcnt[0]*Hrt[i]+dgcnt[0][2]*Gr[i])*gcnt[2]+
                       (gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*gcnt[0]+
                       (gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*gcnt[1]+
                       (gcnt[1]*Hst[i]+dgcnt[1][2]*Gs[i])*gcnt[2]+
                       (gcnt[2]*Htr[i]+dgcnt[2][0]*Gt[i])*gcnt[0]+
                       (gcnt[2]*Hts[i]+dgcnt[2][1]*Gt[i])*gcnt[1]+
                       (gcnt[2]*Htt[i]+dgcnt[2][2]*Gt[i])*gcnt[2]);
    
    return lapv;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of the divergence of a vector function at an integration point
vec3d FESolidDomain::gradivec(FESolidElement& el, vec3d* fn, int n)
{
    vec3d gcnt[3];
    vec3d dgcnt[3][3];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(n); double* Hrs = el.Grs(n); double* Hrt = el.Grt(n);
    double* Hsr = el.Gsr(n); double* Hss = el.Gss(n); double* Hst = el.Gst(n);
    double* Htr = el.Gtr(n); double* Hts = el.Gts(n); double* Htt = el.Gtt(n);
    
    vec3d gdv(0,0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gdv +=
        gcnt[0]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*fn[i])+
        gcnt[1]*((gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*fn[i])+
        gcnt[2]*((gcnt[0]*Hrt[i]+dgcnt[0][2]*Gr[i])*fn[i])+
        gcnt[0]*((gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*fn[i])+
        gcnt[1]*((gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*fn[i])+
        gcnt[2]*((gcnt[1]*Hst[i]+dgcnt[1][2]*Gs[i])*fn[i])+
        gcnt[0]*((gcnt[2]*Htr[i]+dgcnt[2][0]*Gt[i])*fn[i])+
        gcnt[1]*((gcnt[2]*Hts[i]+dgcnt[2][1]*Gt[i])*fn[i])+
        gcnt[2]*((gcnt[2]*Htt[i]+dgcnt[2][2]*Gt[i])*fn[i]);
    
    return gdv;
}

//-----------------------------------------------------------------------------
//! This function calculates the transpose of the spatial gradient of the shape
//! function spatial gradients, gradT(grad H), at an integration point

void FESolidDomain::gradTgradShape(FESolidElement& el, int j, vector<mat3d>& mn)
{
    int N = el.Nodes();
    vec3d g[3], dg[3][3];
    ContraBaseVectors(el, j, g);
    ContraBaseVectorDerivatives(el, j, dg);

    // get the shape functions
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    double* Ht = el.Gt(j);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(j); double* Hrs = el.Grs(j); double* Hrt = el.Grt(j);
    double* Hsr = el.Gsr(j); double* Hss = el.Gss(j); double* Hst = el.Gst(j);
    double* Htr = el.Gtr(j); double* Hts = el.Gts(j); double* Htt = el.Gtt(j);
    
    for (int i=0; i<N; ++i) {
        mn[i] = ((g[0] & dg[0][0]) + (g[1] & dg[0][1]) + (g[2] & dg[0][2]))*Hr[i]
        + ((g[0] & dg[1][0]) + (g[1] & dg[1][1]) + (g[2] & dg[1][2]))*Hs[i]
        + ((g[0] & dg[2][0]) + (g[1] & dg[2][1]) + (g[2] & dg[2][2]))*Ht[i]
        + (g[0] & g[0])*Hrr[i] + (g[1] & g[0])*Hsr[i] + (g[2] & g[0])*Htr[i]
        + (g[0] & g[1])*Hrs[i] + (g[1] & g[1])*Hss[i] + (g[2] & g[1])*Hts[i]
        + (g[0] & g[2])*Hrt[i] + (g[1] & g[2])*Hst[i] + (g[2] & g[2])*Htt[i];
    }
}

//-----------------------------------------------------------------------------
void FESolidDomain::Serialize(DumpStream &ar)
{
	if (ar.IsShallow())
	{
		int NEL = (int) m_Elem.size();
		for (int i=0; i<NEL; ++i)
		{
			FESolidElement& el = m_Elem[i];
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_Node;

			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FESolidElement& el = m_Elem[i];
				int nmat = el.GetMatID();
				ar << el.Type();
			
				ar << nmat;
				ar << el.GetID();
				ar << el.m_node;
				ar << el.m_lnode;

				for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			FEMaterial* pmat = GetMaterial();
			assert(pmat);

			ar >> m_Node;

			FEModel& fem = ar.GetFEModel();
			int n, mat, nid;
			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FESolidElement& el = m_Elem[i];
				ar >> n;

				el.SetType(n);

				ar >> mat; el.SetMatID(mat);
				ar >> nid; el.SetID(nid);
				ar >> el.m_node;
				ar >> el.m_lnode;

				for (int j=0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}
