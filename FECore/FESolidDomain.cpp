#include "stdafx.h"
#include "FESolidDomain.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
//! Domain initialization

bool FESolidDomain::Initialize(FEModel &fem)
{
	int i, j;
	FEMesh& m = *m_pMesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	m_Node.reserve(N);
	for (i=0; i<NE; ++i)
	{
		FESolidElement& e = Element(i);
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) 
			{
				tag[nj] = n++;
				m_Node.push_back(nj);
			}
		}
	}
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
// Reset data
void FESolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
FENode& FESolidDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
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
	if (D <= 0) throw NegativeJacobian(el.m_nID, n, D, &el);

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
	if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);

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
	if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);

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
void FESolidDomain::ShallowCopy(DumpStream& dmp, bool bsave)
{
	int NEL = (int) m_Elem.size();
	for (int i=0; i<NEL; ++i)
	{
		FESolidElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->ShallowCopy(dmp, bsave);
	}
}

//-----------------------------------------------------------------------------
void FESolidDomain::Serialize(DumpFile &ar)
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
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
	}
	else
	{
		ar >> m_Node;

		FEModel& fem = *ar.GetFEModel();
		int n, mat;
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FESolidElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
	}
}
