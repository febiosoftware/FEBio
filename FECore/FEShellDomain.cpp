#include "stdafx.h"
#include "FEShellDomain.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void FEShellDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FEShellElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Init(false);
	}
}

//-----------------------------------------------------------------------------
bool FEShellDomain::Initialize(FEModel &fem)
{
	int i, j;
	FEMesh& m = *m_pMesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	for (i=0; i<NE; ++i)
	{
		FEShellElement& e = Element(i);
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) tag[nj] = n++;
		}
	}
	m_Node.reserve(n);
	for (i=0; i<N; ++i) if (tag[i] >= 0) m_Node.push_back(i);
	assert(m_Node.size() == n);
	return true;
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
FENode& FEShellDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
double FEShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
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
		r[i] = m_pMesh->Node(el.m_node[i]).m_rt;
		D[i] = m_pMesh->Node(el.m_node[i]).m_Dt;
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
//! Calculate the inverse jacobian with respect to the reference frame at 
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEShellDomain::invjac0(FEShellElement& el, double Ji[3][3], int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		D0[i] = m_pMesh->Node(el.m_node[i]).m_D0;
	}

	// calculate jacobian
	double* h0 = &el.m_h0[0];
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Hri = el.Hr(n)[i];
		const double& Hsi = el.Hs(n)[i];
		const double& Hi = el.H(n)[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		const double& dx = D0[i].x;
		const double& dy = D0[i].y;
		const double& dz = D0[i].z;
			
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
//! Calculate the inverse jacobian with respect to the current frame at 
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEShellDomain::invjact(FEShellElement& el, double Ji[3][3], int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d rt[FEElement::MAX_NODES], Dt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
		Dt[i] = m_pMesh->Node(el.m_node[i]).m_Dt;
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
//! Calculate jacobian with respect to reference frame
double FEShellDomain::detJ0(FEShellElement &el, int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		D0[i] = m_pMesh->Node(el.m_node[i]).m_D0;
	}

	// jacobian matrix
	double* h0 = &el.m_h0[0];
	double gt = el.gt(n);
	double J[3][3] = {0};
	for (i=0; i<neln; ++i)
	{
		const double& Hri = el.Hr(n)[i];
		const double& Hsi = el.Hs(n)[i];
		const double& Hi = el.H(n)[i];
		
		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;
		
		const double& dx = D0[i].x;
		const double& dy = D0[i].y;
		const double& dz = D0[i].z;
		
		double za = 0.5*gt*h0[i];
		
		J[0][0] += Hri*x + Hri*za*dx; J[0][1] += Hsi*x + Hsi*za*dx; J[0][2] += 0.5*h0[i]*Hi*dx;
		J[1][0] += Hri*y + Hri*za*dy; J[1][1] += Hsi*y + Hsi*za*dy; J[1][2] += 0.5*h0[i]*Hi*dy;
		J[2][0] += Hri*z + Hri*za*dz; J[2][1] += Hsi*z + Hsi*za*dz; J[2][2] += 0.5*h0[i]*Hi*dz;
	}
			
	// calculate the determinant
	double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

	return det;			
}

//-----------------------------------------------------------------------------
void FEShellDomain::ShallowCopy(DumpStream& dmp, bool bsave)
{
	int NEL = (int) m_Elem.size();
	for (int i=0; i<NEL; ++i)
	{
		FEShellElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->ShallowCopy(dmp, bsave);
	}
}

//-----------------------------------------------------------------------------
void FEShellDomain::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_Node;
		
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEShellElement& el = m_Elem[i];
			ar << el.Type();

			ar << el.GetMatID();
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_h0;

			for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		int n, mat;

		ar >> m_Node;

		FEModel& fem = *ar.GetFEModel();
		FEMaterial* pmat = GetMaterial();
		assert(pmat);

		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEShellElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_h0;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
				el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
	}
}
