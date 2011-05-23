#include "stdafx.h"
#include "FEShellDomain.h"
#include "FEMesh.h"

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
FENode& FEShellDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
double FEShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
{
	int neln = el.Nodes();

	double* Hrn = el.Hr(n);
	double* Hsn = el.Hs(n);
	double* Hn  = el.H(n);
	double NX, NY, NZ, MX, MY, MZ;
	double za;

	vec3d* r = el.rt();
	vec3d* D = el.Dt();

	double g = el.gt(n);

	double Ji[3][3];
	el.invjac0(Ji, n);

	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	for (int i=0; i<neln; ++i)
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
