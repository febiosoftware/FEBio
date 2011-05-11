#include "stdafx.h"
#include "FESolidDomain.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
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
	vec3d x[8];
	for (i=0; i<NE; ++i)
	{
		// get the next element
		FESolidElement& e = Element(i);
		assert(e.Type() == FE_HEX);

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

