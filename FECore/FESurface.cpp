// FESurface.cpp: implementation of the FESurface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
//! Initialize surface node data structure
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FESurface::Init()
{
	// make sure that there is a surface defined
	if (Elements() == 0) return false;

	// get the mesh to which this surface belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	vector<int> tag; tag.assign(mesh.Nodes(), -1);

	// let's find all nodes the surface needs
	int nn = 0;
	int ne = Elements();
	for (int i=0; i<ne; ++i)
	{
		FESurfaceElement& el = Element(i);
		el.m_lid = i;

		for (int j=0; j<el.Nodes(); ++j)
		{
			// get the global node number
			int m = el.m_node[j];
		
			// create a local node number
			if (tag[m] == -1) tag[m] = nn++;

			// set the local node number
			el.m_lnode[j] = tag[m];
		}
	}

	// allocate node index table
	m_node.resize(nn);

	// fill the node index table
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		if (tag[i] >= 0)
		{
			m_node[tag[i]] = i;
		}
	}

	// create the node element list
	m_NEL.Create(*this);
	m_NET.Create(this, 1);

	// see if we can find all elements that the faces belong to
	for (int i=0; i<ne; ++i)
	{
		FESurfaceElement& el = Element(i);
		if (el.m_nelem < 0) el.m_nelem = FindElement(el);
		assert(el.m_nelem >= 0);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Find the element that a face belongs to
//!
int FESurface::FindElement(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *GetMesh();
	FENodeElemList& NEL = mesh.NodeElementList();

	vector<int>& sf = el.m_node;
	int node = el.m_node[0];
	int nval = NEL.Valence(node);
	FEElement** ppe = NEL.ElementList(node);
	for (int i=0; i<nval; ++i)
	{
		FEElement& e = *ppe[i];
		int nfaces = mesh.Faces(e);
		
		int nf[8], nn;
		for (int j=0; j<nfaces; ++j)
		{
			nn = mesh.GetFace(e, j, nf);
			if ((nn == 3) && (el.Nodes() == 3))
			{
				if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return e.m_nID;
			}
			else if ((nn == 4) && (el.Nodes() == 4))
			{
				if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) return e.m_nID;
			}
			else if ((nn == 6) && (el.Nodes() == 6))
			{
				if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return e.m_nID;
			}
		}
	}

	return -1;
}

//-----------------------------------------------------------------------------
void FESurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*MAX_NDOFS);

	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		int* id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[DOF_X];
		lm[3*i+1] = id[DOF_Y];
		lm[3*i+2] = id[DOF_Z];

		// now the pressure dofs
		lm[3*N+i] = id[DOF_P];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[DOF_RU];
		lm[4*N + 3*i+1] = id[DOF_RV];
		lm[4*N + 3*i+2] = id[DOF_RW];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;

		lm[10*N + i] = id[DOF_T];
		
		// concentration dofs
		for (int k=0; k<MAX_CDOFS; ++k)
			lm[(11+k)*N + i] = id[DOF_C+k];
	}
}

//-----------------------------------------------------------------------------
// project onto a triangular face
vec3d project2tri(vec3d* y, vec3d x, double& r, double& s)
{
	// calculate base vectors 
	vec3d e1 = y[1] - y[0];
	vec3d e2 = y[2] - y[0];

	// calculate plane normal
	vec3d n = e1^e2; n.unit();

	// project x onto the plane
	vec3d q = x - n*((x-y[0])*n);

	// set up metric tensor
	double G[2][2];
	G[0][0] = e1*e1;
	G[0][1] = G[1][0] = e1*e2;
	G[1][1] = e2*e2;

	// invert metric tensor
	double D = G[0][0]*G[1][1] - G[0][1]*G[1][0];
	double Gi[2][2];
	Gi[0][0] = G[1][1]/D;
	Gi[1][1] = G[0][0]/D;
	Gi[0][1] = Gi[1][0] = -G[0][1]/D;

	// calculate dual base vectors
	vec3d E1 = e1*Gi[0][0] + e2*Gi[0][1];
	vec3d E2 = e1*Gi[1][0] + e2*Gi[1][1];

	// now we can calculate r and s
	vec3d t = q - y[0];
	r = t*E1;
	s = t*E2;

	return q;
}

//-----------------------------------------------------------------------------
// project onto a quadrilateral surface.
bool project2quad(vec3d* y, vec3d x, double& r, double& s, vec3d& q)
{
	double R[2], u[2], D;
	double gr[4] = {-1, +1, +1, -1};
	double gs[4] = {-1, -1, +1, +1};
	double H[4], Hr[4], Hs[4], Hrs[4];

	int i, j;
	int NMAX = 50, n=0;

	// evaulate scalar products
	double xy[4] = {x*y[0], x*y[1], x*y[2], x*y[3]};
	double yy[4][4];
	yy[0][0] = y[0]*y[0]; yy[1][1] = y[1]*y[1]; yy[2][2] = y[2]*y[2]; yy[3][3] = y[3]*y[3];
	yy[0][1] = yy[1][0] = y[0]*y[1];
	yy[0][2] = yy[2][0] = y[0]*y[2];
	yy[0][3] = yy[3][0] = y[0]*y[3];
	yy[1][2] = yy[2][1] = y[1]*y[2];
	yy[1][3] = yy[3][1] = y[1]*y[3];
	yy[2][3] = yy[3][2] = y[2]*y[3];

	// loop until converged
	bool bconv = false;
	double normu;
	do
	{
		// evaluate shape functions and shape function derivatives.
		for (i=0; i<4; ++i)
		{
			H[i] = 0.25*(1+gr[i]*r)*(1+gs[i]*s);
	
			Hr[i] = 0.25*gr[i]*( 1 + gs[i]*s );
			Hs[i] = 0.25*gs[i]*( 1 + gr[i]*r );

			Hrs[i] = 0.25*gr[i]*gs[i];
		}

		// set up the system of equations
		R[0] = R[1] = 0;
		double A[2][2] = {0};
		for (i=0; i<4; ++i)
		{
			R[0] -= (xy[i])*Hr[i];
			R[1] -= (xy[i])*Hs[i];

			A[0][1] += (xy[i])*Hrs[i];
			A[1][0] += (xy[i])*Hrs[i];

			for (j=0; j<4; ++j)
			{
				double yij = yy[i][j];
				R[0] -= -H[j]*Hr[i]*(yij);
				R[1] -= -H[j]*Hs[i]*(yij);

				A[0][0] -= (yij)*(Hr[i]*Hr[j]);
				A[1][1] -= (yij)*(Hs[i]*Hs[j]);

				A[0][1] -= (yij)*(Hr[i]*Hs[j]+Hrs[i]*H[j]);
				A[1][0] -= (yij)*(Hs[i]*Hr[j]+Hrs[i]*H[j]);
			}
		}
	
		// determinant of A
		D = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// solve for u = A^(-1)*R
		u[0] = (A[1][1]*R[0] - A[0][1]*R[1])/D;
		u[1] = (A[0][0]*R[1] - A[1][0]*R[0])/D;

		// calculate displacement norm
		normu = u[0]*u[0]+u[1]*u[1];

		// check for convergence
		bconv = ((normu < 1e-10));
		if (!bconv && (n <= NMAX))
		{
			// Don't update if converged otherwise the point q
			// does not correspond with the current values for (r,s)
			r += u[0];
			s += u[1];
			++n;
		}
		else break;
	}
	while (1);

	// evaluate q
	q = y[0]*H[0] + y[1]*H[1] + y[2]*H[2] + y[3]*H[3];

	return bconv;
}

//-----------------------------------------------------------------------------
// project onto a 6-node quadratic triangular element
bool project2tri6(vec3d* y, vec3d x, double& r, double& s, vec3d& q)
{
	double R[2], u[2], D;
	double H[6], Hr[6], Hs[6], Hrs[6];
	
	int i, j;
	int NMAX = 50, n=0;
	
	// evaulate scalar products
	double xy[6] = {x*y[0], x*y[1], x*y[2], x*y[3], x*y[4], x*y[5]};
	double yy[6][6];
	yy[0][0] = y[0]*y[0]; yy[1][1] = y[1]*y[1]; yy[2][2] = y[2]*y[2]; yy[3][3] = y[3]*y[3]; yy[4][4] = y[4]*y[4]; yy[5][5] = y[5]*y[5];
	yy[0][1] = yy[1][0] = y[0]*y[1];
	yy[0][2] = yy[2][0] = y[0]*y[2];
	yy[0][3] = yy[3][0] = y[0]*y[3];
	yy[0][4] = yy[4][0] = y[0]*y[4];
	yy[0][5] = yy[5][0] = y[0]*y[5];
	yy[1][2] = yy[2][1] = y[1]*y[2];
	yy[1][3] = yy[3][1] = y[1]*y[3];
	yy[1][4] = yy[4][1] = y[1]*y[4];
	yy[1][5] = yy[5][1] = y[1]*y[5];
	yy[2][3] = yy[3][2] = y[2]*y[3];
	yy[2][4] = yy[4][2] = y[2]*y[4];
	yy[2][5] = yy[5][2] = y[2]*y[5];
	yy[3][4] = yy[4][3] = y[3]*y[4];
	yy[3][5] = yy[5][3] = y[3]*y[5];
	yy[4][5] = yy[5][4] = y[4]*y[5];
	
	// loop until converged
	bool bconv = false;
	double normu;
	do
	{
		// evaluate shape functions and shape function derivatives.
		double r1 = 1.0 - r - s;
		double r2 = r;
		double r3 = s;

		H[0] = r1*(2.0*r1 - 1.0);
		H[1] = r2*(2.0*r2 - 1.0);
		H[2] = r3*(2.0*r3 - 1.0);
		H[3] = 4.0*r1*r2;
		H[4] = 4.0*r2*r3;
		H[5] = 4.0*r3*r1;

		Hr[0] = -3.0 + 4.0*r + 4.0*s;
		Hr[1] =  4.0*r - 1.0;
		Hr[2] =  0.0;
		Hr[3] =  4.0 - 8.0*r - 4.0*s;
		Hr[4] =  4.0*s;
		Hr[5] = -4.0*s;

		Hs[0] = -3.0 + 4.0*s + 4.0*r;
		Hs[1] =  0.0;
		Hs[2] =  4.0*s - 1.0;
		Hs[3] = -4.0*r;
		Hs[4] =  4.0*r;
		Hs[5] =  4.0 - 8.0*s - 4.0*r;

		Hrs[0] =  4.0;
		Hrs[1] =  0.0;
		Hrs[2] =  0.0;
		Hrs[3] = -4.0;
		Hrs[4] =  4.0;
		Hrs[5] = -4.0;

		// set up the system of equations
		R[0] = R[1] = 0;
		double A[2][2] = {0};
		for (i=0; i<6; ++i)
		{
			R[0] -= (xy[i])*Hr[i];
			R[1] -= (xy[i])*Hs[i];
			
			A[0][1] += (xy[i])*Hrs[i];
			A[1][0] += (xy[i])*Hrs[i];
			
			for (j=0; j<6; ++j)
			{
				double yij = yy[i][j];
				R[0] -= -H[j]*Hr[i]*(yij);
				R[1] -= -H[j]*Hs[i]*(yij);
				
				A[0][0] -= (yij)*(Hr[i]*Hr[j]);
				A[1][1] -= (yij)*(Hs[i]*Hs[j]);
				
				A[0][1] -= (yij)*(Hr[i]*Hs[j]+Hrs[i]*H[j]);
				A[1][0] -= (yij)*(Hs[i]*Hr[j]+Hrs[i]*H[j]);
			}
		}
		
		// determinant of A
		D = A[0][0]*A[1][1] - A[0][1]*A[1][0];
		
		// solve for u = A^(-1)*R
		u[0] = (A[1][1]*R[0] - A[0][1]*R[1])/D;
		u[1] = (A[0][0]*R[1] - A[1][0]*R[0])/D;
		
		// calculate displacement norm
		normu = u[0]*u[0]+u[1]*u[1];
		
		// check for convergence
		bconv = ((normu < 1e-10));
		if (!bconv && (n <= NMAX))
		{
			// Don't update if converged otherwise the point q
			// does not correspond with the current values for (r,s)
			r += u[0];
			s += u[1];
			++n;
		}
		else break;
	}
	while (1);
	
	// evaluate q
	q = y[0]*H[0] + y[1]*H[1] + y[2]*H[2] + y[3]*H[3]+ y[4]*H[4]+ y[5]*H[5];
	
	return bconv;
}

//-----------------------------------------------------------------------------
//! This function calculates the projection of x on the surface element el.
//! It does this by finding the solution of the nonlinear equation (x-y)*y,[a]=0,
//! where the comma denotes differentation and a ranges from 1 to 2.
//! The system is solved using the Newton-Raphson method.
//! The surface element may be either a quad or a triangular element.

vec3d FESurface::ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// number of element nodes
	int ne = el.Nodes();

	// get the elements nodal positions
	vec3d y[FEElement::MAX_NODES];
	for (int i=0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;

	// calculate normal projection of x onto element
	vec3d q;
	switch (ne)
	{
	case 3: q = project2tri(y, x, r, s); break;
	case 4: 
		// see if we get lucky
		if (project2quad(y, x, r, s, q) == false)
		{
			// the direct projection failed, so we'll try it more incrementally
			vec3d x0 = (y[0]+y[1]+y[2]+y[3])*0.25;
			r = s = 0;
			bool b = project2quad(y, x0, r, s, q);
			assert(b);
			
			double w = 0.5;
			int l = 1, N = 0, NMAX = 20;
			do
			{
				vec3d xi = x0*(1.0 - w) + x*w;
				b = project2quad(y, xi, r, s, q);
				if (b)
				{
					--l;
					if (l == 0) { x0 = xi; w = 1.0; }
					else w *= 2.0;
				}
				else 
				{
					++l;
					w *= 0.5;
				}
				++N;
			}
			while ((l >= 0) && (N<=NMAX) && (w>0.1));
		}
		break;
	case 6: 
		if (project2tri6(y, x, r, s, q)==false)
		{
//			assert(false);
		}
		break;
	default:
		assert(false);
	}

	return q;
}


//-----------------------------------------------------------------------------
//! This function calculates the projection of x on the surface element el.
//! It does this by finding the solution of the nonlinear equation (x-y)*y,[a]=0,
//! where the comma denotes differentation and a ranges from 1 to 2.
//! The system is solved using the Newton-Raphson method.
//! The surface element may be either a quad or a triangular element.

vec3d FESurface::ProjectToReferenceSurface(FESurfaceElement& el, vec3d x, double& r, double& s)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// number of element nodes
	int ne = el.Nodes();

	// get the elements nodal positions
	vec3d y[FEElement::MAX_NODES];
	for (int i=0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_r0;

	// calculate normal projection of x onto element
	vec3d q;
	switch (ne)
	{
	case 3: q = project2tri(y, x, r, s); break;
	case 4: 
		// see if we get lucky
		if (project2quad(y, x, r, s, q) == false)
		{
			// the direct projection failed, so we'll try it more incrementally
			vec3d x0 = (y[0]+y[1]+y[2]+y[3])*0.25;
			r = s = 0;
			bool b = project2quad(y, x0, r, s, q);
			assert(b);
			
			double w = 0.5;
			int l = 1, N = 0, NMAX = 20;
			do
			{
				vec3d xi = x0*(1.0 - w) + x*w;
				b = project2quad(y, xi, r, s, q);
				if (b)
				{
					--l;
					if (l == 0) { x0 = xi; w = 1.0; }
					else w *= 2.0;
				}
				else 
				{
					++l;
					w *= 0.5;
				}
				++N;
			}
			while ((l >= 0) && (N<=NMAX) && (w>0.1));
		}
		break;
	case 6: 
		if (project2tri6(y, x, r, s, q)==false)
		{
//			assert(false);
		}
		break;
	default:
		assert(false);
	}

	return q;
}

//-----------------------------------------------------------------------------
//! This function calculates the area of a surface element

double FESurface::FaceArea(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// get the number of nodes
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// get the initial nodes
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	// get the integration weights
	double* w = el.GaussWeights();

	double *Gr, *Gs;
	vec3d dxr, dxs;

	double detJ;

	double area = 0;

	int n, k;

	for (n=0; n<nint; ++n)
	{
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate jacobian
		dxr = dxs = vec3d(0,0,0);
		for (k=0; k<neln; ++k)
		{
			dxr.x += Gr[k]*r0[k].x;
			dxr.y += Gr[k]*r0[k].y;
			dxr.z += Gr[k]*r0[k].z;

			dxs.x += Gs[k]*r0[k].x;
			dxs.y += Gs[k]*r0[k].y;
			dxs.z += Gs[k]*r0[k].z;
		}

		detJ = (dxr ^ dxs).norm();

		area += w[n]*detJ;
	}

	return area;
}

//-----------------------------------------------------------------------------
//! Calculate the max element size.
double FESurface::MaxElementSize()
{
	FEMesh& mesh = *m_pMesh;
	double h2 = 0.0;
	int NS = Elements();
	for (int m=0; m<NS; ++m)
	{
		FESurfaceElement& e = Element(m);
		int n = e.Nodes();
		for (int i=0; i<n; ++i)
			for (int j=i+1; j<n; ++j)
			{
				vec3d& a = mesh.Node(e.m_node[i]).m_rt;
				vec3d& b = mesh.Node(e.m_node[j]).m_rt;
				double L2 = (b - a)*(b - a);
				if (L2 > h2) h2 = L2;
			}
	}
	return sqrt(h2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at the point with surface coordinates (r,s)
// TODO: perhaps I should place this function in the element class

mat2d FESurface::Metric0(FESurfaceElement& el, double r, double s)
{
	// nr of element nodes
	int neln = el.Nodes();
	
	// element nodes
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	
	// shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	
	// get the shape function values at this slave node
	el.shape_deriv(Hr, Hs, r, s);
	
	// get the tangent vectors
	vec3d t1(0,0,0);
	vec3d t2(0,0,0);
	for (int k=0; k<neln; ++k)
	{
		t1 += r0[k]*Hr[k];
		t2 += r0[k]*Hs[k];
	}
	
	// calculate metric tensor
	return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at the point with surface coordinates (r,s)
// TODO: perhaps I should place this function in the element class

mat2d FESurface::Metric(FESurfaceElement& el, double r, double s)
{
	// nr of element nodes
	int neln = el.Nodes();
	
	// element nodes
	vec3d rt[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	
	// shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	
	// get the shape function values at this slave node
	el.shape_deriv(Hr, Hs, r, s);
	
	// get the tangent vectors
	vec3d t1(0,0,0);
	vec3d t2(0,0,0);
	for (int k=0; k<neln; ++k)
	{
		t1 += rt[k]*Hr[k];
		t2 += rt[k]*Hs[k];
	}
	
	// calculate metric tensor
	return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Given an element an the natural coordinates of a point in this element, this
//! function returns the global position vector.
vec3d FESurface::Local2Global(FESurfaceElement &el, double r, double s)
{
	// get the mesh
	FEMesh& mesh = *m_pMesh;

	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
	for (int l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).m_rt;

	// calculate the element position
	return el.eval(y, r, s);
}

//-----------------------------------------------------------------------------
//! This function calculates the global location of an integration point
//!

vec3d FESurface::Local2Global(FESurfaceElement &el, int n)
{
	FEMesh& m = *m_pMesh;

	// get the shape functions at this integration point
	double* H = el.H(n);

	// calculate the location
	vec3d r;
	int ne = el.Nodes();
	for (int i=0; i<ne; ++i) r += m.Node(el.m_node[i]).m_rt*H[i];

	return r;
}

//-----------------------------------------------------------------------------
//! This function calculates the noraml of a surface element at integration
//! point n

vec3d FESurface::SurfaceNormal(FESurfaceElement &el, int n)
{
	int i;
	FEMesh& m = *m_pMesh;

	// get the shape function derivatives at this integration point
	double* Hr = el.Gr(n);
	double* Hs = el.Gs(n);

	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
	for (i=0; i<ne; ++i) y[i] = m.Node(el.m_node[i]).m_rt;

	// calculate the tangents
	vec3d xr, xs;
	for (i=0; i<ne; ++i)
	{
		xr += y[i]*Hr[i];
		xs += y[i]*Hs[i];
	}

	// calculate the normal
	vec3d np = xr ^ xs;
	np.unit();

	return np;
}

//-----------------------------------------------------------------------------
//! This function calculates the normal of a surface element at the natural
//! coordinates (r,s)

vec3d FESurface::SurfaceNormal(FESurfaceElement &el, double r, double s)
{
	int l;
	FEMesh& mesh = *m_pMesh;
	
	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
	for (l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).m_rt;
	
	// set up shape functions and derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	el.shape_deriv(Hr, Hs, r, s);
	
	// calculate the element tangents
	vec3d xr(0,0,0), xs;
	for (l=0; l<ne; ++l)
	{
		xr += y[l]*Hr[l];
		xs += y[l]*Hs[l];
	}
	
	// calculate the normal
	vec3d np = xr ^ xs;
	np.unit();
	
	return np;
}

//-----------------------------------------------------------------------------
//! This function checks whether the point with natural coordinates (r,s) is
//! inside the element, within a tolerance of tol.

bool FESurface::IsInsideElement(FESurfaceElement& el, double r, double s, double tol)
{
	int ne = el.Nodes();
	if (ne == 4)
	{
		// check quads
		if ((r >= -1-tol) && (r <= 1+tol) && (s >= -1-tol) && (s <= 1+tol)) return true;
		else return false;
	}
	else
	{
		// check triangles
		if ((r >= -tol) && (s >= -tol) && (r+s <= 1+tol)) return true;
		else return false;
	}
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at an integration point

void FESurface::CoBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	FEMesh& m = *m_pMesh;

	// get the nr of nodes
	int n = el.Nodes();

	// get the shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	el.shape_deriv(Hr, Hs, r, s);

	t[0] = t[1] = vec3d(0,0,0);
	for (int i=0; i<n; ++i)
	{
		t[0] += m.Node(el.m_node[i]).m_rt*Hr[i];
		t[1] += m.Node(el.m_node[i]).m_rt*Hs[i];
	}
}


//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at an integration point

void FESurface::CoBaseVectors(FESurfaceElement& el, int j, vec3d t[2])
{
	FEMesh& m = *m_pMesh;

	// get the nr of nodes
	int n = el.Nodes();

	// get the shape function derivatives
	double* Hr = el.Gr(j);
	double* Hs = el.Gs(j);

	t[0] = t[1] = vec3d(0,0,0);
	for (int i=0; i<n; ++i)
	{
		t[0] += m.Node(el.m_node[i]).m_rt*Hr[i];
		t[1] += m.Node(el.m_node[i]).m_rt*Hs[i];
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at the natural coordinates (r,s)

void FESurface::CoBaseVectors0(FESurfaceElement &el, double r, double s, vec3d t[2])
{
	int i;
	const int MN = FEElement::MAX_NODES;
	vec3d y[MN];
	double H0[MN], H1[MN];
	int n = el.Nodes();
	for (i=0; i<n; ++i) y[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	el.shape_deriv(H0, H1, r, s);
	t[0] = t[1] = vec3d(0,0,0);
	for (i=0; i<n; ++i) 
	{
		t[0] += y[i]*H0[i];
		t[1] += y[i]*H1[i];
	}
}

//-----------------------------------------------------------------------------
void FESurface::ContraBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	vec3d e[2];
	CoBaseVectors(el, r, s, e);
	mat2d M = Metric(el, r, s);
	mat2d Mi = M.inverse();

	t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
	t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}

//-----------------------------------------------------------------------------

void FESurface::ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	vec3d e[2];
	CoBaseVectors0(el, r, s, e);
	mat2d M = Metric0(el, r, s);
	mat2d Mi = M.inverse();

	t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
	t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}


//-----------------------------------------------------------------------------
// This function calculates the intersection of a ray with a triangle
// and returns true if the ray intersects the triangle.
//
bool FESurface::IntersectTri(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	vec3d e[2], E[2];
	mat2d G;

	// create base vectors on triangle
	e[0] = y[1]-y[0];
	e[1] = y[2]-y[0];

	// create triangle normal
	vec3d m = e[0]^e[1]; m.unit();

	double d = n*m;
//	if (d != 0)
	// normals should be pointing opposite to each other for valid contact
	if (d < 0)
	{
		// distance from r to plane of triangle
		g = m*(y[0] - r)/d;

		// intersection point with plane of triangle
		vec3d q = r + n*g;

		// next, we decompose q into its components
		// in the triangle basis
		// we need to create the dual basis
		// first, we calculate the metric tensor
		G[0][0] = e[0]*e[0]; G[0][1] = e[0]*e[1];
		G[1][0] = e[1]*e[0]; G[1][1] = e[1]*e[1];

		// and its inverse
		mat2d Gi = G.inverse();

		// build dual basis
		E[0] = e[0]*Gi[0][0] + e[1]*Gi[0][1];
		E[1] = e[0]*Gi[1][0] + e[1]*Gi[1][1];

		// get the components
		rs[0] = E[0]*(q - y[0]);
		rs[1] = E[1]*(q - y[0]);

		// see if the intersection point is inside the triangle
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a quad
//! and returns true if the ray intersected.
//!
bool FESurface::IntersectQuad(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the two subtriangles
	vec3d x1[3], x2[3];
	x1[0] = y[0]; x2[0] = y[2];
	x1[1] = y[1]; x2[1] = y[3];
	x1[2] = y[3]; x2[2] = y[1];

	bool b = false;
	double rp, sp;

	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = -1.0 + 2.0*rs[0];
		sp = -1.0 + 2.0*rs[1];
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 1.0 - 2.0*rs[0];
		sp = 1.0 - 2.0*rs[1];
	}

	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[4], H1[4], H2[4];

		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 5;
		do
		{
			// shape functions of quad
			H[0] = 0.25*(1 - l1)*(1 - l2);
			H[1] = 0.25*(1 + l1)*(1 - l2);
			H[2] = 0.25*(1 + l1)*(1 + l2);
			H[3] = 0.25*(1 - l1)*(1 + l2);

			// shape function derivatives
			H1[0] = -0.25*(1 - l2); H2[0] = -0.25*(1 - l1);
			H1[1] =  0.25*(1 - l2); H2[1] = -0.25*(1 + l1);
			H1[2] =  0.25*(1 + l2); H2[2] =  0.25*(1 + l1);
			H1[3] = -0.25*(1 + l2); H2[3] =  0.25*(1 - l1);

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3];

			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3];
			F3 = n;

			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

			// calculate solution increment
			dx = -(A.inverse()*F);

			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;

			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));

		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;

		// see if the point is inside the quad
		if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
			(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a surface element.
//! It simply calls the tri or quad intersection function based on the type
//! of element.
//!
bool FESurface::Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	int N = el.Nodes();

	// get the element nodes
	FEMesh& mesh = *m_pMesh;
	vec3d y[FEElement::MAX_NODES];
	for (int i=0; i<N; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;

	// call the correct intersection function
	switch (N)
	{
	case 3: return IntersectTri(y, r, n, rs, g, eps); break;
	case 4: return IntersectQuad(y, r, n, rs, g, eps); break;
	case 6: return IntersectTri6(y, r, n, rs, g, eps); break;
	default:
		assert(false);
	}

	// if we get here, the ray did not intersect the element
	return false;
}


//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a 6-node triangle
//! and returns true if the ray intersected.
//!
bool FESurface::IntersectTri6(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the four subtriangles
	vec3d x1[3], x2[3], x3[3], x4[3];
	x1[0] = y[0]; x2[0] = y[3]; x3[0] = y[5]; x4[0] = y[4];
	x1[1] = y[3]; x2[1] = y[1]; x3[1] = y[4]; x4[1] = y[5];
	x1[2] = y[5]; x2[2] = y[4]; x3[2] = y[2]; x4[2] = y[3];
	
	bool b = false;
	double rp, sp;
	
	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 0.5 + rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x3, r, n, rs, g, eps))
	{
		// we've intersected the third triangle
		b = true;
		rp = rs[0]/2.0;
		sp = 0.5 + rs[1]/2.0;
	}
	else if (IntersectTri(x4, r, n, rs, g, eps))
	{
		// we've intersected the fourth triangle
		b = true;
		rp = 0.5 - rs[0]/2.0;
		sp = 0.5 - rs[1]/2.0;
	}
	
	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[6], H1[6], H2[6];
		
		double l0;
		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 5;
		do
		{
			l0 = 1 - l1 - l2;
			
			// shape functions of quad
			H[0] = l0*(2.0*l0 - 1.0);
			H[1] = l1*(2.0*l1 - 1.0);
			H[2] = l2*(2.0*l2 - 1.0);
			H[3] = 4.0*l0*l1;
			H[4] = 4.0*l1*l2;
			H[5] = 4.0*l2*l0;
			
			// shape function derivatives
			H1[0] = -3.0 + 4.0*l1 + 4.0*l2;
			H1[1] =  4.0*l1 - 1.0;
			H1[2] =  0.0;
			H1[3] =  4.0 - 8.0*l1 - 4.0*l2;
			H1[4] =  4.0*l2;
			H1[5] = -4.0*l2;
			
			H2[0] = -3.0 + 4.0*l2 + 4.0*l1;
			H2[1] =  0.0;
			H2[2] =  4.0*l2 - 1.0;
			H2[3] = -4.0*l1;
			H2[4] =  4.0*l1;
			H2[5] =  4.0 - 8.0*l2 - 4.0*l1;

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3] - y[4]*H[4] - y[5]*H[5];
			
			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3] - y[4]*H1[4] - y[5]*H1[5];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3] - y[4]*H2[4] - y[5]*H2[5];
			F3 = n;
			
			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;
			
			// calculate solution increment
			dx = -(A.inverse()*F);
			
			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;
			
			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));
		
		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;
		
		// see if the point is inside the quad
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}
	
	return false;
}

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point. It searches for the closest patch based on
//! algebraic value of the gap function
//!
FESurfaceElement* FESurface::FindIntersection(vec3d r, vec3d n, double rs[2],
											  bool& binit_nq, double tol,
											  double srad)
{
	double g;
	int j;
	
	// see if we need to initialize the octree structure
	if (binit_nq) m_OT.Init();
	binit_nq = false;
	
	// let's find all the candidate surface elements
	set<int>selist;
	m_OT.FindCandidateSurfaceElements(r, n, selist);
	
	// now that we found candidate surface elements, lets see if we can find 
	// those that intersect the ray, then pick the closest intersection
	int iel;
	set<int>::iterator it;
	bool found = false;
	double rsl[2], gl;
	for (it=selist.begin(); it!=selist.end(); ++it) {
		// get the surface element
		j = *it;
		// project the node on the element
		if (Intersect(m_el[j], r, n, rsl, gl, tol)) {
			if ((!found) && (gl > -srad)) {
				found = true;
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				iel = j;
			} else if ((gl < g) && (gl > -srad)) {
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				iel = j;
			}
		}
	}
	if (found) return &m_el[iel];
	
	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point.  It searches for the closest patch based on
//! the absolute value of the gap function
//!
FESurfaceElement* FESurface::FindIntersection2(vec3d r, vec3d n, double rs[2],
											  bool& binit_nq, double tol,
											  double srad)
{
	double g;
	int j;
	
	// see if we need to initialize the octree structure
	if (binit_nq) m_OT.Init();
	binit_nq = false;
	
	// let's find all the candidate surface elements
	set<int>selist;
	m_OT.FindCandidateSurfaceElements(r, n, selist);
	
	// now that we found candidate surface elements, lets see if we can find 
	// those that intersect the ray, then pick the closest intersection
	int iel;
	set<int>::iterator it;
	bool found = false;
	double rsl[2], gl;
	for (it=selist.begin(); it!=selist.end(); ++it) {
		// get the surface element
		j = *it;
		// project the node on the element
		if (Intersect(m_el[j], r, n, rsl, gl, tol)) {
			if ((!found) && (fabs(gl) < srad)) {
				found = true;
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				iel = j;
			} else if ((fabs(gl) < fabs(g)) && (fabs(gl) < srad)) {
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				iel = j;
			}
		}
	}
	if (found) return &m_el[iel];
	
	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
//! Finds the element that contains the closest point projection of a node.
//! Returns zero if no such element can be found.
// TODO: I need to define a max search radius. For contact problems it is important that only nodes are
//       considered that are within an acceptable distance for contact.
FESurfaceElement* FESurface::ClosestPointProjection(vec3d& x, vec3d& q, vec2d& r, bool binit_nq, double tol)
{
	// get the mesh
	FEMesh& mesh = *GetMesh();

	// see if we need to initialize the NQ structure
	if (binit_nq) m_SNQ.Init();

	// let's find the closest master node
	int mn = m_SNQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = m_node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int nval = m_NET.Valence(mn);
	FEElement** pe = m_NET.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		int N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = ProjectToSurface(el, x, r[0], r[1]);
		if (IsInsideElement(el, r[0], r[1], tol)) return &el;
	}

	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search has failed:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	// TODO: I am not sure yet how to distinguish these cases. One way might be to project the point onto
	//       the edges and see if there is an edge that has a projection that is closer than the closest
	//       node. I am not sure yet if this is a sufficient argument though.
	return 0;
}

//-----------------------------------------------------------------------------
//! Finds the element that contains the closest point projection of a node.
//! Returns zero if no such element can be found.
// TODO: I need to define a max search radius. For contact problems it is important that only nodes are
//       considered that are within an acceptable distance for contact.
FESurfaceElement* FESurface::ClosestReferencePointProjection(vec3d& x, vec3d& q, vec2d& r, bool binit_nq, double tol)
{
	// get the mesh
	FEMesh& mesh = *GetMesh();

	// see if we need to initialize the NQ structure
	if (binit_nq) m_SNQ.InitReference();

	// let's find the closest master node
	int mn = m_SNQ.FindReference(x);

	// mn is a local index, so get the global node number too
	int m = m_node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_r0;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		int N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = ProjectToReferenceSurface(el, x, r[0], r[1]);
		if (IsInsideElement(el, r[0], r[1], tol)) return &el;
	}

	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search has failed:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	// TODO: I am not sure yet how to distinguish these cases. One way might be to project the point onto
	//       the edges and see if there is an edge that has a projection that is closer than the closest
	//       node. I am not sure yet if this is a sufficient argument though.
	return 0;
}

//-----------------------------------------------------------------------------
void FESurface::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		int ne = Elements();
		ar << ne;

		for (int k=0; k<ne; ++k)
		{
			FESurfaceElement& el = Element(k);
			ar << el.Type();
			ar << el.GetMatID() << el.m_nID << el.m_nrigid;
			ar << el.m_node;
			ar << el.m_lnode;
			ar << el.m_nelem;
		}
	}
	else
	{
		int ne=0;
		ar >> ne;
		create(ne);

		for (int k=0; k<ne; ++k)
		{
			FESurfaceElement& el = Element(k);

			int n, mat;
			ar >> n;
			el.SetType(n);

			ar >> mat >> el.m_nID >> el.m_nrigid;
			ar >> el.m_node;
			ar >> el.m_lnode;
			ar >> el.m_nelem;
			el.SetMatID(mat);
		}

		// initialize surface
		Init();
	}
}
