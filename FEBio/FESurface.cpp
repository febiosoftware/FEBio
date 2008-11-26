// FESurface.cpp: implementation of the FESurface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESurface.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! Initialize surface node data structure
//! Note that it is assumed that the element array is already created
//! and initialized.

void FESurface::Init()
{
	int i, j, m;

	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pmesh;

	// This array is used to keep tags on each node
	vector<int> tag(mesh.Nodes());
	tag.set(-1);

	// let's find all nodes the surface needs
	int nn = 0;
	int ne = Elements();
	for (i=0; i<ne; ++i)
	{
		FESurfaceElement& el = Element(i);

		for (j=0; j<el.Nodes(); ++j)
		{
			// get the global node number
			m = el.m_node[j];
		
			// create a local node number
			if (tag[m] == -1) tag[m] = nn++;

			// set the local node number
			el.m_lnode[j] = tag[m];
		}
	}

	// allocate node index table
	node.create(nn);

	// fill the node index table
	for (i=0; i<mesh.Nodes(); ++i)
	{
		if (tag[i] >= 0)
		{
			node[tag[i]] = i;
		}
	}

	// create the node element list
	m_NEL.Create(*this);
}

//-----------------------------------------------------------------------------
//! This function calculates the projection of x on the surface element el.
//! It does this by finding the solution of the nonlinear equation (x-y)*y,[a]=0,
//! where the comma denotes differentation and a ranges from 1 to 2.
//! The system is solved using the Newton-Raphson method.
//! The surface element may be either a quad or a triangular element.

vec3d FESurface::ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s)
{
	double R[2], u[2], D;

	vec3d q(0,0,0), y[4];

	double gr[4] = {-1, +1, +1, -1};
	double gs[4] = {-1, -1, +1, +1};
	double H[4], Hr[4], Hs[4], Hrs[4];
	double normu;

	int i, j;
	int NMAX = 5, n=0;

	// number of element nodes
	int ne = el.Nodes();

	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pmesh;

	// get the elements nodal positions
	for (i=0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;

	// loop until converged
	do
	{
		double A[2][2] = {0};

		if (ne == 4)
		{
			// do quadrilaterals
			for (i=0; i<4; ++i)
			{
				H[i] = 0.25*(1+gr[i]*r)*(1+gs[i]*s);
	
				Hr[i] = 0.25*gr[i]*( 1 + gs[i]*s );
				Hs[i] = 0.25*gs[i]*( 1 + gr[i]*r );

				Hrs[i] = 0.25*gr[i]*gs[i];
			}
		}
		else
		{
			// do triangles
			H[0] = 1 - r - s;
			H[1] = r;
			H[2] = s;
			Hr[0] = -1; Hs[0] = -1;
			Hr[1] =  1; Hs[1] =  0;
			Hr[2] =  0; Hs[2] =  1;
			Hrs[0] = Hrs[1] = Hrs[2] = 0;
		}

		// set up the system of equations
		q = vec3d(0,0,0);
		R[0] = R[1] = 0;
		for (i=0; i<ne; ++i)
		{
			R[0] -= (x*y[i])*Hr[i];
			R[1] -= (x*y[i])*Hs[i];

			A[0][1] += (x*y[i])*Hrs[i];
			A[1][0] += (x*y[i])*Hrs[i];

			for (j=0; j<ne; ++j)
			{
				R[0] -= -H[i]*Hr[j]*(y[i]*y[j]);
				R[1] -= -H[i]*Hs[j]*(y[i]*y[j]);

				A[0][0] += -(y[i]*y[j])*(Hr[i]*Hr[j]);
				A[1][1] += -(y[i]*y[j])*(Hs[i]*Hs[j]);

				A[0][1] += -(y[i]*y[j])*(Hs[j]*Hr[i]+H[i]*Hrs[j]);
				A[1][0] += -(y[i]*y[j])*(Hr[j]*Hs[i]+H[i]*Hrs[j]);
			}
		
			q += y[i]*H[i];
		}
	
		// determinant of A
		D = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// solve for u = A^(-1)*R
		u[0] = (A[1][1]*R[0] - A[0][1]*R[1])/D;
		u[1] = (A[0][0]*R[1] - A[1][0]*R[0])/D;

		normu = sqrt(u[0]*u[0]+u[1]*u[1]);
	
		r += u[0];
		s += u[1];

		++n;
	}
	while ((normu > 1e-5) && (n < NMAX));

	return q;
}

//-----------------------------------------------------------------------------
//! This function calculates the area of a surface element

double FESurface::FaceArea(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pmesh;

	// unpack surface element data
	mesh.UnpackElement(el);

	// get the number of nodes
	int neln = el.Nodes();

	// get the initial nodes
	vec3d* r0 = el.r0();

	// get the integration weights
	double* w = el.GaussWeights();

	double *Gr, *Gs;
	vec3d dxr, dxs;

	double detJ;

	double area = 0;

	int n, k;

	for (n=0; n<neln; ++n)
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
//! Calculates the metric tensor at the point with surface coordinates (r,s)
//! Note that we assume that the surface element is unpacked.
// TODO: perhaps I should place this function in the element class

mat2d FESurface::Metric0(FESurfaceElement& el, double r, double s)
{
	// make sure the element is unpacked
	assert(el.IsUnpacked());

	// nr of element nodes
	int neln = el.Nodes();

	// shape function derivatives
	double Hr[4], Hs[4];

	// get the shape function values at this slave node
	if (neln == 4)
	{
		Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
		Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
		Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
		Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
	}
	else if (neln == 3)
	{
		Hr[0] = -1; Hs[0] = -1;
		Hr[1] =  1; Hs[1] =  0;
		Hr[2] =  0; Hs[2] =  1;
	}
	else assert(false);

	// reference nodal coordinates
	vec3d* r0 = el.r0();

	// get the tangent vectors
	vec3d t1(0,0,0);
	vec3d t2(0,0,0);
	for (int k=0; k<neln; ++k)
	{
		t1.x += Hr[k]*r0[k].x;
		t1.y += Hr[k]*r0[k].y;
		t1.z += Hr[k]*r0[k].z;
		
		t2.x += Hs[k]*r0[k].x;
		t2.y += Hs[k]*r0[k].y;
		t2.z += Hs[k]*r0[k].z;
	}

	// calculate metric tensor
	return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}
