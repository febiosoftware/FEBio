#include "stdafx.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! Finds the (master) element that contains the projection of a (slave) node

FEElement* FEContactSurface::FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool binit_nq)
{
	// get the mesh
	FEMesh& mesh = *m_pmesh;

	// see if we need to initialize the NQ structure
	if (binit_nq) m_NQ.Init();

	// let's find the closest master node
	int mn = m_NQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	const double eps = 0.01;
	int N;

	// loop over all master elements that contain the node mn
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = ProjectToSurface(el, x, r[0], r[1]);

		// see if we have a winner
		if (N == 4)
		{
			if ((r[0] >= -1-eps) && (r[0] <= 1+eps) && (r[1] >= -1-eps) && (r[1] <= 1+eps))
			{
				// we have a winner!
				return pe[j];
			}
		}
		else if (N == 3)
		{
			if ((r[0] >= -eps) && (r[1] >= -eps) && (r[0]+r[1] <= 1 + eps))
			{
				// we have a winner
				return pe[j];
			}
		} 
		else assert(false);
	}

	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

void FEContactSurface::Init()
{
	int i, j, n;

	// always intialize base class first!
	FESurface::Init();

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	gap.create(nn);		// gap funtion
	nu.create(nn);		// node normal 
	pme.create(nn);		// penetrated master element
	rs.create(nn);		// natural coords of projected slave node on master element
	rsp.create(nn);
	Lm.create(nn);
	M.create(nn);
	Lt.create(nn);
	off.create(nn);
	eps.create(nn);

	// set initial values
	gap.zero();
	pme.set(0);
	Lm.zero();
	off.zero();
	eps.set(1.0);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pmesh;
	vector<double> tag(m.Nodes());
	tag.zero();
	for (i=0; i<m.ShellElements(); ++i)
	{
		FEShellElement& el = m.ShellElement(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
	}
	for (i=0; i<nn; ++i) off[i] = tag[node[i]];
}

//-----------------------------------------------------------------------------
//! Update surface data

void FEContactSurface::Update()
{
}
