/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "FEContactPotential.h"
#include <FECore/FENode.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEBox.h>
#include <stdexcept>

vec3d MaterialPointPosition(FESurfaceElement& el, int n)
{
	FESurface* surf = dynamic_cast<FESurface*>(el.GetMeshPartition());
	vec3d r(0, 0, 0);
	double* H = el.H(n);
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d ri = surf->Node(el.m_lnode[i]).m_rt;
		r += ri * H[i];
	}
	return r;
}

void UpdateSurface(FESurface& surface)
{
// This assumes we are inside a omp parallel region!
#pragma omp for
	for (int i = 0; i < surface.Elements(); ++i)
	{
		FESurfaceElement& el = surface.Element(i);

		vec3d re[FEElement::MAX_NODES];
		surface.GetNodalCoordinates(el, re);

		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FECPContactPoint& mp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
			mp.m_rt = MaterialPointPosition(el, n);
			
			// kinematics at integration points
			mp.dxr = el.eval_deriv1(re, n);
			mp.dxs = el.eval_deriv2(re, n);

			mp.m_Jt = (mp.dxr ^ mp.dxs).norm();

			mp.m_gap = 0.0;
			mp.m_tc = vec3d(0, 0, 0);
		}
	}
}

FEContactPotentialSurface::FEContactPotentialSurface(FEModel* fem) : FEContactSurface(fem)
{

}

FEMaterialPoint* FEContactPotentialSurface::CreateMaterialPoint()
{
	return new FECPContactPoint;
}

void FEContactPotentialSurface::GetContactTraction(int nelem, vec3d& tc)
{
	FESurfaceElement& el = Element(nelem);
	tc = vec3d(0, 0, 0);
	for (int n = 0; n < el.GaussPoints(); ++n)
	{
		FECPContactPoint& mp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
		tc += mp.m_tc;
	}
	tc /= (double)el.GaussPoints();
}

double FEContactPotentialSurface::GetContactArea()
{
	double area = 0.0;
	for (int i = 0; i < Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		double* gw = el.GaussWeights();
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FECPContactPoint& mp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
			if (mp.m_tc.norm2() != 0.0)
			{
				vec3d dA = mp.dxr ^ mp.dxs;
				double da = dA.norm();
				area += da * gw[n];
			}
		}
	}
	return area;
}


BEGIN_FECORE_CLASS(FEContactPotential, FEContactInterface)
	ADD_PARAMETER(m_kc, "kc");
	ADD_PARAMETER(m_p, "p");
	ADD_PARAMETER(m_Rin, "R_in");
	ADD_PARAMETER(m_Rout, "R_out");
	ADD_PARAMETER(m_Rmin, "R0_min");
	ADD_PARAMETER(m_wtol, "w_tol");
END_FECORE_CLASS();

FEContactPotential::FEContactPotential(FEModel* fem) : FEContactInterface(fem), m_surf1(fem), m_surf2(fem)
{
	m_surf1.SetContactInterface(this);
	m_surf2.SetContactInterface(this);

	m_kc = 0.0;
	m_p = 4;
	m_Rin = 1.0;
	m_Rout = 2.0;
	m_Rmin = 0.0;
	m_wtol = 0.0;
}

//! return the primary surface
FESurface* FEContactPotential::GetPrimarySurface()
{
	return &m_surf1;
}

//! return the secondary surface
FESurface* FEContactPotential::GetSecondarySurface()
{
	return &m_surf2;
}

//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
bool FEContactPotential::UseNodalIntegration()
{
	return false;
}

static bool is_neighbor(FESurfaceElement& e1, FESurfaceElement& e2)
{
	int n1 = e1.Nodes();
	int n2 = e2.Nodes();
	for (int i = 0; i < n1; ++i)
	for (int j = 0; j < n2; ++j)
	{
		if (e1.m_node[i] == e2.m_node[j]) return true;
	}
	return false;
}

struct BOX
{
public:
	vec3d r0, r1;

public:
	BOX() {}
	BOX(const vec3d& p0, const vec3d& p1) : r0(p0), r1(p1) {}

	void add(const vec3d& r)
	{
		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;

		if (r.x > r1.x) r1.x = r.x;
		if (r.y > r1.y) r1.y = r.y;
		if (r.z > r1.z) r1.z = r.z;
	}

	void inflate(double l)
	{
		r0.x -= l; r0.y -= l; r0.z -= l;
		r1.x += l; r1.y += l; r1.z += l;
	}

	bool isInside(const vec3d& r) const
	{
		return (
			(r.x >= r0.x) && (r.x <= r1.x) &&
			(r.y >= r0.y) && (r.y <= r1.y) &&
			(r.z >= r0.z) && (r.z <= r1.z));
	}

	double MaxExtent() const
	{
		double sx = r1.x - r0.x;
		double sy = r1.y - r0.y;
		double sz = r1.z - r0.z;

		if ((sx >= sy) && (sx >= sz)) return sx;
		if ((sy >= sx) && (sy >= sz)) return sy;
		if ((sz >= sx) && (sz >= sy)) return sz;
		// we should never reach here
		return 0;
	}

	double width () const { return r1.x - r0.x; }
	double height() const { return r1.y - r0.y; }
	double depth () const { return r1.z - r0.z; }
};

class Grid
{
public:
	class Cell
	{
	public:
		Cell() {}
		Cell(const Cell& c) : m_box(c.m_box), m_elemList(c.m_elemList), m_nbr(c.m_nbr) {}
		void operator = (const Cell& c) { 
			m_box = c.m_box; 
			m_elemList = c.m_elemList; 
			m_nbr = c.m_nbr;
		}

		void add(FESurfaceElement* el)
		{
			m_elemList.insert(el);
		}

		bool empty() const { return m_elemList.empty(); }

	public:
		BOX m_box;
		set<FESurfaceElement*>	m_elemList;
		vector<Cell*>	m_nbr;
	};

	Cell* FindCell(const vec3d& r)
	{
		if (box.isInside(r) == false) return nullptr;

		int ix = (int)(m_nx * (r.x - box.r0.x) / box.width());
		int iy = (int)(m_ny * (r.y - box.r0.y) / box.height());
		int iz = (int)(m_nz * (r.z - box.r0.z) / box.depth());
		if (ix == m_nx) ix--;
		if (iy == m_ny) iy--;
		if (iz == m_nz) iz--;

		int icell = iz * (m_nx * m_ny) + iy * m_nx + ix;
		Cell* c = m_cell + icell;
		assert(c->m_box.isInside(r));
		return m_cell + icell;
	}

	int GetCellNeighborHood(const vec3d& r, Cell** cellList)
	{
		Cell* c = FindCell(r);
		if (c == nullptr) return 0;

		int n = 0;
		for (int i = 0; i < 27; ++i)
		{
			Cell* ci = c->m_nbr[i];
			if ((ci != nullptr) && (ci->empty() == false)) cellList[n++] = ci;
		}
		return n;
	}

public:
	Grid() { m_nx = m_ny = m_nz = 0; m_cell = nullptr; }

	bool Build(FESurface& s, int boxDivs, double minBoxSize)
	{
		// update the bounding box
		for (int i = 0; i < s.Nodes(); ++i)
		{
			FENode& node = s.Node(i);
			vec3d ri = node.m_rt;
			if (i == 0) box.r0 = box.r1 = ri;
			else box.add(ri);
		}

		// inflate a little, just to be sure
		double R = box.MaxExtent();
		box.inflate(minBoxSize);

		// determine the sizes
		double W = box.width();
		double H = box.height();
		double D = box.depth();

		double boxSize = box.MaxExtent() / boxDivs;
		if (boxSize < minBoxSize) boxSize = minBoxSize;

		m_nx = (int)(W / boxSize); if (m_nx < 1) m_nx = 1;
		m_ny = (int)(H / boxSize); if (m_ny < 1) m_ny = 1;
		m_nz = (int)(D / boxSize); if (m_nz < 1) m_nz = 1;

		// allocate the grid
		int ncells = m_nx * m_ny * m_nz;
		m_cell = new Cell[ncells];

		// helper lookup-table for cell neighbors
		const int LUT[27][3] = {
			{-1,-1,-1},{ 0,-1,-1},{1,-1,-1},
			{-1, 0,-1},{ 0, 0,-1},{1, 0,-1},
			{-1, 1,-1},{ 0, 1,-1},{1, 1,-1},
			{-1,-1, 0},{ 0,-1, 0},{1,-1, 0},
			{-1, 0, 0},{ 0, 0, 0},{1, 0, 0},
			{-1, 1, 0},{ 0, 1, 0},{1, 1, 0},
			{-1,-1, 1},{ 0,-1, 1},{1,-1, 1},
			{-1, 0, 1},{ 0, 0, 1},{1, 0, 1},
			{-1, 1, 1},{ 0, 1, 1},{1, 1, 1} };

		// construct cells
		Cell* c = m_cell;
		double X0 = box.r0.x;
		double Y0 = box.r0.y;
		double Z0 = box.r0.z;
		for (int k = 0; k < m_nz; ++k)
		{
			double z0 = Z0 + k * D / m_nz;
			double z1 = Z0 + (k+1) * D / m_nz;
			for (int j = 0; j < m_ny; ++j)
			{
				double y0 = Y0 + j * H / m_ny;
				double y1 = Y0 + (j + 1) * H / m_ny;
				for (int i = 0; i < m_nx; ++i, ++c)
				{
					double x0 = X0 + i * W / m_nx;
					double x1 = X0 + (i + 1) * W / m_nx;

					c->m_box = BOX(vec3d(x0, y0, z0), vec3d(x1, y1, z1));
					double R = c->m_box.MaxExtent();
					c->m_box.inflate(R * 1e-4);

					// assign neighbors
					c->m_nbr.assign(27, nullptr);
					for (int n = 0; n < 27; ++n)
					{
						const int* l = LUT[n];
						c->m_nbr[n] = GetCell(i + l[0], j + l[1], k + l[2]);
					}
				}
			}
		}

		// assign elements to grid cells
		for (int i = 0; i < s.Elements(); ++i)
		{
			FESurfaceElement& el = s.Element(i);
			int nint = el.GaussPoints();
			for (int n = 0; n < nint; ++n)
			{
				FECPContactPoint& mp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
				Cell* c = FindCell(mp.m_rt); assert(c);
				if (c == nullptr) return false;
				c->add(&el);
			}
		}

		return true;
	}

	Cell* GetCell(int i, int j, int k) const
	{
		if ((i < 0) || (i >= m_nx)) return nullptr;
		if ((j < 0) || (j >= m_ny)) return nullptr;
		if ((k < 0) || (k >= m_nz)) return nullptr;
		return m_cell + (k * (m_nx * m_ny) + j * m_nx + i);
	}

	~Grid() { delete [] m_cell; }

protected:
	BOX		box;
	int		m_nx, m_ny, m_nz;
	Cell*	m_cell;
};

// initialization
bool FEContactPotential::Init()
{
	if (FEContactInterface::Init() == false) return false;
	BuildNeighborTable();
	return true;
}

void FEContactPotential::BuildNeighborTable()
{
	m_elemNeighbors.resize(m_surf1.Elements());
	for (int i = 0; i < m_surf1.Elements(); ++i)
	{
		FESurfaceElement& el1 = m_surf1.Element(i);

		set<FESurfaceElement*>& nbrList = m_elemNeighbors[i];
		nbrList.clear();

		for (int j = 0; j < m_surf2.Elements(); ++j)
		{
			FESurfaceElement& el2 = m_surf2.Element(j);
			if (is_neighbor(el1, el2))
			{
				nbrList.insert(&el2);
			}
		}
	}
}

// update
void FEContactPotential::Update()
{
	FEContactInterface::Update();

	// update the constants
	m_c1 = 0.5 * m_p * m_kc / ((m_Rout - m_Rin) * pow(m_Rin, m_p + 1.0));
	m_c2 = m_kc / pow(m_Rin, m_p) - m_c1 * (m_Rin - m_Rout) * (m_Rin - m_Rout);

	// Update the surfaces
#pragma omp parallel 
	{
		UpdateSurface(m_surf1);
		UpdateSurface(m_surf2);
	}

	// build the grid
	int ndivs = (int)pow(m_surf2.Elements(), 0.33333);
	if (ndivs < 2) ndivs = 2;
	Grid g;
	if (g.Build(m_surf2, ndivs, m_Rout) == false)
	{
		throw std::runtime_error("Failed to build grid in FEContactPotential::Update");
	}

	// build the list of active elements
	m_activeElements.resize(m_surf1.Elements());
#pragma omp parallel for shared(g) schedule(dynamic)
	for (int i = 0; i < m_surf1.Elements(); ++i)
	{
		FESurfaceElement& el1 = m_surf1.Element(i);

		set<FESurfaceElement*>& activeElems = m_activeElements[i];
		activeElems.clear();

		// list of elements to exclude. This will be neighbors and elements
		// already processed
		set<FESurfaceElement*> excludeList = m_elemNeighbors[i];

		for (int n = 0; n < el1.GaussPoints(); ++n)
		{
			FECPContactPoint& mp1 = static_cast<FECPContactPoint&>(*el1.GetMaterialPoint(n));
			mp1.m_gap = 0.0;
			vec3d r1 = mp1.m_rt;
			vec3d R1 = mp1.m_r0;
			vec3d n1 = mp1.dxr ^ mp1.dxs; n1.unit();

			// find the grid cell this point is in and loop over the cell's neighborhood
			Grid::Cell* c[27] = { nullptr };
			int nc = g.GetCellNeighborHood(r1, &c[0]);
			for (int l = 0; l < nc; ++l)
			{
				Grid::Cell* cl = c[l];
				for (FESurfaceElement* el2 : cl->m_elemList)
				{
					// make sure we did not process this element yet
					// and the element is not a neighbor (which can be the case for self-contact)
					if (excludeList.find(el2) == excludeList.end())
					{
						// Next, we see if any integration point of el2 is close to the current 
						// integration point of el1. 
						vec3d r12;
						for (int m = 0; m < el2->GaussPoints(); ++m)
						{
							FEMaterialPoint* mp2 = el2->GetMaterialPoint(m);
							vec3d r2 = mp2->m_rt;
							vec3d R2 = mp2->m_r0;

							r12.x = r1.x - r2.x;
							r12.y = r1.y - r2.y;
							r12.z = r1.z - r2.z;
//							vec3d r12 = r1 - mp2.m_rt;
							if ((r12.x < m_Rout) && (r12.x > -m_Rout) &&
								(r12.y < m_Rout) && (r12.y > -m_Rout) &&
								(r12.z < m_Rout) && (r12.z > -m_Rout) &&
								(r12.norm2() < m_Rout * m_Rout))
							{
								double L12 = (R2 - R1).norm2();
								double l12 = r12.unit();
								if ((fabs(r12 * n1) >= m_wtol) && (L12 >= m_Rmin))
								{
									// we found one, so insert it to the list of active elements
									activeElems.insert(el2);

									// also insert it to the exclude list
									excludeList.insert(el2);

//									fprintf(stderr, "r12: %lg, %lg, %lg\n", r12.x, r12.y, r12.z);
//									fprintf(stderr, "n1 : %lg, %lg, %lg\n", n1.x, n1.y, n1.z);
//									fprintf(stderr, "l12 = %lg\n", l12);

									if ((mp1.m_gap == 0.0) || (l12 < mp1.m_gap))
									{
										mp1.m_gap = l12;
									}
									break;
								}
							}
						}

					}
				}
			}
		}
	}
}

// Build the matrix profile
void FEContactPotential::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// connect every element of surface 1 to surface 2
	for (int i = 0; i < m_surf1.Elements(); ++i)
	{
		FESurfaceElement& el1 = m_surf1.Element(i);

		// add the dofs of element 1
		vector<int> lm;
		for (int j = 0; j < el1.Nodes(); ++j)
		{
			FENode& node = m_surf1.Node(el1.m_lnode[j]);
			lm.push_back(node.m_ID[0]);
			lm.push_back(node.m_ID[1]);
			lm.push_back(node.m_ID[2]);
		}

		// add all active dofs of surface 2
		set<FESurfaceElement*>& activeElems = m_activeElements[i];
		for (FESurfaceElement* el2 : activeElems)
		{
			for (int j = 0; j < el2->Nodes(); ++j)
			{
				FENode& node = m_surf2.Node(el2->m_lnode[j]);
				lm.push_back(node.m_ID[0]);
				lm.push_back(node.m_ID[1]);
				lm.push_back(node.m_ID[2]);
			}
		}
		M.build_add(lm);
	}
}

double FEContactPotential::PotentialDerive(double r)
{
	double f = 0.0;
	if (r < m_Rin)
	{
		f = (m_kc / pow(r, m_p)) - m_c2;
	}
	else if (r < m_Rout)
	{
		f = m_c1 * (r - m_Rout) * (r - m_Rout);
	}
	return f;
}

double FEContactPotential::PotentialDerive2(double r)
{
	double f = 0.0;
	if (r < m_Rin)
	{
		f = -m_p*m_kc / pow(r, m_p + 1.0);
	}
	else if (r < m_Rout)
	{
		f = 2.0*m_c1 * (r - m_Rout);
	}
	return f;
}

// The LoadVector function evaluates the "forces" that contribute to the residual of the system
void FEContactPotential::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	const int ndof = 3;

	// clear all contact tractions
#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < m_surf1.Elements(); ++i)
		{
			FESurfaceElement& el = m_surf1.Element(i);
			for (int n = 0; n < el.GaussPoints(); ++n)
			{
				FECPContactPoint& cp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
				cp.m_Ln = 0.0;
				cp.m_tc = vec3d(0, 0, 0);
			}
		}

		#pragma omp for
		for (int i = 0; i < m_surf2.Elements(); ++i)
		{
			FESurfaceElement& el = m_surf2.Element(i);
			for (int n = 0; n < el.GaussPoints(); ++n)
			{
				FECPContactPoint& cp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
				cp.m_Ln = 0.0;
				cp.m_tc = vec3d(0, 0, 0);
			}
		}
	}

	// loop over all elements of surf 1
#pragma omp parallel for shared(R)
	for (int i = 0; i < m_surf1.Elements(); ++i)
	{
		FESurfaceElement& eli = m_surf1.Element(i);
		int na = eli.Nodes();

		vector<double> fe;
		vector<int> lm;

		// loop over all elements of surf 2
		set<FESurfaceElement*>& activeElems = m_activeElements[i];
		for (FESurfaceElement* elj : activeElems)
		{
			int nb = elj->Nodes();

			// evaluate contribution to force vector
			fe.assign((na + nb) * ndof, 0.0);
			ElementForce(eli, *elj, fe);

			// setup the lm vector
			lm.resize(3 * (na + nb));
			for (int a = 0; a < na; ++a)
			{
				FENode& node = m_surf1.Node(eli.m_lnode[a]);
				lm[3*a    ] = node.m_ID[0];
				lm[3*a + 1] = node.m_ID[1];
				lm[3*a + 2] = node.m_ID[2];
			}
			for (int b = 0; b < nb; ++b)
			{
				FENode& node = m_surf2.Node(elj->m_lnode[b]);
				lm[3*na + 3*b    ] = node.m_ID[0];
				lm[3*na + 3*b + 1] = node.m_ID[1];
				lm[3*na + 3*b + 2] = node.m_ID[2];
			}

			// assemble
			R.Assemble(lm, fe);
		}
	}

	// update contact pressures (only needed for plot output)
#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < m_surf1.Elements(); ++i)
		{
			FESurfaceElement& el = m_surf1.Element(i);
			for (int n = 0; n < el.GaussPoints(); ++n)
			{
				FECPContactPoint& cp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
				cp.m_Ln = cp.m_tc.norm();
			}
		}

		#pragma omp for
		for (int i = 0; i < m_surf2.Elements(); ++i)
		{
			FESurfaceElement& el = m_surf2.Element(i);
			for (int n = 0; n < el.GaussPoints(); ++n)
			{
				FECPContactPoint& cp = static_cast<FECPContactPoint&>(*el.GetMaterialPoint(n));
				cp.m_Ln = cp.m_tc.norm();
			}
		}
	}
}

void FEContactPotential::ElementForce(FESurfaceElement& el1, FESurfaceElement& el2, vector<double>& fe)
{
	int na = el1.Nodes();
	int nb = el2.Nodes();
	double* w1 = el1.GaussWeights();
	double* w2 = el2.GaussWeights();

	for (int n = 0; n < el1.GaussPoints(); ++n)
	{
		const double* H1 = el1.H(n);
		FECPContactPoint& mp1 = static_cast<FECPContactPoint&>(*el1.GetMaterialPoint(n));
		double Jw1 = mp1.m_Jt*w1[n];
		vec3d r1 = mp1.m_rt;

		for (int m = 0; m < el2.GaussPoints(); ++m)
		{
			const double* H2 = el2.H(m);
			FECPContactPoint& mp2 = static_cast<FECPContactPoint&>(*el2.GetMaterialPoint(m));
			double Jw2 = mp2.m_Jt*w2[m];

			double Jw12 = Jw1 * Jw2;

			// get the position of the integration point
			vec3d r2 = mp2.m_rt;

			// get the normalized direction vector
			vec3d e12 = r1 - r2;
			double r12 = e12.unit();

			// calculate the potential derivative
			double df = PotentialDerive(r12);

			// add to the force
			if (df != 0.0)
			{
				vec3d F = e12 * df;
				mp1.m_tc += F * Jw2;
				mp2.m_tc -= F * Jw1;

				for (int a = 0; a < na; ++a)
				{
					fe[3 * a    ] += F.x * (H1[a] * Jw12);
					fe[3 * a + 1] += F.y * (H1[a] * Jw12);
					fe[3 * a + 2] += F.z * (H1[a] * Jw12);
				}

				for (int b = 0; b < nb; ++b)
				{
					fe[3 * (na + b)    ] -= F.x * (H2[b] * Jw12);
					fe[3 * (na + b) + 1] -= F.y * (H2[b] * Jw12);
					fe[3 * (na + b) + 2] -= F.z * (H2[b] * Jw12);
				}
			}
		}
	}
}

// Evaluates the contriubtion to the stiffness matrix
void FEContactPotential::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	const int ndof = 3;

	// loop over all elements of surf 1
#pragma omp parallel for shared(LS)
	for (int i = 0; i < m_surf1.Elements(); ++i)
	{
		FESurfaceElement& eli = m_surf1.Element(i);
		int na = eli.Nodes();

		set<FESurfaceElement*>& activeElems = m_activeElements[i];
		for (FESurfaceElement* elj : activeElems)
		{
			int nb = elj->Nodes();

			FEElementMatrix ke((na + nb) * ndof, (na + nb) * ndof);
			ke.zero();

			ElementStiffness(eli, *elj, ke);
			
			vector<int> lm(3 * (na+nb));
			vector<int> en(na + nb);
			for (int a = 0; a < na; ++a)
			{
				FENode& node = m_surf1.Node(eli.m_lnode[a]);
				lm[3 * a    ] = node.m_ID[0];
				lm[3 * a + 1] = node.m_ID[1];
				lm[3 * a + 2] = node.m_ID[2];
				en[a] = eli.m_node[a];
			}
			for (int b = 0; b < nb; ++b)
			{
				FENode& node = m_surf2.Node(elj->m_lnode[b]);
				lm[3 * na + 3 * b    ] = node.m_ID[0];
				lm[3 * na + 3 * b + 1] = node.m_ID[1];
				lm[3 * na + 3 * b + 2] = node.m_ID[2];
				en[na + b] = elj->m_node[b];
			}

			ke.SetIndices(lm);
			ke.SetNodes(en);

			LS.Assemble(ke);
		}
	}
}

void FEContactPotential::ElementStiffness(FESurfaceElement& el1, FESurfaceElement& el2, matrix& ke)
{
	double* w1 = el1.GaussWeights();
	double* w2 = el2.GaussWeights();

	int na = el1.Nodes();
	int nb = el2.Nodes();

	for (int n = 0; n < el1.GaussPoints(); ++n)
	{
		const double* H1 = el1.H(n);
		FECPContactPoint& mp1 = static_cast<FECPContactPoint&>(*el1.GetMaterialPoint(n));
		double Jw1 = mp1.m_Jt*w1[n];

		vec3d r1 = mp1.m_rt;

		for (int m = 0; m < el2.GaussPoints(); ++m)
		{
			const double* H2 = el2.H(m);
			FECPContactPoint& mp2 = static_cast<FECPContactPoint&>(*el2.GetMaterialPoint(m));
			double Jw2 = mp2.m_Jt*w2[m];

			vec3d r2 = mp2.m_rt;

			// get the normalized direction vector
			vec3d e12 = r1 - r2;
			double r12 = e12.unit();

			// calculate the potential's derivatives
			double df  = PotentialDerive(r12);
			double d2f = PotentialDerive2(r12);

			// add to the stiffness
			mat3dd I(1.0);
			mat3d ExE = dyad(e12);
			mat3d P = (I - ExE) / r12;
			mat3d M = ExE * d2f + P * df;

			// 1-1 contribution
			for (int a = 0; a < na; ++a)
			{
				for (int b = 0; b < na; ++b)
				{
					mat3d K2 = M * Jw2;
					ke.add(3 * a, 3 * b, K2 * (H1[a] * H1[b])*Jw1);
				}
			}

			// 1-2 contribution
			for (int a = 0; a < na; ++a)
			{
				for (int b = 0; b < nb; ++b)
				{
					mat3d K2b = M * (H2[b] * Jw2);
					ke.sub(3 * a, 3 * na + 3 * b, K2b * (H1[a] * Jw1));
				}
			}

			// 2-1 contribution
			for (int b = 0; b < nb; ++b)
			{
				for (int a = 0; a < na; ++a)
				{
					mat3d K1a = M * (H1[a] * Jw1);
					ke.sub(3 * na + 3*b, 3 * a, K1a * (H2[b] * Jw2));
				}
			}

			// 2-2 contribution
			for (int a = 0; a < nb; ++a)
			{
				for (int b = 0; b < nb; ++b)
				{
					mat3d K1 = M * Jw1;
					ke.add(3*na + 3 * a, 3*na + 3 * b, K1 * (H2[a] * H2[b]) * Jw2);
				}
			}
		}
	}

	ke *= -1.0;
}

void FEContactPotential::Serialize(DumpStream& ar)
{
	FEContactInterface::Serialize(ar);

	m_surf1.Serialize(ar);
	m_surf2.Serialize(ar);

	BuildNeighborTable();
}
