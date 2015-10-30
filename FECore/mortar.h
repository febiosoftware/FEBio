#pragma once
#include "FECore/vec3d.h"
#include "FESurface.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// Structure for representing 2D points
struct POINT2D
{
	double x, y;
};

//-----------------------------------------------------------------------------
// Calculates the intersection between two convex polygons
int ConvexIntersect(POINT2D* P, int n, POINT2D* Q, int m, POINT2D* R);

//-----------------------------------------------------------------------------
class Patch
{
public:
	struct FACET
	{
		FACET(){}
		FACET(vec3d x[3]) { r[0] = x[0]; r[1] = x[1]; r[2] = x[2]; }
		vec3d	r[3];	//!< position of nodes

		// Evaluate the spatial position using iso-parametric coordinates
		vec3d Position(double rp, double sp)
		{
			return (r[0]*(1.0 - rp - sp) + r[1]*rp + r[2]*sp);
		}

		//! calculate the area of the patch
		double Area()
		{
			vec3d n = (r[1] - r[0])^(r[2] - r[0]);
			return n.norm()*0.5;
		}
	};

public:
	//! Clear the patch
	void Clear() { m_tri.clear(); }

	//! Add a facet to the patch
	void Add(vec3d x[3]) { m_tri.push_back(FACET(x)); }

	//! retrieve a facet
	FACET& Facet(int i) { return m_tri[i]; }

	//! facet count
	int Size() { return (int) m_tri.size(); }

	//! See if the facet is empty
	bool Empty() { return m_tri.empty(); }

private:
	vector<FACET>	m_tri;	//!< triangular patches
};

bool CalculateMortarIntersection(FESurface& ss, FESurface& ms, int k, int l, Patch& patch);
