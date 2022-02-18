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
#include <assert.h>
#include "mortar.h"
#include <math.h>
#include "FEMesh.h"

//-----------------------------------------------------------------------------
// subtract operator for POINT2D
POINT2D operator - (POINT2D a, POINT2D b)
{
	POINT2D p = {a.x - b.x, a.y - b.y};
	return p;
}

//-----------------------------------------------------------------------------
// dot product for POINT2D
double operator * (POINT2D a, POINT2D b)
{
	return a.x*b.x + a.y*b.y;
}

//-----------------------------------------------------------------------------
// This function returns twice the area of the triangle defined by [a,b,c].
double Area2(POINT2D a, POINT2D b, POINT2D c)
{
	return ((b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y));
}

//-----------------------------------------------------------------------------
// This function returns the sign of the area, defined by three POINT2Ds. 
int AreaSign(POINT2D a, POINT2D b, POINT2D c, const double eps)
{
	double d = ((b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y));
	if (d > eps) return 1;
	else if (d < -eps) return -1;
	return 0;
}

//-----------------------------------------------------------------------------
// Checks collinearity of three POINT2Ds
bool Collinear(POINT2D a, POINT2D b, POINT2D c, const double eps)
{
	return fabs(Area2(a, b, c)) <= eps;
}

//-----------------------------------------------------------------------------
// This function checks if c is between a, b. 
// NOTE: This function assumes that the three points are collinear
bool Between(POINT2D a, POINT2D b, POINT2D c)
{
	if (a.x != b.x)
	{
		return (((a.x <= c.x)&&(c.x <= b.x))||
				((a.x >= c.x)&&(c.x >= b.x)));
	}
	else
	{
		return (((a.y <= c.y)&&(c.y <= b.y))||
				((a.y >= c.y)&&(c.y >= b.y)));
	}
}

//-----------------------------------------------------------------------------
// This function checks if the two segments (a,b) and (c,d) are collinear
// and if so, returns a POINT2D of intersection. The return code can be used 
// to identify the intersection type
// 0 = no overlap
// 3 = overlap
// NOTE: This function assumes that the segments are parallel
int ParallelInt(POINT2D a, POINT2D b, POINT2D c, POINT2D d, POINT2D& p, const double eps)
{
	// make sure POINT2Ds are collinear
	if (!Collinear(a,b,c, eps)) return 0;

	if (Between(a, b, c)) { p = c; return 3; }
	if (Between(a, b, d)) { p = d; return 3; }
	if (Between(c, d, a)) { p = a; return 3; }
	if (Between(c, d, b)) { p = b; return 3; }

	return 0;
}

//-----------------------------------------------------------------------------
// This function calculates the segment-segment intersection. The first segment
// is defined by POINT2Ds (a,b) and the second segment by (c,d). The intersection
// POINT2D (if it exists) is returned in p. The return code can be used to identify
// the intersection type:
// 0 = no intersection
// 1 = proper intersection
// 2 = vertex intersection (segments coincide at vertex)
// 3 = edge intersections (segments cooincide)
int SegSegInt(POINT2D a, POINT2D b, POINT2D c, POINT2D d, POINT2D& p, const double eps)
{
	// return code
	int nret = -1;

	// calculate the denominator
	double denom = a.x*(d.y - c.y) + b.x*(c.y - d.y) + d.x*(b.y - a.y) + c.x*(a.y - b.y);

	// if the denominator is zero, the segments are parallel; handle separately
	if (fabs(denom) <= eps) return ParallelInt(a, b, c, d, p, eps);

	double num = a.x*(d.y - c.y) + c.x*(a.y - d.y) + d.x*(c.y - a.y);
	if ((fabs(num)<=eps)||(num==denom)) nret = 2;
	double s = num / denom;

	num = -(a.x*(c.y - b.y) + b.x*(a.y - c.y) + c.x*(b.y - a.y));
	if ((fabs(num)<=eps)||(num==denom)) nret = 2;
	double t = num / denom;

	// check for proper intersection
	if ((0.0 < s)&&(s < 1.0)&&(0.0 < t)&&(t < 1.0)) nret = 1;
	else if ((0.0 > s)||(s > 1.0) || (0.0 > t) || (t > 1.0)) nret = 0;

	// find the intersection POINT2D
	p.x = a.x + s*(b.x - a.x);
	p.y = a.y + s*(b.y - a.y);

	assert(nret >= 0);
	return nret;
}

//-----------------------------------------------------------------------------
// See if a point is inside a (counter-clockwise) polygon
bool PointInConvexPoly(POINT2D p, POINT2D* P, int n)
{
	for (int i=0; i<n; ++i)
	{
		int i1 = (i+1)%n;
		if (Area2(P[i],P[i1],p) < 0) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// see if a convex polygon P is inside the convex polygon Q
bool ConvexPolyInConvexPoly(POINT2D* P, int n, POINT2D* Q, int m)
{
	for (int i=0; i<n; ++i)
	{
		if (PointInConvexPoly(P[i], Q, m) == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
int InOut(int inFlag, int aHB, int bHA)
{
	if      (aHB > 0) return -1;
	else if (bHA > 0) return  1;
	else return inFlag;
}

int Advance(int a, int *aa, int n)
{
	(*aa)++;
	return (a+1)%n;
}

//-----------------------------------------------------------------------------
// Calculates the intersection between two convex polygons
int ConvexIntersect(POINT2D* P, int n, POINT2D* Q, int m, POINT2D* R)
{
	int		a, b;			// indices on P and Q (resp)
	int		a1, b1;			// a-1, b-1 (resp)
	POINT2D	A, B;			// directed edges on P, Q (resp)
	int		cross;			// sign of z-component of A x B
	int		bHA, aHB;		// b in H(a); a in H(b)
	POINT2D	O = {0,0};		// origin
	POINT2D	p;				// intersection points
	int		inFlag;			// -1,0,1 = InP, Unknown, InQ
	int		aa, ba;			// number of advances on a & b indices (after 1st iter)
	bool	FirstPOINT2D;		// is this the first POINT2D (used to initialize)
	POINT2D	p0;				// the first POINT2D
	int		code;			// SegSegInt return code
	int		nr;				// number of POINT2Ds in intersection

	const double eps = 1e-7;

	// Initialize the variables
	nr = 0;
	a = 0; b = 0; aa = 0; ba = 0;
	inFlag = 0; FirstPOINT2D = true;

	do 
	{
		// computation of key variables
		a1 = (a + n - 1) % n;
		b1 = (b + m - 1) % m;

		A = P[a] - P[a1];
		B = Q[b] - Q[b1];

		cross = AreaSign(O, A, B, eps);
		aHB = AreaSign(Q[b1], Q[b], P[a], eps);
		bHA = AreaSign(P[a1], P[a], Q[b], eps);

		// if A,B intersect, update inflag
		code = SegSegInt(P[a1], P[a], Q[b1], Q[b], p, eps);
		if ((code==1)||(code==2))
		{
			if ((inFlag == 0) && FirstPOINT2D)
			{
				aa = ba = 0;
				FirstPOINT2D = false;
				p0 = p;
			}
			R[nr++] = p;
			inFlag = InOut(inFlag, aHB, bHA);
		}

		// --- Advance rules ---
		if ((code == 3) && (A*B < 0))
		{
			// P,Q intersect in this edge only
			break;
		}

		if ((cross == 0) && (aHB < 0) && (bHA < 0))
		{
//			printf("P,Q are disjoint\n"); 
			break;
		}
		else if ((cross == 0) && (aHB == 0) && (bHA == 0))
		{
			if (inFlag == -1)
				b = Advance(b, &ba, m);
			else
				a = Advance(a, &aa, n);
		}
		else if (cross >= 0)
		{
			if (bHA > 0) 
			{
				if (inFlag == -1) R[nr++] = P[a];
				a = Advance(a, &aa, n);
			}
			else 
			{
				if (inFlag == 1) R[nr++] = Q[b];
				b = Advance(b, &ba, m);
			}
		}
		else // if cross < 0
		{
			if (aHB > 0)
			{
				if (inFlag == 1) R[nr++] = Q[b];
				b = Advance(b, &ba, m);
			}
			else
			{
				if (inFlag == -1) R[nr++] = P[a];
				a = Advance(a, &aa, n);
			}
		}

		// Quit when both adv, indices have cycled, or one has cycled twice
	}
	while ( ((aa<n) || (ba < m)) && (aa < 2*n) && (ba < 2*m) );

	// deal with special cases (not implemented)
	if (inFlag == 0)
	{
		// The boundaries of P and Q do not cross
		// see if one polygon is inside the other one
		if (ConvexPolyInConvexPoly(P, n, Q, m))
		{
			// P is inside Q so copy P to R
			for (int i=0; i<n; ++i) R[i] = P[i];
			nr = n;
		}
		else if (ConvexPolyInConvexPoly(Q, m, P, n))
		{
			// Q is inside P so copy Q to R
			for (int i=0; i<m; ++i) R[i] = Q[i];
			nr = m;
		}
	}

	return nr;
}

//-----------------------------------------------------------------------------
// This function calculates the intersection of a line and a segment. The line
// is defined by POINT2Ds (a,b) and the segment by (c,d). The intersection
// POINT2D (if it exists) is returned in p. The return code can be used to identify
// the intersection type:
// 0 = no intersection
// 1 = proper intersection
// 2 = vertex intersection (segment coincide at vertex)
// 3 = edge intersections (segment cooincide with line)
int LineSegInt(POINT2D a, POINT2D b, POINT2D c, POINT2D d, POINT2D& p, double eps)
{
	// return code
	int nret = -1;

	// calculate the denominator
	double denom = a.x*(d.y - c.y) + b.x*(c.y - d.y) + d.x*(b.y - a.y) + c.x*(a.y - b.y);

	// if the denominator is zero, the segments are parallel; handle separately
	if (fabs(denom) <= eps) return ParallelInt(a, b, c, d, p, eps);

	double num = a.x*(d.y - c.y) + c.x*(a.y - d.y) + d.x*(c.y - a.y);
	if ((fabs(num) <= eps)||(num==denom)) nret = 2;
	double s = num / denom;

	num = -(a.x*(c.y - b.y) + b.x*(a.y - c.y) + c.x*(b.y - a.y));
	if ((fabs(num) <= eps)||(num==denom)) nret = 2;
	double t = num / denom;

	// check for proper intersection
	if (nret == -1)
	{
		if      ((0.0 < t) && (t < 1.0)) nret = 1;
		else if ((0.0 > t) || (t > 1.0)) nret = 0;
	}

	// find the intersection POINT2D
	p.x = a.x + s*(b.x - a.x);
	p.y = a.y + s*(b.y - a.y);

	assert(nret >= 0);
	return nret;
}

//-----------------------------------------------------------------------------
int ConvexIntersectSH(POINT2D* P, int n, POINT2D* Q, int m, POINT2D* R)
{
	// temporary buffer
	POINT2D tmp[11];

	// copy Q to R
	for (int i=0; i<m; ++i) R[i] = Q[i];
	int nr = m;
	if (nr == 0) return 0;

	const double eps = 1e-7;

	POINT2D p;
	// loop over all the edges of the clipping polygon P
	for (int a=0; a<n; ++a)
	{
		// an edge is defined by the current point and the next one
		int a1 = (a+1)%n;

		// reset the tmp buffer
		int ntmp = 0;

		// see if the last point lies inside the edge A
		// we'll use this to decide if the edges cross
		int b1 = nr-1;
		int b1HA = AreaSign(P[a], P[a1], R[b1], eps);

		// loop over all the input points
		for (int b=0; b<nr; ++b)
		{
			// see if this point lies inside the edge
			int bHA = AreaSign(P[a], P[a1], R[b], eps);
			if (bHA == 0)
			{
				// point b lies on the edge so just add it
				tmp[ntmp++] = R[b];
			}
			else if (bHA > 0)
			{
				// point b lies inside the edge A
				// so check whether an edge was crossed
				if (b1HA < 0)
				{
					int ncode = LineSegInt(P[a],P[a1],R[b1],R[b], p, eps);
					assert((ncode==1)||(ncode==2));
					// add intersection point
					tmp[ntmp++] = p;
				}

				// add point b
				tmp[ntmp++] = R[b];
			}
			else if (bHA < 0)
			{
				// point b lies outside edge, but we still may have an intersection
				if (b1HA > 0)
				{
					int ncode = LineSegInt(P[a],P[a1],R[b1],R[b], p, eps);
					assert((ncode==1)||(ncode==2));
					// add intersection point
					tmp[ntmp++] = p;
				}
			}

			// advance
			b1 = b;
			b1HA = bHA;
		}

		// copy tmp buffer to output buffer
		for (int i=0; i<ntmp; i++) R[i] = tmp[i];
		nr = ntmp;

		// if we did not add any points, the polygons do not intersect
		if (nr == 0) break;
	}

	return nr;
}

double tri_area(vec3d r[3])
{
	return ((r[1]-r[0])^(r[2]-r[0])).norm()*0.5;
}

bool CalculateMortarIntersection(FESurface& ss, FESurface& ms, int k, int l, Patch& patch)
{
	// clear the patch
	patch.Clear();

	// get the surface elements
	FESurfaceElement& es = ss.Element(k);
	FESurfaceElement& em = ms.Element(l);

	// get the nodal coordinates
	const int M = FEElement::MAX_NODES;
	vec3d rs[M], rm[M];
	int ns = es.Nodes(), nm = em.Nodes();
	for (int i=0; i<ns; ++i) rs[i] = ss.Node(es.m_lnode[i]).m_rt;
	for (int i=0; i<nm; ++i) rm[i] = ms.Node(em.m_lnode[i]).m_rt;

	// setup an orthonormal coordinate system
	vec3d c = rs[0];
	vec3d e1 = rs[   1] - rs[0]; e1.unit();
	vec3d e2 = rs[ns-1] - rs[0]; e2.unit();
	vec3d e3 = e1^e2; e3.unit();
	e2 = e3^e1;

	// project all points onto this plane
	POINT2D P[M], Q[M], R[M];
	for (int i=0; i<ns; ++i)
	{
		vec3d r = rs[i] - c;
		vec3d q = r - e3*(r*e3);
		P[i].x = q*e1;
		P[i].y = q*e2;
	}

	// now we do the nodes
	// Note that we loop backwards since the element will
	// in general have opposite winding
	for (int i=0; i<nm; ++i)
	{
		vec3d r = rm[nm-i-1] - c;
		vec3d q = r - e3*(r*e3);
		Q[i].x = q*e1;
		Q[i].y = q*e2;
	}

	// now we calculate the intersection
	int nr = ConvexIntersectSH(P, ns, Q, nm, R);
/*	if (nr > 0)
	{
		for (int i=0; i<nr; ++i)
		{
			POINT2D& r = R[i];
			vec3d x = e1*r.x + e2*r.y + c;

			double qr, qs;
			vec3d ps = ss.ProjectToSurface(es, x, qr, qs);
			const double eps = 1e-4;
			if ((qr < -eps) || (qs < -eps) || (qr+qs > 1+eps))
			{
				// error
				int a = 0;
			}

			vec3d pm = ms.ProjectToSurface(em, x, qr, qs);
			if ((qr < -eps) || (qs < -eps) || (qr+qs > 1+eps))
			{
				// error
				int a = 0;
			}

		}
	}
*/
	if (nr >= 3)
	{
		// evaluate the center of the patch
		POINT2D d;
		d.x = d.y = 0;
		for (int k=0; k<nr; ++k)
		{
			d.x += R[k].x;
			d.y += R[k].y;
		}
		d.x /= nr; d.y /= nr;

		// the center point is always the same
		vec3d r[3];
		r[0] = e1*d.x + e2*d.y + c;

		// calculate the other patch points
		for (int k=0; k<nr; ++k)
		{
			int k1 = (k+1)%nr;
			r[1] = e1*R[k ].x + e2*R[k ].y + c;
			r[2] = e1*R[k1].x + e2*R[k1].y + c;

			patch.Add(r);
		}
	}

	// return
	return (patch.Empty() == false);
}

void CalculateMortarSurface(FESurface& ss, FESurface& ms, MortarSurface& mortar)
{
	// loop over all non-mortar facets
	int NSF = ss.Elements();
	int NMF = ms.Elements();
	for (int i=0; i<NSF; ++i)
	{
		// get the non-mortar surface element
		FESurfaceElement& se = ss.Element(i);

		// loop over all the mortar surface elements
		for (int j=0; j<NMF; ++j)
		{
			// get the next surface element
			FESurfaceElement& me = ms.Element(j);

			// calculate the patch of triangles, representing the intersection
			// of the non-mortar facet with the mortar facet
			Patch patch(i,j);
			CalculateMortarIntersection(ss, ms, i, j, patch);
			mortar.AddPatch(patch);
		}
	}
}

bool ExportMortar(MortarSurface& mortar, const char* szfile)
{
	FILE* fp = fopen(szfile, "wt");
	if (fp == 0) return false;

	fprintf(fp, "solid %s\n", "mortar");
	int NP = mortar.Patches();
	for (int i=0; i<NP; ++i)
	{
		Patch& patch = mortar.GetPatch(i);
		int np = patch.Size();
		for (int j=0; j<np; j++)
		{
			Patch::FACET& tri = patch.Facet(j);

			vec3d e1 = tri.r[1] - tri.r[0];
			vec3d e2 = tri.r[2] - tri.r[0];
			vec3d n = e1^e2; n.unit();

			fprintf(fp, "facet normal %g %g %g\n", n.x, n.y, n.z);
			fprintf(fp, "outer loop\n");
			for (int k=0; k<3; ++k)
			{
				vec3d& r = tri.r[k];
				fprintf(fp, "vertex %g %g %g\n", r.x, r.y, r.z);
			}
			fprintf(fp, "endloop\n");
			fprintf(fp, "endfacet\n");
		}
	}
	fprintf(fp, "endsolid\n");

	fclose(fp);

	return true;
}
