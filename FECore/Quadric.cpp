/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

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

#include "Quadric.h"
#include "matrix.h"
#include <math.h>
using namespace std;
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

//-------------------------------------------------------------------------------
// constructor
Quadric::Quadric(Quadric* q)
{
    for (int i=0; i<10; ++i) m_c[i] = q->m_c[i];
    m_p = q->m_p;
}

//-------------------------------------------------------------------------------
// copy constructor
Quadric::Quadric(const Quadric& q)
{
    for (int i=0; i<10; ++i) m_c[i] = q.m_c[i];
    m_p = q.m_p;
}

//-------------------------------------------------------------------------------
// destructor
Quadric::~Quadric()
{
    for (int i=0; i<10; ++i) m_c[i] = 0;
    m_p.clear();
}

//--------------------------------------------------------------------------------------
bool Quadric::GetQuadricCoeficients()
{
    // get number of points in point cloud
    int N = (int)m_p.size();

    matrix A(10,10);
    matrix d(1,10);
    
    //Populate D
    A.zero();
    for (int j = 0; j < N; ++j)
    {
        vec3d p = m_p[j];
        d(0,0) = p.x*p.x;
        d(0,1) = p.y*p.y;
        d(0,2) = p.z*p.z;
        d(0,3) = p.y*p.z;
        d(0,4) = p.x*p.z;
        d(0,5) = p.x*p.y;
        d(0,6) = p.x;
        d(0,7) = p.y;
        d(0,8) = p.z;
        d(0,9) = 1;
        matrix D(10,10);
        D = d.transpose()*d;
        A += D;
    }

    // find eigenvalues and eigenvectors of A
    vector<double> Aval;
    matrix Avec(10,10);
    bool good = A.eigen_vectors(Avec, Aval);

    if (good) {
        // find smallest eigenvalue
        int imin=0;
        for (int i=1; i<10; ++i)
            if (fabs(Aval[i]) < fabs(Aval[imin])) imin = i;
            //if (Aval[i] < Aval[imin]) imin = i;
        
        // store the coefficients of the quadric surface
        for (int j=0; j<10; ++j)
            m_c[j] = Avec(j,imin);
        
        return true;
    }
    
    return false;
}

//--------------------------------------------------------------------------------------
// Evaluate the surface normal at point p
vec3d Quadric::SurfaceNormal(const vec3d p)
{
    // extract quadric surface coefficients
    // F(x,y,z) = a x^2 + b y^2 + c z^2 + e y z + f z x + g x y + l x + m y + n z + d = 0
    double a = m_c[0];
    double b = m_c[1];
    double c = m_c[2];
    double e = m_c[3];
    double f = m_c[4];
    double g = m_c[5];
    double l = m_c[6];
    double m = m_c[7];
    double n = m_c[8];
    
    // evaluate first derivatives
    double Fx = 2*a*p.x + f*p.z + g*p.y + l;
    double Fy = 2*b*p.y + e*p.z + g*p.x + m;
    double Fz = 2*c*p.z + e*p.y + f*p.x + n;
    
    vec3d xn(Fx,Fy,Fz);
    xn.unit();
    
    return xn;
}

//--------------------------------------------------------------------------------------
void Quadric::SurfaceCurvature(const vec3d p, const vec3d pn, vec2d& kappa, vec3d* v)
{
    double eps = 1e-9;
    
    // find magnitude of quadric gradient components
    // to determine if we need to do permutation of axes
    double f1 = fabs(2*m_c[0]*p.x + 2*m_c[4]*p.z + 2*m_c[5]*p.y + 2*m_c[6]);
    double f2 = fabs(2*m_c[1]*p.y + 2*m_c[3]*p.z + 2*m_c[5]*p.x + 2*m_c[7]);
    double f3 = fabs(2*m_c[2]*p.z + 2*m_c[3]*p.y + 2*m_c[4]*p.x + 2*m_c[8]);
    
    // pick direction with maximum gradient component magnitude
    int e1 = 0, e2 = 1, e3 = 2;
    double vmax = fmax(fmax(f1,f2),f3);
    if (vmax == f1) { e3 = 0; e1 = 1; e2 = 2; }
    else if (vmax == f2) { e3 = 1; e1 = 2; e2 = 0; }
    
    // extract quadric surface coefficients (and clean up roundoff errors)
    // F(x,y,z) = a x^2 + b y^2 + c z^2 + e y z + f z x + g x y + l x + m y + n z + d = 0
    double a = (fabs(m_c[e1]) > eps) ? m_c[e1] : 0;
    double b = (fabs(m_c[e2]) > eps) ? m_c[e2] : 0;
    double c = (fabs(m_c[e3]) > eps) ? m_c[e3] : 0;
    double e = (fabs(m_c[e1+3]) > eps) ? m_c[e1+3] : 0;
    double f = (fabs(m_c[e2+3]) > eps) ? m_c[e2+3] : 0;
    double g = (fabs(m_c[e3+3]) > eps) ? m_c[e3+3] : 0;
    double l = (fabs(m_c[e1+6]) > eps) ? m_c[e1+6] : 0;
    double m = (fabs(m_c[e2+6]) > eps) ? m_c[e2+6] : 0;
    double n = (fabs(m_c[e3+6]) > eps) ? m_c[e3+6] : 0;
    
    // evaluate quadric surface derivatives and normal
    double x = p(e1), y = p(e2), z = p(e3);
    double den = n + f*x + e*y + 2*c*z;
    double zx = -(l + 2*a*x + g*y + f*z)/den;
    double zy = -(m + g*x + 2*b*y + e*z)/den;
    double zxx = -2*(a + (f + c*zx)*zx)/den;
    double zyy = -2*(b + (e + c*zy)*zy)/den;
    double zxy = -(g + e*zx + f*zy + 2*c*zx*zy)/den;
    vec3d xu, xv;
    xu(e1) = 1; xu(e2) = 0; xu(e3) = zx;
    xv(e1) = 0; xv(e2) = 1; xv(e3) = zy;
    vec3d xuu, xvv, xuv;
    xuu(e1) = 0; xuu(e2) = 0; xuu(e3) = zxx;
    xvv(e1) = 0; xvv(e2) = 0; xvv(e3) = zyy;
    xuv(e1) = 0; xuv(e2) = 0; xuv(e3) = zxy;
    vec3d xn = (xu ^ xv).normalized();
    
    // get coefficients of fundamental forms
    double E = xu*xu, F = xu*xv, G = xv*xv;
    double L = xuu*xn, M = xuv*xn, N = xvv*xn;

    // evaluate mean and gaussian curvatures
    double tmp = E*G - F*F;
    // gaussian
    double kg = (L*N - M*M)/tmp;
    // mean
    double km = (2*F*M - E*N - G*L)/2./tmp;
    // evaluate principal curvatures
    double d1 = sqrt(km*km - kg);
    // max curvature
    double kmax = km + d1;
    double kmin = km - d1;

    // evaluate principal directions of curvature
    double a2 = F*N - G*M, b2 = E*N - G*L, c2 = E*M - F*L;
    double d2 = b2*b2 - 4*a2*c2; if (fabs(d2) < 0) d2 = 0;
    d2 = sqrt(d2);  // d2=0 represents an umbilical point
    double thmax = atan2(-b2 + d2, 2*a2);
    double thmin = (d2 != 0) ? atan2(-b2 - d2, 2*a2) : thmax + M_PI/2;
    vec3d xmax = xu*cos(thmax) + xv*sin(thmax); xmax.unit();
    vec3d xmin = xu*cos(thmin) + xv*sin(thmin); xmin.unit();
    
    // check quadric normal versus face normal
    if (xn*pn >= 0) {
        kappa.x() = kmax; kappa.y() = kmin;
        v[0] = xmax; v[1] = xmin;
    }
    else {
        kappa.x() = kmin; kappa.y() = kmax;
        v[0] = xmin; v[1] = xmax;
    }
    // fix handedness if neeeded
    if ((v[0]^v[1])*pn < 0) v[0] = -v[0];
}

//--------------------------------------------------------------------------------------
// Find ray-quadric surface intersections x: p is point on ray, n is normal along ray
// There are three possible solutions: 0 roots, 1 root, and 2 roots
// When 2 roots are found, sort results from closest to farthest
void Quadric::RayQuadricIntersection(const vec3d p, const vec3d n, vector<vec3d>* x, vector<double>* t)
{
    double a = n.x*n.x*m_c[0]+n.y*n.y*m_c[1]+n.z*n.z*m_c[2]+n.y*n.z*m_c[3]+n.x*n.z*m_c[4]+n.x*n.y*m_c[5];
    double b = n.x*(2*p.x*m_c[0]+p.z*m_c[4]+p.y*m_c[5]+m_c[6])
    +n.y*(2*p.y*m_c[1]+p.z*m_c[3]+p.x*m_c[5]+m_c[7])
    +n.z*(2*p.z*m_c[2]+p.y*m_c[3]+p.x*m_c[4]+m_c[8]);
    double c = p.x*p.x*m_c[0]+p.y*p.y*m_c[1]+p.z*p.z*m_c[2]
    +p.x*(p.z*m_c[4]+p.y*m_c[5]+m_c[6])
    +p.y*(p.z*m_c[3]+m_c[7])+p.z*m_c[8]+m_c[9];
    
    x->clear();
    t->clear();
    double d = b*b - 4*a*c;
    if (d < 0) return;
    else if ((a == 0) && (b != 0)) {
        double t1 = -c/b;
        t->push_back(t1);
        x->push_back(p + n*t1);
        return;
    }
    else if ((d == 0) && (a != 0)) {
        double t1 = -b/(2*a);
        t->push_back(t1);
        x->push_back(p + n*t1);
        return;
    }
    else if (a != 0)  {
        d = sqrt(d);
        double t1 = (-b - d)/(2*a);
        double t2 = (-b + d)/(2*a);
        if (fabs(t1) < fabs(t2)) {
            t->push_back(t1);
            t->push_back(t2);
            x->push_back(p + n*t1);
            x->push_back(p + n*t2);
        }
        else {
            t->push_back(t2);
            t->push_back(t1);
            x->push_back(p + n*t2);
            x->push_back(p + n*t1);
        }
    }
}

//--------------------------------------------------------------------------------------
// This routine finds a closest point approximation (not the exact solution)
vec3d Quadric::ClosestPoint(const vec3d p)
{
    vector<vec3d> xsol;
    vector<double> tsol;

    vec3d n1 = vec3d(1,0,0);
    vec3d n2 = vec3d(0,1,0);
    vec3d n3 = vec3d(0,0,1);
    vector<vec3d> x1, x2, x3;
    vector<double> t1, t2, t3;
    RayQuadricIntersection(p, n1, &x1, &t1);
    RayQuadricIntersection(p, n2, &x2, &t2);
    RayQuadricIntersection(p, n3, &x3, &t3);
    if (t1.size() > 0) {
        xsol.push_back(x1[0]);
        tsol.push_back(t1[0]);
    }
    if (t2.size() > 0) {
        xsol.push_back(x2[0]);
        tsol.push_back(t2[0]);
    }
    if (t3.size() > 0) {
        xsol.push_back(x3[0]);
        tsol.push_back(t3[0]);
    }

    int N = (int)tsol.size();
    if (N > 0) {
        int imin = 0;
        double tmin = fabs(tsol[imin]);
        for (int i = 1; i<N; ++i) {
            if (fabs(tsol[i]) < tmin) {
                imin = i;
                tmin = fabs(tsol[i]);
            }
        }
        return xsol[imin];
    }
    return p;
}

//--------------------------------------------------------------------------------------
// This routine finds a closest point approximation (use norm of face)
vec3d Quadric::ClosestPoint(const vec3d p, const vec3d norm)
{
    vec3d xsol;
    double tsol;
    
    vector<vec3d> x;
    vector<double> t;
    RayQuadricIntersection(p, norm, &x, &t);

    if (t.size() > 0) {
        xsol = x[0];
        tsol = t[0];
        
        return xsol;
    }
    else {
        ClosestPoint(p);
    }
    return p;
}
