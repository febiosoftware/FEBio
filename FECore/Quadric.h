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

#pragma once
#include <vector>
#include "vec3d.h"
#include "fecore_api.h"

class FECORE_API Quadric
{
public:
    Quadric() {}
    Quadric(std::vector<vec3d>& p) { m_p = p; }
    Quadric(Quadric* q);
    Quadric(const Quadric& q);
    ~Quadric();
    
public:
    // assign a point cloud to this bivariate spline object
    void SetPoints(std::vector<vec3d>& p) { m_p = p; }
    
    // fit the point cloud to get quadric surface coefficients
    bool GetQuadricCoeficients();
    
    // Evaluate the surface normal at point p
    vec3d SurfaceNormal(const vec3d p);
    
    // Evaluate the surface principal curvatures kappa and directions v at point p
    void SurfaceCurvature(const vec3d p, const vec3d n, vec2d& kappa, vec3d* v);
    
    // Find ray-quadric surface intersections x: p is point on ray, n is normal along ray
    void RayQuadricIntersection(const vec3d p, const vec3d n, std::vector<vec3d>* x, std::vector<double>* t = nullptr);
    
    // Find the point on the quadric closest to the point p
    vec3d ClosestPoint(const vec3d p);
    
    // Find the point on the quadric closest to the point p
    vec3d ClosestPoint(const vec3d p, const vec3d norm);

public:
    double          m_c[10];    // quadric surface coefficients
    std::vector<vec3d>  m_p;    // point coordinates
};
