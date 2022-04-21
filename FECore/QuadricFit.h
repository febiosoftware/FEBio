/*This file is part of the FEBio Studio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio-Studio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include <string>
#include "Quadric.h"
#include "quatd.h"
#include "fecore_api.h"

class FECORE_API QuadricFit
{
public:
    enum Q_TYPE { Q_ELLIPSOID, Q_ELLIPTIC_PARABOLOID, Q_HYPERBOLIC_PARABOLOID,
        Q_ELLIPTIC_HYPERBOLOID_1, Q_ELLIPTIC_HYPERBOLOID_2, Q_ELLIPTIC_CONE,
        Q_ELLIPTIC_CYLINDER, Q_HYPERBOLIC_CYLINDER, Q_PARABOLIC_CYLINDER,
        Q_SPHEROID, Q_SPHERE, Q_CIRCULAR_PARABOLOID, Q_CIRCULAR_HYPERBOLOID_1,
        Q_CIRCULAR_HYPERBOLOID_2, Q_CIRCULAR_CONE, Q_CIRCULAR_CYLINDER, Q_UNKNOWN
    };

public:
    QuadricFit();
    QuadricFit(const QuadricFit& qf);
    virtual ~QuadricFit();
    
    bool Fit(std::vector<vec3d>& pc);
    Q_TYPE GetType();
    std::string GetStringType(Q_TYPE qtype);

protected:
    vec3d Transform(vec3d& rc, quatd& q, const vec3d& p)
    {
        vec3d r = p - rc;
        q.RotateVector(r);
        return r;
    }
    
    void GetOrientation();
    bool isSame(const double& a, const double&b);
    

public:
    vec3d   m_rc;   // center of quadric
    vec3d   m_ax[3];// quadric axes
    vec3d   m_c2;   // coefficients of square terms
    vec3d   m_v;    // coefficients of linear terms
    double  m_c;    // constant
    double  m_eps;  // tolerance
    
    Quadric*    m_quad; // quadric object
};
