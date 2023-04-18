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

// onchoidFit.cpp: implementation of the QuadricFit class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "QuadricFit.h"
#include <stdio.h>

using std::swap;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------
// constructor
QuadricFit::QuadricFit()
{
    m_rc = m_c2  = vec3d(0,0,0);
    m_c = 0;
    m_ax[0] = m_ax[1] = m_ax [2] = vec3d(0,0,0);
    m_quad = new Quadric();
    m_eps = 1e-3;
}

//-------------------------------------------------------------------------------
// copy constructor
QuadricFit::QuadricFit(const QuadricFit& qf)
{
    m_rc = qf.m_rc;
    m_c2 = qf.m_c2;
    m_c = qf.m_c;
    m_ax[0] = qf.m_ax[0];
    m_ax[1] = qf.m_ax[1];
    m_ax[2] = qf.m_ax[2];
    m_quad = new Quadric(*qf.m_quad);
    m_eps = qf.m_eps;
}

//-------------------------------------------------------------------------------
// destructor
QuadricFit::~QuadricFit()
{
    
}

//-------------------------------------------------------------------------------
void QuadricFit::GetOrientation()
{
    
    mat3ds C(m_quad->m_c[0], m_quad->m_c[1], m_quad->m_c[2],
             m_quad->m_c[5], m_quad->m_c[3], m_quad->m_c[4]);
    vec3d v = vec3d(m_quad->m_c[6], m_quad->m_c[7], m_quad->m_c[8]);
    double c = m_quad->m_c[9];

    double eval[3];
    vec3d evec[3];
    C.eigen2(eval, evec);   // eigenvalues sorted from smallest to largest
    
    // rectify sign of eigenvalues
    double emag = sqrt(eval[0]*eval[0] + eval[1]*eval[1] + eval[2]*eval[2]);
    if (eval[2] <= emag*m_eps) {
        eval[0] = -eval[0];
        eval[1] = -eval[1];
        eval[2] = -eval[2];
        v = -v;
        c = -c;
    }
    
    // re-sort eigenvalues and eigenvectors from largest to smalles
    if (eval[0] < eval[2]) { swap(eval[0],eval[2]); swap(evec[0], evec[2]); }
    if (eval[0] < eval[1]) { swap(eval[0],eval[1]); swap(evec[0], evec[1]); }
    if (eval[1] < eval[2]) { swap(eval[1],eval[2]); swap(evec[1], evec[2]); }

    // estimate the relative magnitudes of the eigenvalues
    if (fabs(eval[0]) < m_eps*emag) eval[0] = 0;
    if (fabs(eval[1]) < m_eps*emag) eval[1] = 0;
    if (fabs(eval[2]) < m_eps*emag) eval[2] = 0;
    
    // normalize vectors
    evec[0].unit(); evec[1].unit(); evec[2].unit();
    // clean up vector
    if (fabs(evec[0].x) < m_eps) evec[0].x = 0; if (fabs(evec[0].y) < m_eps) evec[0].y = 0; if (fabs(evec[0].z) < m_eps) evec[0].z = 0;
    if (fabs(evec[1].x) < m_eps) evec[1].x = 0; if (fabs(evec[1].y) < m_eps) evec[1].y = 0; if (fabs(evec[1].z) < m_eps) evec[1].z = 0;
    if (fabs(evec[2].x) < m_eps) evec[2].x = 0; if (fabs(evec[2].y) < m_eps) evec[2].y = 0; if (fabs(evec[2].z) < m_eps) evec[2].z = 0;
    evec[0].unit(); evec[1].unit(); evec[2].unit();
    // check if basis is right-handed or not
    double rh = (evec[0] ^ evec[1])*evec[2];
    if (rh < 0) evec[1] = evec[2] ^ evec[0];
    
    // estimate the relative magnitudes of the components of v
    if (fabs(v.x) < m_eps*emag) v.x = 0;
    if (fabs(v.y) < m_eps*emag) v.y = 0;
    if (fabs(v.z) < m_eps*emag) v.z = 0;

    mat3d Q(evec[0].x, evec[1].x, evec[2].x,
            evec[0].y, evec[1].y, evec[2].y,
            evec[0].z, evec[1].z, evec[2].z);
    
    m_ax[0] = evec[0];
    m_ax[1] = evec[1];
    m_ax[2] = evec[2];
    
    m_c2 = vec3d(eval[0], eval[1], eval[2]);
    vec3d d = Q.transpose()*v;
    m_v = v;

    vec3d X0;
    if (m_c2.x*m_c2.y*m_c2.z != 0) {
        X0.x = -d.x/(2*m_c2.x);
        X0.y = -d.y/(2*m_c2.y);
        X0.z = -d.z/(2*m_c2.z);
        m_c = c - m_c2.x*pow(X0.x,2) - m_c2.y*pow(X0.y,2) - m_c2.z*pow(X0.z,2);
        if (fabs(m_c) < m_eps) m_c = 0;
        if (m_c != 0) {
            m_c2 /= fabs(m_c);
            m_c /= fabs(m_c);
        }
    }
    else if (m_c2.x*m_c2.y != 0) {
        X0.x = -d.x/(2*m_c2.x);
        X0.y = -d.y/(2*m_c2.y);
        X0.z = 0;
        m_c = c - m_c2.x*pow(X0.x,2) - m_c2.y*pow(X0.y,2);
        if (m_c != 0) {
            m_c2 /= fabs(m_c);
            m_c /= fabs(m_c);
        }
    }
    else {
        X0.x = -d.x/(2*m_c2.x);
        X0.y = 0;
        X0.z = 0;
        m_c = c - m_c2.x*pow(X0.x,2);
        if (m_c != 0) {
            m_c2 /= fabs(m_c);
            m_c /= fabs(m_c);
        }
    }
    
    m_rc = Q*X0;
    // clean up
    double rc = m_rc.norm();
    if (fabs(m_rc.x) < m_eps*rc) m_rc.x = 0;
    if (fabs(m_rc.y) < m_eps*rc) m_rc.y = 0;
    if (fabs(m_rc.z) < m_eps*rc) m_rc.z = 0;

}

//-------------------------------------------------------------------------------
// complete fit algorithm
//bool QuadricFit::Fit(PointCloud3d* pc)
bool QuadricFit::Fit(std::vector<vec3d>& pc)
{
//    m_quad->SetPointCloud3d(pc);
    m_quad->m_p = pc;
    if (m_quad->GetQuadricCoeficients())
    {
        GetOrientation();
        return true;
    }
    else
        return false;
}

//-------------------------------------------------------------------------------
bool QuadricFit::isSame(const double& a, const double&b)
{
    if (fabs(b-a) < m_eps*fabs(a+b)/2) return true;
    else return false;
}

//-------------------------------------------------------------------------------
QuadricFit::Q_TYPE QuadricFit::GetType()
{
    // determine quadric type
    Q_TYPE type = Q_UNKNOWN;
    double A = m_c2.x, B =  m_c2.y, C = m_c2.z;
    
    if ((A > 0) && (B > 0)) {
        if ((C > 0) && isSame(m_c,-1)) {
            type = Q_ELLIPSOID;
            if (isSame(A,B)) {
                type = Q_SPHEROID;
                if (isSame(B,C)) type = Q_SPHERE;
            }
        }
        else if (C == 0) {
            if (isSame(m_v.z,-1) && (m_c == 0)) {
                type = Q_ELLIPTIC_PARABOLOID;
                if (isSame(A,B))
                    type = Q_CIRCULAR_PARABOLOID;
            }
            else if ((m_v.z == 0) && isSame(m_c,-1)) {
                type = Q_ELLIPTIC_CYLINDER;
                if (isSame(A,B))
                    type = Q_CIRCULAR_CYLINDER;
            }
        }
        else if (C < 0) {
            if  (isSame(m_c,-1)) {
                type = Q_ELLIPTIC_HYPERBOLOID_1;
                if (isSame(A,B)) {
                    type = Q_CIRCULAR_HYPERBOLOID_1;
                }
            }
            else if (isSame(m_c,1)) {
                type = Q_ELLIPTIC_HYPERBOLOID_2;
                if (isSame(A,B)) {
                    type = Q_CIRCULAR_HYPERBOLOID_2;
                }
            }
            else if  (m_c == 0) {
                type = Q_ELLIPTIC_CONE;
                if (isSame(A,B)) {
                    type = Q_CIRCULAR_CONE;
                }
            }
        }
    }
    else if ((A > 0) && (B < 0)) {
        if ((C == 0) && (m_c == 0) && isSame(m_v.z,-1)) {
            type = Q_HYPERBOLIC_PARABOLOID;
        }
        else if ((C == 0) && isSame(m_c,-1) && (m_v.z == 0)) {
            type = Q_HYPERBOLIC_CYLINDER;
        }
    }
    else if ((A > 0) && (B == 0)) {
        type = Q_PARABOLIC_CYLINDER;
    }

    return type;
}

//-------------------------------------------------------------------------------
std::string QuadricFit::GetStringType(Q_TYPE qtype)
{
    std::string quadric;
    switch (qtype) {
        case QuadricFit::Q_ELLIPSOID: quadric = std::string("Ellipsoid"); break;
        case QuadricFit::Q_ELLIPTIC_PARABOLOID: quadric = std::string("Elliptic Paraboloid"); break;
        case QuadricFit::Q_HYPERBOLIC_PARABOLOID: quadric = std::string("Hyperbolic Paraboloid"); break;
        case QuadricFit::Q_ELLIPTIC_HYPERBOLOID_1: quadric = std::string("Elliptic Hyperboloid of One Sheet"); break;
        case QuadricFit::Q_ELLIPTIC_HYPERBOLOID_2: quadric = std::string("Elliptic Hyperboloid of Two Sheets"); break;
        case QuadricFit::Q_ELLIPTIC_CONE: quadric = std::string("Elliptic Cone"); break;
        case QuadricFit::Q_ELLIPTIC_CYLINDER: quadric = std::string("Elliptic Cylinder"); break;
        case QuadricFit::Q_PARABOLIC_CYLINDER: quadric = std::string("Parabolic Cylinder"); break;
        case QuadricFit::Q_SPHEROID: quadric = std::string("Spheroid"); break;
        case QuadricFit::Q_SPHERE: quadric = std::string("Sphere"); break;
        case QuadricFit::Q_CIRCULAR_PARABOLOID: quadric = std::string("Circular Paraboloid"); break;
        case QuadricFit::Q_CIRCULAR_HYPERBOLOID_1: quadric = std::string("Circular Hyperboloid of One Sheet"); break;
        case QuadricFit::Q_CIRCULAR_HYPERBOLOID_2: quadric = std::string("Circular Hyperboloid of Two Sheets"); break;
        case QuadricFit::Q_CIRCULAR_CONE: quadric = std::string("Circular Cone"); break;
        case QuadricFit::Q_CIRCULAR_CYLINDER: quadric = std::string("Circular Cylinder"); break;
        case QuadricFit::Q_UNKNOWN: quadric = std::string("Unknown"); break;
        default: break;
    }
    return quadric;
}
