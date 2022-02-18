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
#include "FEDomain2D.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
bool FEDomain2D::Create(int nelems, FE_Element_Spec espec)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEElement2D& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE) 
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
void FEDomain2D::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

//-----------------------------------------------------------------------------
void FEDomain2D::Reset()
{
	ForEachMaterialPoint([](FEMaterialPoint& mp) {
		mp.Init();
	});
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEDomain2D::invjac0(FEElement2D& el, double Ji[2][2], int n)
{
    // initial nodal coordinates
    int neln = el.Nodes();
    double x[FEElement::MAX_NODES];
    double y[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        vec3d& ri = m_pMesh->Node(el.m_node[i]).m_r0;
		x[i] = ri.x;
		y[i] = ri.y;
    }
    
	// calculate Jacobian
	double J[2][2] = {0};
	for (int i=0; i<neln; ++i)
	{
		const double& Gri = el.Hr(n)[i];
		const double& Gsi = el.Hs(n)[i];
		
		J[0][0] += Gri*x[i]; J[0][1] += Gsi*x[i];
		J[1][0] += Gri*y[i]; J[1][1] += Gsi*y[i];
	}
		
	// calculate the determinant
	double det =  J[0][0]*J[1][1] - J[0][1]*J[1][0];
		
	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);

	// calculate the inverse jacobian
	double deti = 1.0 / det;
	Ji[0][0] =  deti*J[1][1];
	Ji[1][0] = -deti*J[1][0];
	Ji[0][1] = -deti*J[0][1];
	Ji[1][1] =  deti*J[0][0];

    return det;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEDomain2D::invjact(FEElement2D& el, double Ji[2][2], int n)
{
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    double x[FEElement::MAX_NODES];
    double y[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        vec3d& ri = m_pMesh->Node(el.m_node[i]).m_rt;
		x[i] = ri.x;
		y[i] = ri.y;
    }
    
	// calculate Jacobian
	double J[2][2] = {0};
	for (int i=0; i<neln; ++i)
	{
		const double& Gri = el.Hr(n)[i];
		const double& Gsi = el.Hs(n)[i];
		
		J[0][0] += Gri*x[i]; J[0][1] += Gsi*x[i];
		J[1][0] += Gri*y[i]; J[1][1] += Gsi*y[i];
	}
		
	// calculate the determinant
	double det =  J[0][0]*J[1][1] - J[0][1]*J[1][0];
		
	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);

	// calculate the inverse jacobian
	double deti = 1.0 / det;
	Ji[0][0] =  deti*J[1][1];
	Ji[1][0] = -deti*J[1][0];
	Ji[0][1] = -deti*J[0][1];
	Ji[1][1] =  deti*J[0][0];    
    return det;
}

//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
vec2d FEDomain2D::gradient(FEElement2D& el, double* fn, int n)
{
    double Ji[2][2];
    invjact(el, Ji, n);
				
    double* Grn = el.Hr(n);
    double* Gsn = el.Hs(n);
    
    double Gx, Gy;
    
    vec2d gradf(0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i];
        
        // calculate gradient
        gradf.x() += Gx*fn[i];
        gradf.y() += Gy*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
vec2d FEDomain2D::gradient(FEElement2D& el, vector<double>& fn, int n)
{
    double Ji[2][2];
    invjact(el, Ji, n);
				
    double* Grn = el.Hr(n);
    double* Gsn = el.Hs(n);
    
    vec2d gradf(0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        double Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i];
        double Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i];
        
        // calculate pressure gradient
        gradf.x() += Gx*fn[i];
        gradf.y() += Gy*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
mat2d FEDomain2D::gradient(FEElement2D& el, vec2d* fn, int n)
{
    double Ji[2][2];
    invjact(el, Ji, n);
				
    vec2d g1(Ji[0][0],Ji[0][1]);
    vec2d g2(Ji[1][0],Ji[1][1]);
    
    double* Gr = el.Hr(n);
    double* Gs = el.Hs(n);
    
    mat2d gradf;
    gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        vec2d tmp = g1*Gr[i] + g2*Gs[i];
        gradf += dyad(fn[i], tmp);
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
mat3d FEDomain2D::gradient(FEElement2D& el, vec3d* fn, int n)
{
    double Ji[2][2];
    invjact(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],0.0);
    vec3d g2(Ji[1][0],Ji[1][1],0.0);
    
    double* Gr = el.Hr(n);
    double* Gs = el.Hs(n);
    
    mat3d gradf;
    gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gradf += fn[i] & (g1*Gr[i] + g2*Gs[i]);
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FEDomain2D::detJ0(FEElement2D &el, int n)
{
    // initial nodal coordinates
    int neln = el.Nodes();
    double x[FEElement::MAX_NODES];
    double y[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        vec3d& ri = m_pMesh->Node(el.m_node[i]).m_r0;
		x[i] = ri.x;
		y[i] = ri.y;
    }
    
	// calculate Jacobian
	double J[2][2] = {0};
	for (int i=0; i<neln; ++i)
	{
		const double& Gri = el.Hr(n)[i];
		const double& Gsi = el.Hs(n)[i];
		J[0][0] += Gri*x[i]; J[0][1] += Gsi*x[i];
		J[1][0] += Gri*y[i]; J[1][1] += Gsi*y[i];
	}
		
	// calculate the determinant
	double det =  J[0][0]*J[1][1] - J[0][1]*J[1][0];
    
    return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to current frame
double FEDomain2D::detJt(FEElement2D &el, int n)
{
    // initial nodal coordinates
    int neln = el.Nodes();
    double x[FEElement::MAX_NODES];
    double y[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        vec3d& ri = m_pMesh->Node(el.m_node[i]).m_rt;
        x[i] = ri.x;
        y[i] = ri.y;
    }
    
    // calculate Jacobian
    double J[2][2] = {0};
    for (int i=0; i<neln; ++i)
    {
        const double& Gri = el.Hr(n)[i];
        const double& Gsi = el.Hs(n)[i];
        J[0][0] += Gri*x[i]; J[0][1] += Gsi*x[i];
        J[1][0] += Gri*y[i]; J[1][1] += Gsi*y[i];
    }
    
    // calculate the determinant
    double det =  J[0][0]*J[1][1] - J[0][1]*J[1][0];
    
    return det;
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant basis vectors of a 2D element
//! at an integration point

void FEDomain2D::CoBaseVectors(FEElement2D& el, int j, vec2d g[2])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hr = el.Hr(j);
    double* Hs = el.Hs(j);
    
    g[0] = g[1] = vec2d(0,0);
    for (int i=0; i<n; ++i)
    {
        vec2d rt = vec2d(m.Node(el.m_node[i]).m_rt.x,m.Node(el.m_node[i]).m_rt.y);
        g[0] += rt*Hr[i];
        g[1] += rt*Hs[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the contravariant basis vectors of a 2D element
//! at an integration point

void FEDomain2D::ContraBaseVectors(FEElement2D& el, int j, vec2d gcnt[2])
{
    vec2d gcov[2];
    CoBaseVectors(el, j, gcov);
    
    mat2d J = mat2d(gcov[0].x(), gcov[1].x(),
                    gcov[0].y(), gcov[1].y());
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec2d(Ji(0,0),Ji(0,1));
    gcnt[1] = vec2d(Ji(1,0),Ji(1,1));
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of covariant basis
//! vectors of a 2D element at an integration point

void FEDomain2D::CoBaseVectorDerivatives(FEElement2D& el, int j, vec2d dg[2][2])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hrr = el.Hrr(j); double* Hrs = el.Hrs(j);
    double* Hsr = el.Hsr(j); double* Hss = el.Hss(j);
    
    dg[0][0] = dg[0][1] = vec2d(0,0);  // derivatives of g[0]
    dg[1][0] = dg[1][1] = vec2d(0,0);  // derivatives of g[1]
    
    for (int i=0; i<n; ++i)
    {
        vec2d rt = vec2d(m.Node(el.m_node[i]).m_rt.x,m.Node(el.m_node[i]).m_rt.y);
        dg[0][0] += rt*Hrr[i]; dg[0][1] += rt*Hsr[i];
        dg[1][0] += rt*Hrs[i]; dg[1][1] += rt*Hss[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of contravariant basis
//! vectors of a solid element at an integration point

void FEDomain2D::ContraBaseVectorDerivatives(FEElement2D& el, int j, vec2d dgcnt[2][2])
{
    vec2d gcnt[2];
    vec2d dgcov[2][2];
    ContraBaseVectors(el, j, gcnt);
    CoBaseVectorDerivatives(el, j, dgcov);
    
    // derivatives of gcnt[0]
    dgcnt[0][0] = -gcnt[0]*(gcnt[0]*dgcov[0][0])-gcnt[1]*(gcnt[0]*dgcov[1][0]);
    dgcnt[0][1] = -gcnt[0]*(gcnt[0]*dgcov[0][1])-gcnt[1]*(gcnt[0]*dgcov[1][1]);
    
    // derivatives of gcnt[1]
    dgcnt[1][0] = -gcnt[0]*(gcnt[1]*dgcov[0][0])-gcnt[1]*(gcnt[1]*dgcov[1][0]);
    dgcnt[1][1] = -gcnt[0]*(gcnt[1]*dgcov[0][1])-gcnt[1]*(gcnt[1]*dgcov[1][1]);
}

//-----------------------------------------------------------------------------
//! calculate the laplacian of a vector function at an integration point
vec2d FEDomain2D::lapvec(FEElement2D& el, vec2d* fn, int n)
{
    vec2d gcnt[2];
    vec2d dgcnt[2][2];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Hr(n);
    double* Gs = el.Hs(n);
    
    // get the shape function derivatives
    double* Hrr = el.Hrr(n); double* Hrs = el.Hrs(n);
    double* Hsr = el.Hsr(n); double* Hss = el.Hss(n);
    
    vec2d lapv(0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        lapv += fn[i]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*gcnt[0]+
                       (gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*gcnt[1]+
                       (gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*gcnt[0]+
                       (gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*gcnt[1]);
    
    return lapv;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of the divergence of a vector function at an integration point
vec2d FEDomain2D::gradivec(FEElement2D& el, vec2d* fn, int n)
{
    vec2d gcnt[2];
    vec2d dgcnt[2][2];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Hr(n);
    double* Gs = el.Hs(n);
    
    // get the shape function derivatives
    double* Hrr = el.Hrr(n); double* Hrs = el.Hrs(n);
    double* Hsr = el.Hsr(n); double* Hss = el.Hss(n);
    
    vec2d gdv(0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gdv +=
        gcnt[0]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*fn[i])+
        gcnt[1]*((gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*fn[i])+
        gcnt[0]*((gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*fn[i])+
        gcnt[1]*((gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*fn[i]);
    
    return gdv;
}
