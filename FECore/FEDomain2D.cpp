//
//  FEDomain2D.cpp
//  FECore
//
//  Created by Gerard Ateshian on 12/17/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEDomain2D.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void FEDomain2D::InitElements()
{
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEElement2D& el = m_Elem[i];
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Init(false);
    }
}

//-----------------------------------------------------------------------------
void FEDomain2D::Reset()
{
    for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
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
	if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);

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
	if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);

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
        gradf += dyad(fn[i], g1*Gr[i] + g2*Gs[i]);
    
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
void FEDomain2D::ShallowCopy(DumpStream& dmp, bool bsave)
{
    int NEL = (int) m_Elem.size();
    for (int i=0; i<NEL; ++i)
    {
        FEElement2D& el = m_Elem[i];
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->ShallowCopy(dmp, bsave);
    }
}

//-----------------------------------------------------------------------------
void FEDomain2D::Serialize(DumpFile &ar)
{
    if (ar.IsSaving())
    {
        ar << m_Node;
        
        for (size_t i=0; i<m_Elem.size(); ++i)
        {
            FEElement2D& el = m_Elem[i];
            ar << el.Type();
            
            ar << el.GetMatID();
            ar << el.m_nID;
            ar << el.m_node;
			ar << el.m_lnode;
            
            for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
        }
    }
    else
    {
        int n, mat;
        
        ar >> m_Node;
        
        FEModel& fem = *ar.GetFEModel();
        FEMaterial* pmat = GetMaterial();
        assert(pmat);
        
        for (size_t i=0; i<m_Elem.size(); ++i)
        {
            FEElement2D& el = m_Elem[i];
            ar >> n;
            
            el.SetType(n);
            
            ar >> mat; el.SetMatID(mat);
            ar >> el.m_nID;
            ar >> el.m_node;
			ar >> el.m_lnode;
            
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
                el.GetMaterialPoint(j)->Serialize(ar);
            }
        }
    }
}
