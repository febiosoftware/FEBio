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
//! A 2D element is assumed to have no variation through the thickness.
double FEDomain2D::invjac0(FEElement2D& el, double Ji[3][3], int n)
{
    // number of nodes
    int neln = el.Nodes();
    
    // initial nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    }
    
    // evaluate covariant basis vectors
    vec3d g1, g2, g3;
    g1 = g2 = vec3d(0,0,0);
    double J[3][3];
    for (int i=0; i<neln; ++i)
    {
        const double& Hri = el.Hr(n)[i];
        const double& Hsi = el.Hs(n)[i];
        
        // in-plane covariant basis vectors
        g1 += r0[i]*Hri;
        g2 += r0[i]*Hsi;
    }
    // normal vector
    g3 = g1 ^ g2; g3.unit();
    
    // calculate jacobian
    J[0][0] = g1.x; J[0][1] = g2.x; J[0][2] = g3.x;
    J[1][0] = g1.y; J[1][1] = g2.y; J[1][2] = g3.y;
    J[2][0] = g1.z; J[2][1] = g2.z; J[2][2] = g3.z;
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);
    
    // calculate the inverse of the jacobian
    double deti = 1.0 / det;
    
    Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    
    Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
    Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    return det;
}


//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
//! A 2D element is assumed to have no variation through the thickness.
double FEDomain2D::invjact(FEElement2D& el, double Ji[3][3], int n)
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    // initial nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    }
    
    // evaluate covariant basis vectors
    vec3d g1, g2, g3;
    g1 = g2 = vec3d(0,0,0);
    double J[3][3];
    for (i=0; i<neln; ++i)
    {
        const double& Hri = el.Hr(n)[i];
        const double& Hsi = el.Hs(n)[i];
        
        // in-plane covariant basis vectors
        g1 += r0[i]*Hri;
        g2 += r0[i]*Hsi;
    }
    // normal vector
    g3 = g1 ^ g2; g3.unit();
    
    // calculate jacobian
    J[0][0] = g1.x; J[0][1] = g2.x; J[0][2] = g3.x;
    J[1][0] = g1.y; J[1][1] = g2.y; J[1][2] = g3.y;
    J[2][0] = g1.z; J[2][1] = g2.z; J[2][2] = g3.z;
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.m_nID, n+1, det);
    
    // calculate the inverse of the jacobian
    double deti = 1.0 / det;
    
    Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    
    Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
    Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    return det;
}


//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
vec3d FEDomain2D::gradient(FEElement2D& el, double* fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    double* Grn = el.Hr(n);
    double* Gsn = el.Hs(n);
    
    double Gx, Gy, Gz;
    
    vec3d gradf;
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i];
        Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i];
        
        // calculate pressure gradient
        gradf.x += Gx*fn[i];
        gradf.y += Gy*fn[i];
        gradf.z += Gz*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
vec3d FEDomain2D::gradient(FEElement2D& el, vector<double>& fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    double* Grn = el.Hr(n);
    double* Gsn = el.Hs(n);
    
    double Gx, Gy, Gz;
    
    vec3d gradf;
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i];
        Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i];
        
        // calculate pressure gradient
        gradf.x += Gx*fn[i];
        gradf.y += Gy*fn[i];
        gradf.z += Gz*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
//! A 2D element is assumed to have no variation through the thickness.
mat3d FEDomain2D::gradient(FEElement2D& el, vec3d* fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
    vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
    vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
    
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
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    // initial nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    }
    
    // evaluate covariant basis vectors
    vec3d g1, g2, g3;
    g1 = g2 = vec3d(0,0,0);
    double J[3][3];
    for (i=0; i<neln; ++i)
    {
        const double& Hri = el.Hr(n)[i];
        const double& Hsi = el.Hs(n)[i];
        
        // in-plane covariant basis vectors
        g1 += r0[i]*Hri;
        g2 += r0[i]*Hsi;
    }
    // normal vector
    g3 = g1 ^ g2; g3.unit();
    
    // calculate jacobian
    J[0][0] = g1.x; J[0][1] = g2.x; J[0][2] = g3.x;
    J[1][0] = g1.y; J[1][1] = g2.y; J[1][2] = g3.y;
    J[2][0] = g1.z; J[2][1] = g2.z; J[2][2] = g3.z;
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
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
