//
//  FEFergusonShellDomain.hpp
//  FECore
//
//  Created by Gerard Ateshian on 1/29/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FEFergusonShellDomain : public FEDomain
{
public:
    //! constructor
    FEFergusonShellDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_FERGUSON, pm) {}
    
    //! create storage for elements
    void create(int nsize) { m_Elem.resize(nsize); }
    
    //! return nr of elements
    int Elements() { return (int)m_Elem.size(); }
    
    //! element access
    FEFergusonShellElement& Element(int n) { return m_Elem[n]; }
    FEElement& ElementRef(int n) { return m_Elem[n]; }
    
    int GetElementType() { return m_Elem[0].Type(); }
    
    //! Initialize elements
    void InitElements();
    
    //! Reset element data
    void Reset();
    
    // inverse jacobian with respect to reference frame
    double invjac0(FEFergusonShellElement& el, double J[3][3], int n);
    
    // jacobian with respect to reference frame
    double detJ0(FEFergusonShellElement& el, int n);
    
    //! calculates referential covariant basis vectors at an integration point
    void CoBaseVectors0(FEFergusonShellElement& el, int j, vec3d g[3]);
    
    //! calculates referential contravariant basis vectors at an integration point
    void ContraBaseVectors0(FEFergusonShellElement& el, int j, vec3d g[3]);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEFergusonShellElement& el, int j, vec3d g[3]);
    
    //! calculates mid-shell surface data at an integration point
    void MidShellSurface(FEFergusonShellElement& el, int n,
                         double& L, double& Lr, double& Ls,
                         double& H, double& Hr, double& Hs,
                         vec3d& ghr, vec3d& ghs,
                         vec3d& ghrr, vec3d& ghrs, vec3d& ghss,
                          vec3d& t,vec3d& tr, vec3d& ts);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEFergusonShellElement& el, int j, vec3d g[3]);
    
    //! calculate jacobian in current frame
    double detJt(FEFergusonShellElement& el, int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEFergusonShellElement& el, double J[3][3], int n);
    
    // calculate deformation gradient
    double defgrad(FEFergusonShellElement& el, mat3d& F, int n);
    
    //! calculates nodal covariant basis vectors in reference frame
    void NodalCoBaseVectors0(FEFergusonShellElement& el, vec3d* r0, vec3d* gr, vec3d* gs, vec3d* t);
    
    //! calculates nodal covariant basis vectors in current frame
    void NodalCoBaseVectors(FEFergusonShellElement& el, vec3d* rt, vec3d* gr, vec3d* gs, vec3d* t, double* lam);
    
    //! calculates nodal covariant basis torsional tangents in current frame
    void NodalCoBaseTorsionTangents(FEFergusonShellElement& el, vec3d* rt, vec3d* t, double* lam, mat3d* Gr, mat3d* Gs);
    
    //! Serialize domain data to archive
    void Serialize(DumpStream& ar);
    
protected:
    vector<FEFergusonShellElement>	m_Elem;	//!< array of elements
    int					m_dofU;
    int					m_dofV;
    int					m_dofW;
};
