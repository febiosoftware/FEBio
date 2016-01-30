//
//  FEElasticFergusonShellDomain.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 1/29/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#pragma once
#include "FECore/FEFergusonShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticFergusonShellDomain : public FEFergusonShellDomain, public FEElasticDomain
{
public:
    FEElasticFergusonShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticFergusonShellDomain& operator = (FEElasticFergusonShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }
    
    //! Initialize domain
    bool Initialize(FEModel& fem);
    
    //! Activate the domain
    void Activate();
    
    //! Unpack shell element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat);
    
public: // overrides from FEElasticDomain
    
    //! calculates the residual
    //	void Residual(FESolver* psolver, vector<double>& R);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R);
    
    //! Calculates inertial forces for dynamic problems | todo implement this (removed assert DSR)
    void InertialForces(FEGlobalVector& R, vector<double>& F) { }
    
    //! calculate body force
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf);
    
    // update stresses
    void Update();
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver);
    
    // inertial stiffness \todo implement this (removed assert DSR)
    void MassMatrix(FESolver* psolver, double scale) { }
    
    // body force stiffness \todo implement this (removed assert DSR)
    void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf) { }
    
public:
    
    // --- S T I F F N E S S ---
    
    //! calculates the shell element stiffness matrix
    void ElementStiffness(int iel, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for shell elements
    void ElementInternalForce(FEFergusonShellElement& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
    void ElementBodyForce(FEModel& fem, FEFergusonShellElement& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
    void ElementBodyForce(FEBodyForce& BF, FEFergusonShellElement& el, vector<double>& fe);
    
protected:
    FESolidMaterial*	m_pMat;
};
