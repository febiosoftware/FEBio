//
//  FEBiphasicDomain.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEBiphasicDomain_hpp
#define FEBiphasicDomain_hpp

#include <vector>
#include <FEBioMech/FEElasticDomain.h>
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FESolver;

//-----------------------------------------------------------------------------
//! Abstract interface class for biphasic domains.

//! A biphasic domain is used by the biphasic solver.
//! This interface defines the functions that have to be implemented by a
//! biphasic domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEBIOMIX_API FEBiphasicDomain : public FEElasticDomain
{
public:
    FEBiphasicDomain(FEModel* pfem);
    virtual ~FEBiphasicDomain(){}
    
    // --- R E S I D U A L ---
    
    //! internal work for steady-state case
    virtual void InternalForcesSS(FEGlobalVector& R) = 0;
    
    // --- S T I F F N E S S   M A T R I X ---
    
    //! calculates the global stiffness matrix for this domain
    virtual void StiffnessMatrix(FESolver* psolver, bool bsymm) = 0;
    
    //! calculates the global stiffness matrix (steady-state case)
    virtual void StiffnessMatrixSS(FESolver* psolver, bool bsymm) = 0;
    
public: // biphasic domain "properties"
    virtual vec3d FluidFlux(FEMaterialPoint& mp) = 0;
    
protected:
    FEBiphasic*	m_pMat;
    int			m_dofP;		//!< pressure dof index
    int			m_dofQ;     //!< shell extra pressure dof index
    int			m_dofVX;
    int			m_dofVY;
    int			m_dofVZ;
};

#endif /* FEBiphasicDomain_hpp */
