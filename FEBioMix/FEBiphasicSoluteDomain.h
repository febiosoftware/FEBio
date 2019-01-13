#pragma once
#include <vector>
using namespace std;
#include "FEBioMech/FEElasticDomain.h"
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FESolver;

//-----------------------------------------------------------------------------
//! Abstract interface class for biphasic-solute domains.

//! A biphasic domain is used by the biphasic solver.
//! This interface defines the functions that have to be implemented by a
//! biphasic domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEBIOMIX_API FEBiphasicSoluteDomain : public FEElasticDomain
{
public:
    FEBiphasicSoluteDomain(FEModel* pfem);
    virtual ~FEBiphasicSoluteDomain(){}
    
    // --- R E S I D U A L ---
    
    //! internal work for steady-state case
    virtual void InternalForcesSS(FEGlobalVector& R) = 0;
    
    // --- S T I F F N E S S   M A T R I X ---
    
    //! calculates the global stiffness matrix for this domain
    virtual void StiffnessMatrix(FESolver* psolver, bool bsymm) = 0;
    
    //! calculates the global stiffness matrix (steady-state case)
    virtual void StiffnessMatrixSS(FESolver* psolver, bool bsymm) = 0;
    
protected:
    FEBiphasicSolute*   m_pMat;
    int                 m_dofP;		//!< pressure dof index
    int                 m_dofQ;     //!< shell extra pressure dof index
    int                 m_dofC;		//!< concentration dof index
    int                 m_dofD;		//!< shell extra concentration dof index
    int                 m_dofVX;
    int                 m_dofVY;
    int                 m_dofVZ;
};
