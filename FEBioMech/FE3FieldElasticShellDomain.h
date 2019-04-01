#pragma once
#include "FEElasticShellDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! shell element. Results indicate that using this class produces poorer convergence
//! with shells than the standard FEElasticShellDomain.  This class is included
//! only for development purposes.
class FE3FieldElasticShellDomain : public FEElasticShellDomain
{
protected:
    struct ELEM_DATA
    {
        double    eJ;        // average element jacobian
        double    ep;        // average pressure
        double    Lk;        // Lagrangian multiplier

		void Serialize(DumpStream& ar);
    };
    
public:
    //! constructor
    FE3FieldElasticShellDomain(FEModel* pfem) : FEElasticShellDomain(pfem) {}
    
    //! \todo Do I really use this?
    FE3FieldElasticShellDomain& operator = (FE3FieldElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }
    
    //! initialize class
	bool Init() override;
    
    //! Reset data
    void Reset() override;
    
    //! augmentation
    bool Augment(int naug);
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
public: // overridden from FEElasticDomain
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // calculate stiffness matrix
    void StiffnessMatrix(FESolver* psolver) override;
    
protected:
    //! Dilatational stiffness component for nearly-incompressible materials
    void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke);
    
    //! material and geometrical stiffness components
    void ElementStiffness(int iel, matrix& ke);
    
    //! update the stress of an element
    void UpdateElementStress(int iel);
    
protected:
    vector<ELEM_DATA>    m_Data;
};
