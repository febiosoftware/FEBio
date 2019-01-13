//
//  FEBackFlowStabilization.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEBackFlowStabilization_hpp
#define FEBackFlowStabilization_hpp

#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Backflow stabilization prescribes a normal traction that opposes
//! backflow on a boundary surface.
class FEBIOFLUID_API FEBackFlowStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FEBackFlowStabilization(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! initialization
    bool Init() override;
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const FETimeInfo& tp);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const FETimeInfo& tp);
    
protected:
    double			m_beta;     //!< backflow stabilization coefficient
    double          m_rho;      //!< fluid density
    
    // degrees of freedom
    int     m_dofX, m_dofY, m_dofZ;
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofWXP, m_dofWYP, m_dofWZP;
    
    DECLARE_FECORE_CLASS();
};


#endif /* FEBackFlowStabilization_hpp */
