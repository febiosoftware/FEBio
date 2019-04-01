#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Tangential flow stabilization prescribes a shear traction that opposes
//! tangential fluid velocity on a boundary surface, in the presence of normal
//! flow.  This can help stabilize inflow/outflow conditions.
class FEBIOFLUID_API FETangentialFlowStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FETangentialFlowStabilization(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! initialization
    bool Init() override;
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const double alpha);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const double alpha);
    
protected:
    double			m_beta;     //!< damping coefficient
    double          m_rho;      //!< fluid density
    
    // degrees of freedom
    int     m_dofX, m_dofY, m_dofZ;
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofWXP, m_dofWYP, m_dofWZP;
    
    DECLARE_FECORE_CLASS();
};
