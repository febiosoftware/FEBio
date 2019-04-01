#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! FEFluidResistanceBC is a fluid surface that has a normal
//! pressure proportional to the flow rate (resistance).
//!
class FEBIOFLUID_API FEFluidResistanceBC : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidResistanceBC(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override { m_alpha = tp.alpha; m_alphaf = tp.alphaf; }
    
    //! set the dilatation
    void Update() override;
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;
    
    //! activate
    void Activate() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

private:
    double			m_R;        //!< flow resistance
	double          m_p0;       //!< fluid pressure offset

private:
	double          m_alpha;
    double          m_alphaf;
    FEFluid*        m_pfluid;   //!< pointer to fluid
    
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofWXP, m_dofWYP, m_dofWZP;
    int		m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
