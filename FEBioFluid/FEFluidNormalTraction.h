#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelParam.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FEFluidNormalTraction is a fluid surface that has a normal
//! viscous traction prescribed on it.
//!
class FEBIOFLUID_API FEFluidNormalTraction : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalTraction(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
private:
    FEParamDouble	m_traction;	//!< magnitude of traction load
    
    int		m_dofWX;
    int		m_dofWY;
    int		m_dofWZ;
    
    DECLARE_FECORE_CLASS();
};
