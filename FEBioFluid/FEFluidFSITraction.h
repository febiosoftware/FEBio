#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! This surface load represents the traction applied on the solid at the
//! interface between a fluid and solid in an FSI analysis.
class FEBIOFLUID_API FEFluidFSITraction : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidFSITraction(FEModel* pfem);
    
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
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const FETimeInfo& tp, const int iel);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const FETimeInfo& tp, const int iel);
    
protected:
	vector<double>      m_K;        //!< fluid bulk modulus
	vector<double>      m_s;        //!< scale factor
	vector<bool>        m_bself;    //!< flag if fluid pressure is applied on its own FSI mesh
	vector<FEElement*>  m_elem;     //!< list of fluid-FSI elements

    // degrees of freedom
    int		m_dofX, m_dofY, m_dofZ;
    int     m_dofSX, m_dofSY, m_dofSZ;
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofEF, m_dofEFP;
};
