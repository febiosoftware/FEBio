//
//  FEFluidRCRBC.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 3/30/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include "FEFluidMaterial.h"

//-----------------------------------------------------------------------------
//! FEFluidResistanceBC is a fluid surface that has a normal
//! pressure proportional to the flow rate (resistance).
//!
class FEBIOFLUID_API FEFluidRCRBC : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidRCRBC(FEModel* pfem);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override {}
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
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
    double          m_R;        //!< flow resistance
    double          m_Rd;       //!< distal resistance
    double          m_p0;       //!< initial fluid pressure
    double          m_C;        //!< capacitance
    double          m_pd;       //!< downstream pressure
    int             m_Bern;     //!< Use Bernoulli's Relation (Q*|Q|)
    
    vector<double> m_stepHist;  //!< history of time step size (recorded each step)
    vector<double> m_timeHist;  //!< history of time at each step
    vector<double> m_flowHist;  //!< history of flow rate at each step
    
private:
    double              m_alpha;
    double              m_alphaf;
    FEFluidMaterial*    m_pfluid;   //!< pointer to fluid
    
    FEDofList   m_dofW;
    int         m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
