/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

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
    bool            m_Bern;     //!< Use Bernoulli's Relation (Q*|Q|)
    
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
