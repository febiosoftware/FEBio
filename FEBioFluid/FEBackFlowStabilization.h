/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
