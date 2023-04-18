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
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! This surface load represents the traction applied on the solid at the
//! interface between a fluid and solid in an FSI analysis.
class FEBIOFLUID_API FEBiphasicFSITraction : public FESurfaceLoad
{
public:
    //! constructor
    FEBiphasicFSITraction(FEModel* pfem);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! initialization
    bool Init() override;
    
private:
    double GetFluidDilatation(FESurfaceMaterialPoint& mp, double alpha);
    mat3ds GetFluidStress(FESurfaceMaterialPoint& mp);
    
protected:
    vector<double>      m_s;        //!< scale factor
    vector<FEElement*>  m_elem;     //!< list of fluid-FSI elements
    
    // degrees of freedom
    FEDofList    m_dofU, m_dofSU, m_dofW;
    int        m_dofEF;
    
protected:
    bool                m_bshellb;  //!< flag for prescribing traction on shell bottom
    
    DECLARE_FECORE_CLASS();
};
