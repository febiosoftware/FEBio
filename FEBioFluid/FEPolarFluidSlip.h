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
#include <FECore/FEModelParam.h>
#include <FECore/FEPrescribedBC.h>
#include "febiofluid_api.h"

class FEFluidMaterial;

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FEPolarFluidSlip : public FEPrescribedSurface
{
public:
    //! constructor
    FEPolarFluidSlip(FEModel* pfem);
    
    //! set the dilatation
    void Update() override;
    void UpdateModel() override;
    
    //! initialize
    bool Init() override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;
    
    void PrepStep(std::vector<double>& ui, bool brel) override;
    
    // return the value for node i, dof j
    void GetNodalValues(int nodelid, std::vector<double>& val) override;
    
    // copy data from another class
    void CopyFrom(FEBoundaryCondition* pbc) override;
    
private:
    void UpdateAngularVelocity();
    
protected:
    FEDofList       m_dofW, m_dofG;
    FESurface*      m_psurf;
    
public:
    double          m_m0;       //!< threshold of couple traction
    double          m_ksi;      //!< slope of linear angular velocity-couple traction variation
    bool            m_brel;     //!< relative flag
    vector<vec3d>   m_v;        //!< nodal values of linear velocity
    vector<vec3d>   m_g;        //!< nodal values of angular velocity
    
    DECLARE_FECORE_CLASS();
};
