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
#include "febiofluid_api.h"
#include "FEFluidSolutes.h"

//-----------------------------------------------------------------------------
//! Backflow stabilization prescribes a normal traction that opposes
//! backflow on a boundary surface.
class FEBIOFLUID_API FEFluidSolutesNaturalFlux : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidSolutesNaturalFlux(FEModel* pfem);
    
    //! Initialization
    bool Init() override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    void SetSolute(int isol) { m_isol = isol; }
    
    //! calculate flux stiffness
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    //! calculate residual
    void LoadVector(FEGlobalVector& R) override;
    
    //! update
    void Update() override;
    
protected:
    int         m_isol;      //!< solute index
    bool        m_bup;          //!< flag to call Update function
    
    // degrees of freedom
    FEDofList    m_dofW;
    FEDofList    m_dofEF;
    FEDofList    m_dofC;
    
    DECLARE_FECORE_CLASS();
};
