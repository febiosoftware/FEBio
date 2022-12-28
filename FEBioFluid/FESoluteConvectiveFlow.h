/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEOctreeSearch.h>
#include "FECore/FENormalProjection.h"
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FESoluteConvectiveFlow is a fluid domain where solute is transported
//! purely via convection.  The solute concentration can be obtained from
//! the fluid velocity and dilatation results.
// TODO: This is not a surface load!
class FEBIOFLUID_API FESoluteConvectiveFlow : public FESurfaceLoad
{
public:
    //! constructor
    FESoluteConvectiveFlow(FEModel* pfem);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS) override {}
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R) override {};
    
    //! set the dilatation
    void Update() override;
    
    //! initialize
    bool Init() override;
    
    //! activate
    void Activate() override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;

private:
    int         m_sol;      //!< solute id
    FEDofList   m_dofW;
    int         m_dofEF;
    int         m_dofC;
    vector<bool>    m_bexclude;
    FENormalProjection* m_np;
    FEOctreeSearch* m_octree;
    
    DECLARE_FECORE_CLASS();
};
