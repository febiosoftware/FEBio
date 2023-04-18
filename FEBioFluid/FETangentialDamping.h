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
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Tangential damping prescribes a shear traction that opposes tangential
//! fluid velocity on a boundary surface.  This can help stabilize inflow
//! conditions.
class FEBIOFLUID_API FETangentialDamping : public FESurfaceLoad
{
public:
    //! constructor
    FETangentialDamping(FEModel* pfem);
    
    //! Initialization
    bool Init() override;
    
    //! data serialization
    void Serialize(DumpStream& ar) override;

    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R) override;
    
protected:
	vec3d FluidVelocity(FESurfaceMaterialPoint& mp, double alpha);

protected:
    double			m_eps;      //!< damping coefficient (penalty)
    
    // degrees of freedom
	FEDofList	m_dofW;
    
    DECLARE_FECORE_CLASS();
};
