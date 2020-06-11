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
#include "FEThermoFluid.h"
#include <FEBioMech/FEAugLagLinearConstraint.h>
#include <FECore/FENodeSet.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEConstraintNormalFlow class implements a fluid surface with zero
//! tangential velocity as a linear constraint.

class FEBIOFLUID_API FEThermoFluidPressureLoad : public FENodeConstraintSet
{
public:
    //! constructor
    FEThermoFluidPressureLoad(FEModel* pfem);
    
    //! destructor
    ~FEThermoFluidPressureLoad() {}
    
    //! Activation
    void Activate() override;
    
    //! initialization
    bool Init() override;
    
    //! Get the surface
    FENodeSet* GetNodeSet() override { return &m_nset; }
    
protected:
    FENodeSet       m_nset;
    int             m_dofT;
    int             m_dofEF;
    double          m_alpha;
    FEThermoFluid*  m_tfluid;
    
    double constraint(FEAugLagLinearConstraint& LC) override;

public:
    double  m_p0;       // prescribed pressure
    int     m_matID;    // material associated with node set
    
    DECLARE_FECORE_CLASS();
};
