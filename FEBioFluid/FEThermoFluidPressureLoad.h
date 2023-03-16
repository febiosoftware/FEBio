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
#include "FEThermoFluid.h"
#include <FECore/FESurfaceConstraint.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEConstraintNormalFlow class implements a fluid surface with zero
//! tangential velocity as a linear constraint.

class FEBIOFLUID_API FEThermoFluidPressureLoad : public FESurfaceConstraint
{
public:
    //! constructor
    FEThermoFluidPressureLoad(FEModel* pfem);
    
    //! destructor
    ~FEThermoFluidPressureLoad() {}
    
    //! calculate traction stiffness (there is none for this load)
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate load vector (there is none for this load)
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! set the dilatation
    void Update() override;
    
    //! initialize
    bool Init() override;
    
    // allocate equations
    int InitEquations(int neq) override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;

    // return the surface
    FESurface* GetSurface() override { return &m_surf; }

protected:
    void UnpackLM(vector<int>& lm, int n);
    
    // Build the matrix profile
    void BuildMatrixProfile(FEGlobalMatrix& M) override;
    
    void Update(const std::vector<double>& Ui, const std::vector<double>& ui) override;
    void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) override;
    
    void PrepStep() override;
    
protected:
    int             m_dofT;
    int             m_dofEF;
    vector<int>     m_EQ;
    vector<double>  m_Lm, m_Lmp;
    FESurface       m_surf;
    
    FEThermoFluid*    m_pfluid; //!< pointer to thermo-fluid material

public:
    FEParamDouble   m_p;        // prescribed pressure
    
    DECLARE_FECORE_CLASS();
};
