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
#include <FECore/FEAugLagLinearConstraint.h>
#include <FECore/FESurface.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEConstraintFrictionlessWall class implements a frictionless fluid wall
//! as a linear constraint on the components of the fluid velocity.

class FEBIOFLUID_API FEConstraintFrictionlessWall : public FESurfaceConstraint
{
public:
    //! constructor
    FEConstraintFrictionlessWall(FEModel* pfem);
    
    //! destructor
    ~FEConstraintFrictionlessWall() {}
    
    //! Activation
    void Activate() override;

    //! initialization
    bool Init() override;
    
    // allocate equations
    int InitEquations(int neq) override;
    
protected:
    void Update(const std::vector<double>& Ui, const std::vector<double>& ui) override;
    void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) override;
    
    void PrepStep() override;
    
    //! Get the surface
    FESurface* GetSurface() override { return &m_surf; }
    
public:
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;

    //! add the linear constraint contributions to the residual
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

    //! add the linear constraint contributions to the stiffness matrix
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

    //! do the augmentation
    bool Augment(int naug, const FETimeInfo& tp) override;

    //! build connectivity for matrix profile
    void BuildMatrixProfile(FEGlobalMatrix& M) override;

protected:
    FESurface	m_surf;
    FELinearConstraintSet   m_lc;
    bool        m_binit;

    DECLARE_FECORE_CLASS();
};
