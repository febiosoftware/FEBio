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
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
class FECORE_API FESurfacePairConstraintNL : public FENLConstraint
{
public:
    //! constructor
    FESurfacePairConstraintNL(FEModel* pfem);
    
public:
    //! return the primary surface
    virtual FESurface* GetPrimarySurface() = 0;
    
    //! return the secondary surface
    virtual FESurface* GetSecondarySurface () = 0;
    
    //! temporary construct to determine if contact interface uses nodal integration rule (or facet)
    virtual bool UseNodalIntegration() = 0;
    
    //! create a copy of this interface
    virtual void CopyFrom(FESurfacePairConstraintNL* pci) {}
    
    virtual int InitEquations(int neq);
    
    void Update(vector<double>& ui);
    
public:
    
    using FENLConstraint::Update;
};
