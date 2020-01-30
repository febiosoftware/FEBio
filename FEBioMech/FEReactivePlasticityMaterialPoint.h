/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FECore/FEMaterialPoint.h"
#include "FEReactivePlasticity.h"
#include <vector>

class FEReactivePlasticity;

//-----------------------------------------------------------------------------
// Define a material point that stores the plasticity variables.
class FEReactivePlasticityMaterialPoint : public FEMaterialPoint
{
public:
    //! constructor
    FEReactivePlasticityMaterialPoint(FEMaterialPoint *pt, FEReactivePlasticity* pmat) : FEMaterialPoint(pt) { m_pMat = pmat; }
    
    FEMaterialPoint* Copy();
    
    //! Initialize material point data
    void Init();
    
    //! Update material point data
    void Update(const FETimeInfo& timeInfo);
    
    //! Serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! Evaluate net mass fraction of yielded bonds
    double YieldedBonds();
    
public:
    vector<mat3d>           m_Fusi;     //!< inverse of plastic deformation gradient at previous yield
    vector<mat3d>           m_Fvsi;     //!< trial value of plastic deformation gradient at current yield
    vector<double>          m_w;        //!< mass fraction of yielded bonds
    vector<double>          m_Kv;       //!< value of yield measure at current yield
    vector<double>          m_Ku;       //!< value of yield measure at previous yield
    vector<double>          m_gp;       //!< current value of octahedral plastic shear strain
    mat3d                   m_Fp;       //!< deformation gradient at previous time
    FEReactivePlasticity*   m_pMat;     //!< parent material
};
