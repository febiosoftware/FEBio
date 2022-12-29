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
#include <FECore/FEMaterialPoint.h>
#include "FEReactivePlasticity.h"
#include <vector>

//-----------------------------------------------------------------------------
// Define a material point that stores the plasticity variables.
class FEReactivePlasticityMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEReactivePlasticityMaterialPoint(FEMaterialPointData*pt, FEElasticMaterial* pmat) : FEMaterialPointData(pt) { m_pMat = pmat; }

	FEMaterialPointData* Copy();
    
    //! Initialize material point data
    void Init();
    
    //! Update material point data
    void Update(const FETimeInfo& timeInfo);
    
    //! Serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! Evaluate net mass fraction of yielded bonds
    double YieldedBonds() const;
    
public:
    vector<mat3d>           m_Fusi;     //!< inverse of plastic deformation gradient at previous yield
    vector<mat3d>           m_Fvsi;     //!< trial value of plastic deformation gradient at current yield
    vector<double>          m_w;        //!< mass fraction of yielded bonds
    vector<double>          m_Kv;       //!< value of yield measure at current yield
    vector<double>          m_Ku;       //!< value of yield measure at previous yield
    vector<double>          m_gp;       //!< current value of octahedral plastic shear strain
    vector<double>          m_gpp;      //!< previous value of octahedral plastic shear strain
    vector<double>          m_gc;       //!< cumulative value of octahedral plastic shear strain
    vector<bool>            m_byld;     //!< flag on which bonds have already yielded at start of current time
    mat3d                   m_Fp;       //!< deformation gradient at previous time
    double                  m_Rhat;     //!< reactive heat supply density
    FEElasticMaterial*      m_pMat;     //!< parent material
};
