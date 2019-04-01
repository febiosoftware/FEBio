/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEReactiveViscoelastic.h"
#include "FEUncoupledReactiveViscoelastic.h"
#include <deque>

class FEReactiveViscoelasticMaterial;
class FEUncoupledReactiveViscoelasticMaterial;

//-----------------------------------------------------------------------------
//! Material point data for reactive viscoelastic materials
class FEReactiveVEMaterialPoint : public FEMaterialPoint
{
public:
    //! olverloaded constructors
    FEReactiveVEMaterialPoint(FEMaterialPoint *pt, FEReactiveViscoelasticMaterial *pe) : FEMaterialPoint(pt) { m_pRve = pe; m_pRuc = 0; }
    FEReactiveVEMaterialPoint(FEMaterialPoint *pt, FEUncoupledReactiveViscoelasticMaterial *pe) : FEMaterialPoint(pt) { m_pRve = 0; m_pRuc = pe; }
    
    //! copy material point data
    FEMaterialPoint* Copy();
    
    //! Initialize material point data
    void Init();

    //! Update material point data
    void Update(const FETimeInfo& timeInfo);
    
    //! Serialize data to archive
    void Serialize(DumpStream& ar);
    
public:
    // multigenerational material data
    deque <mat3d>  m_Fi;	//!< inverse of relative deformation gradient
    deque <double> m_Ji;	//!< determinant of Fi (store for efficiency)
    deque <double> m_v;     //!< time when generation starts breaking
    deque <double> m_w;     //!< mass fraction when generation starts breaking
    FEReactiveViscoelasticMaterial*  m_pRve; //!< pointer to parent material
    FEUncoupledReactiveViscoelasticMaterial*  m_pRuc; //!< pointer to parent material
};
