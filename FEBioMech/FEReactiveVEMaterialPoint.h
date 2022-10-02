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
#include "FECore/FEMaterialPoint.h"
#include "FEReactiveViscoelastic.h"
#include "FEUncoupledReactiveViscoelastic.h"
#include <deque>

class FEReactiveViscoelasticMaterial;
class FEUncoupledReactiveViscoelasticMaterial;

//-----------------------------------------------------------------------------
//! Material point data array for reactive viscoelastic materials
//!
class FEReactiveViscoelasticMaterialPoint : public FEMaterialPointArray
{
public:
    //! constructor
    FEReactiveViscoelasticMaterialPoint();
    
    //! Copy material point data
    FEMaterialPointData* Copy();
};

//-----------------------------------------------------------------------------
//! Material point data for reactive viscoelastic materials
class FEReactiveVEMaterialPoint : public FEMaterialPointData
{
public:
    //! olverloaded constructors
    FEReactiveVEMaterialPoint(FEMaterialPointData*pt) : FEMaterialPointData(pt) {}
    
    //! copy material point data
	FEMaterialPointData* Copy() override;
    
    //! Initialize material point data
    void Init() override;
    
    //! Update material point data
    void Update(const FETimeInfo& timeInfo) override;

    //! Serialize data to archive
    void Serialize(DumpStream& ar) override;
    
public:
    // multigenerational material data
    deque <mat3ds> m_Uv;	//!< right stretch tensor at tv (when generation u starts breaking)
    deque <double> m_Jv;	//!< determinant of Uv (store for efficiency)
    deque <double> m_v;     //!< time tv when generation starts breaking
    deque <double> m_f;     //!< mass fraction when generation starts breaking
    
public:
    // weak bond recruitment parameters
    double m_Et;            //!< trial strain value at time t
    double m_Em;            //!< max strain value up to time t
    deque <double> m_wv;    //!< total mass fraction of weak bonds
};
