/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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

//-----------------------------------------------------------------------------
// Define a material point that stores fiber CDF data and integrates it
class FEFiberCDFMaterialPoint : public FEMaterialPointData
{
public:
    FEFiberCDFMaterialPoint(FEMaterialPointData* pt);
    
    FEMaterialPointData* Copy() override;
    
    void Init() override;
    
    void Serialize(DumpStream& ar) override;
    
    void Update(const FETimeInfo& timeInfo) override;
    
public:
    void SetFiberStrain(const double In_1);
    
    void Integrate(const double cdf);
    
public:
    double  m_In_1_t, m_In_1_p; //!< In - 1 (square of stretch ratio along fiber - 1) at current and previous time
    double  m_sed_t, m_sed_p;   //!< normalized sed at current and previous time
    double  m_ds_t, m_ds_p;     //!< first derivative of normalized sed w.r.t. In at current and previous time
    double  m_d2s_t, m_d2s_p;   //!< second derivative of normalized sed w.r.t. In at current and previous time
};

