/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "FEDamageMaterialPoint.h"
#include <deque>
#include <FECore/FEMaterialPoint.h>

//-----------------------------------------------------------------------------
// structure for fatigue bonds
class FatigueBond
{
public:
    // constructor
    FatigueBond();
    
    // copy constructor
    FatigueBond(const FatigueBond& fb);
    FatigueBond(FatigueBond& fb);
    // update
    void Update();
    
public:
    double  m_wft;      //!< fatigued bond fraction at current time
    double  m_wfp;      //!< fatigued bond fraction at previous time
    double  m_Xfmax;    //!< max damage criterion for fatigued bonds
    double  m_Xftrl;    //!< trial value of Xfmax
    double  m_Fft;      //!< fatigue bond damage CDF at current time
    double  m_Ffp;      //!< fatigue bond damage CDF at previous time
    double  m_time;     //!< fatigue bond generation time
    bool    m_erase;    //!< flag for erasing a generation
};

//-----------------------------------------------------------------------------
// Define a material point that stores the fatigue and damage variables.
class FEReactiveFatigueMaterialPoint : public FEDamageMaterialPoint
{
public:
    // default constructor
    FEReactiveFatigueMaterialPoint(FEMaterialPointData*pt);
    
	// copy constructors
    FEReactiveFatigueMaterialPoint(const FEReactiveFatigueMaterialPoint& rfmp);
    FEReactiveFatigueMaterialPoint(FEReactiveFatigueMaterialPoint& rfmp);
    
	FEMaterialPointData* Copy();
    
    void Init();
    void Update(const FETimeInfo& timeInfo);
    
    void Serialize(DumpStream& ar);
    
public:
    double      m_wit;          //!< intact bond mass fraction at current time
    double      m_wip;          //!< intact bond mass fraction at previous time
    
    double      m_Ximax;        //!< max damage criterion for intact bonds
    double      m_Xitrl;        //!< trial value of Ximax
    double      m_Xip;          //!< Xi at previous time
    
    double      m_aXit;         //!< rate of change of Xi at current time
    double      m_aXip;         //!< rate of change of Xi at previous time
    
    double      m_Fit;          //!< intact bond damage CDF at current time
    double      m_Fip;          //!< intact bond damage CDF at previous time
    
    double      m_wbt;          //!< broken (damaged) bond fraction at current time
    double      m_wbp;          //!< broken (damaged) bond fraction at previous time
    
    double      m_wft;          //!< fatigue bond fraction at current time
    
    
    std::deque <FatigueBond> m_fb;   //!< generations of fatigued bonds
};

