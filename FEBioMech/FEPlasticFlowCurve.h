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
#include <FECore/FEMaterial.h>
#include <FECore/FEPointFunction.h>

//-----------------------------------------------------------------------------
// Virtual base class for plastic flow curve

class FEPlasticFlowCurve : public FEMaterialProperty
{
public:
    FEPlasticFlowCurve(FEModel* pfem) : FEMaterialProperty(pfem) { m_binit = false; }

    //! return
    vector<double> BondYieldMeasures() { return Ky; }
    vector<double> BondMassFractions() { return w; }
    size_t BondFamilies() { return Ky.size(); }
    
protected:
    vector<double> Ky;
    vector<double> w;
    bool    m_binit; //!< flag indicating whether the material has been initialized

    FECORE_BASE_CLASS(FEPlasticFlowCurve);
};

//-----------------------------------------------------------------------------
// Flow curve from reactive plasticity paper

class FEPlasticFlowCurvePaper : public FEPlasticFlowCurve
{
public:
    FEPlasticFlowCurvePaper(FEModel* pfem);
    
    bool Init() override;
    
public:
    double      m_wmin;     // initial fraction of yielding bonds
    double      m_wmax;     // final fraction of yielding bonds
    double      m_we;       // fraction of unyielding bonds
    double      m_Ymin;     // initial yield measure
    double      m_Ymax;     // yield measure when all bonds have yielded
    int         m_n;        // number of yield levels
    double      m_bias;     // biasing factor for intervals in yield measures and bond fractions

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// User-defined flow curve

class FEPlasticFlowCurveUser : public FEPlasticFlowCurve
{
public:
    FEPlasticFlowCurveUser(FEModel* pfem);
    
    bool Init() override;
    
public:
    FEPointFunction*    m_Y;    //!< true stress-true strain flow curve

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Math equation for flow curve

class FEPlasticFlowCurveMath : public FEPlasticFlowCurve
{
public:
    FEPlasticFlowCurveMath(FEModel* pfem);
    
    bool Init() override;
    
public:
    int                 m_n;    //!< number of yield levels
    double              m_emin; //!< min true strain
    double              m_emax; //!< max true strain
    std::string         m_Ymath;    //!< true stress-true strain flow curve
    
    DECLARE_FECORE_CLASS();
};

