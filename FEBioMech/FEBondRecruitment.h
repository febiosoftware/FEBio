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
#include <FECore/FEFunction1D.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
// Virtual base class for damage bond recruitment functions

class FEBIOMECH_API FEBondRecruitment : public FEMaterialProperty
{
public:
    FEBondRecruitment(FEModel* pfem) : FEMaterialProperty(pfem) {}
    
    //! bond recruitment function
    virtual double brf(FEMaterialPoint& mp, const double X) = 0;
    
public:
    
    FECORE_BASE_CLASS(FEBondRecruitment)
};

//-----------------------------------------------------------------------------
// User-specified load curve for damage bond recruitment function

class FEBondRecruitmentUser : public FEBondRecruitment
{
public:
    FEBondRecruitmentUser(FEModel* pfem);
    ~FEBondRecruitmentUser() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEFunction1D*    m_brf;           //!< user-defined BRF
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Power law bond recruitment function

class FEBondRecruitmentPower : public FEBondRecruitment
{
public:
    FEBondRecruitmentPower(FEModel* pfem);
    ~FEBondRecruitmentPower() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble    m_alpha;            //!< power exponent alpha
    FEParamDouble    m_mu0;              //!< constant coeff
    FEParamDouble    m_mu1;              //!< coeff of power
    FEParamDouble    m_s;                //!< scale factor for argument

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Exponential law bond recruitment function

class FEBondRecruitmentExp : public FEBondRecruitment
{
public:
    FEBondRecruitmentExp(FEModel* pfem);
    ~FEBondRecruitmentExp() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble    m_alpha;            //!< power exponent alpha
    FEParamDouble    m_mu0;              //!< constant coeff
    FEParamDouble    m_mu1;              //!< coeff of power
    FEParamDouble    m_s;                //!< scale factor for argument

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Quadratic polynomial bond recruitment function

class FEBondRecruitmentPoly : public FEBondRecruitment
{
public:
    FEBondRecruitmentPoly(FEModel* pfem);
    ~FEBondRecruitmentPoly() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble    m_mu0;              //!< constant coeff
    FEParamDouble    m_mu1;              //!< coeff of linear term
    FEParamDouble    m_mu2;              //!< coeff of quadratic term

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Log-normal damage bond recruitment function

class FEBondRecruitmentLogNormal : public FEBondRecruitment
{
public:
    FEBondRecruitmentLogNormal(FEModel* pfem);
    ~FEBondRecruitmentLogNormal() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble   m_mu;       //!< mean on log scale
    FEParamDouble   m_sigma;    //!< standard deviation on log scale
    FEParamDouble   m_max;      //!< maximum increase in recruitment
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Weibull damage bond recruitment function

class FEBondRecruitmentWeibull : public FEBondRecruitment
{
public:
    FEBondRecruitmentWeibull(FEModel* pfem);
    ~FEBondRecruitmentWeibull() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble   m_alpha;    //!< exponent alpha
    FEParamDouble   m_mu;       //!< mean mu
    FEParamDouble   m_ploc;     //!< location parameter
    FEParamDouble   m_max;      //!< maximum increase in recruitment

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial bond recruitment function

class FEBondRecruitmentPQP : public FEBondRecruitment
{
public:
    FEBondRecruitmentPQP(FEModel* pfem);
    ~FEBondRecruitmentPQP() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
    bool Validate() override;
    
public:
    FEParamDouble   m_mumin;    //!< mu threshold
    FEParamDouble   m_mumax;    //!< mu cap
    FEParamDouble   m_max;      //!< maximum increase in recruitment

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Gamma damage bond recruitment function

class FEBondRecruitmentGamma : public FEBondRecruitment
{
public:
    FEBondRecruitmentGamma(FEModel* pfem);
    ~FEBondRecruitmentGamma() {}
    
    //! bond recruitment function
    double brf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble   m_alpha;    //!< exponent alpha
    FEParamDouble   m_mu;       //!< pdf expected mean mu
    FEParamDouble   m_max;      //!< maximum increase in recruitment

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

