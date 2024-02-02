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
// Virtual base class for damage cumulative distribution functions

class FEBIOMECH_API FEBondRecruitment : public FEMaterialProperty
{
public:
    FEBondRecruitment(FEModel* pfem) : FEMaterialProperty(pfem) {}
    
    //! cumulative distribution function
    virtual double cdf(FEMaterialPoint& mp, const double X) = 0;
    
    //! probability density function
    virtual double pdf(FEMaterialPoint& mp, const double X) = 0;

public:
    
    FECORE_BASE_CLASS(FEBondRecruitment)
};

//-----------------------------------------------------------------------------
// User-specified load curve for damage cumulative distribution function

class FEBondRecruitmentUser : public FEBondRecruitment
{
public:
    FEBondRecruitmentUser(FEModel* pfem);
    ~FEBondRecruitmentUser() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEFunction1D*    m_cdf;           //!< user-defined CDF
    
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
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;
    
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
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;
    
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
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble    m_mu0;              //!< constant coeff
    FEParamDouble    m_mu1;              //!< coeff of linear term
    FEParamDouble    m_mu2;              //!< coeff of quadratic term

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

