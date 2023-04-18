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

class FEBIOMECH_API FEDamageCDF : public FEMaterialProperty
{
public:
    FEDamageCDF(FEModel* pfem) : FEMaterialProperty(pfem) { m_Dmax = 1; }
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
    //! cumulative distribution function
    virtual double cdf(FEMaterialPoint& mp, const double X) = 0;
    
    //! probability density function
    virtual double pdf(FEMaterialPoint& mp, const double X) = 0;

public:
    double  m_Dmax;              //!< maximum allowable damage
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
    FECORE_BASE_CLASS(FEDamageCDF)
};

//-----------------------------------------------------------------------------
// Simo damage cumulative distribution function
// Simo, CMAME 60 (1987), 153-173

class FEDamageCDFSimo : public FEDamageCDF
{
public:
	FEDamageCDFSimo(FEModel* pfem);
	~FEDamageCDFSimo() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;

public:
    FEParamDouble	m_alpha;			//!< parameter alpha
    FEParamDouble	m_beta;             //!< parameter beta
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Log-normal damage cumulative distribution function

class FEDamageCDFLogNormal : public FEDamageCDF
{
public:
	FEDamageCDFLogNormal(FEModel* pfem);
	~FEDamageCDFLogNormal() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;

public:
    FEParamDouble	m_mu;               //!< mean on log scale
    FEParamDouble	m_sigma;            //!< standard deviation on log scale
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Weibull damage cumulative distribution function

class FEDamageCDFWeibull : public FEDamageCDF
{
public:
	FEDamageCDFWeibull(FEModel* pfem);
	~FEDamageCDFWeibull() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;

public:
    FEParamDouble	m_alpha;            //!< exponent alpha
    FEParamDouble	m_mu;               //!< mean mu
    FEParamDouble  m_ploc;             //!< location parameter
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Step cumulative distribution function (sudden fracture)

class FEDamageCDFStep : public FEDamageCDF
{
public:
	FEDamageCDFStep(FEModel* pfem);
	~FEDamageCDFStep() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;

public:
    FEParamDouble	m_mu;               //!< threshold mu
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial cumulative distribution function

class FEDamageCDFPQP : public FEDamageCDF
{
public:
	FEDamageCDFPQP(FEModel* pfem);
	~FEDamageCDFPQP() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;

	bool Validate() override;
    
public:
    FEParamDouble	m_mumin;            //!< mu threshold
    FEParamDouble	m_mumax;            //!< mu cap
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Gamma damage cumulative distribution function

class FEDamageCDFGamma : public FEDamageCDF
{
public:
    FEDamageCDFGamma(FEModel* pfem);
    ~FEDamageCDFGamma() {}
    
    //! cumulative distribution function
    double cdf(FEMaterialPoint& mp, const double X) override;
    
    //! probability density function
    double pdf(FEMaterialPoint& mp, const double X) override;
    
public:
    FEParamDouble    m_alpha;            //!< exponent alpha
    FEParamDouble    m_mu;               //!< pdf expected mean mu
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// User-specified load curve for damage cumulative distribution function

class FEDamageCDFUser : public FEDamageCDF
{
public:
    FEDamageCDFUser(FEModel* pfem);
    ~FEDamageCDFUser() {}
    
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
// Power law cumulative distribution function

class FEDamageCDFPower : public FEDamageCDF
{
public:
    FEDamageCDFPower(FEModel* pfem);
    ~FEDamageCDFPower() {}
    
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
// Exponential law cumulative distribution function

class FEDamageCDFExp : public FEDamageCDF
{
public:
    FEDamageCDFExp(FEModel* pfem);
    ~FEDamageCDFExp() {}
    
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
// Quadratic polynomial cumulative distribution function

class FEDamageCDFPoly2 : public FEDamageCDF
{
public:
    FEDamageCDFPoly2(FEModel* pfem);
    ~FEDamageCDFPoly2() {}
    
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

