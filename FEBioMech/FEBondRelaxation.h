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
//! Base class for bond relaxation of reactive viscoelastic materials.
//! These materials need to define a relaxation function.
//!
class FEBIOMECH_API FEBondRelaxation : public FEMaterialProperty
{
public:
	FEBondRelaxation(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~FEBondRelaxation() {}
    
	//! relaxation
	virtual double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) = 0;

    FECORE_BASE_CLASS(FEBondRelaxation)
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation with constant relaxation time

class FEBondRelaxationExponential : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationExponential(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
	FEParamDouble   m_tau;      //!< relaxation time
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationExpDistortion : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationExpDistortion(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble  m_tau0;     //!< relaxation time
    FEParamDouble  m_tau1;     //!< relaxation time coeff. of 2nd term
    FEParamDouble  m_alpha;    //!< exponent of 2nd term for tau

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationExpDistUser : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationExpDistUser(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! performs initialization
    bool Init() override;

public:
    FEFunction1D*   m_tau;      //!< relaxation time
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements Fung's relaxation

class FEBondRelaxationFung : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationFung(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! data initialization and checking
    bool Validate() override;
    
public:
    FEParamDouble   m_tau1;     //!< lower relaxation time
    FEParamDouble   m_tau2;     //!< upper relaxation time
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements Park's relaxation with constant relaxation time

class FEBondRelaxationPark : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPark(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau;      //!< relaxation time
    FEParamDouble   m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements Park's relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationParkDistortion : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationParkDistortion(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau0;     //!< relaxation time
    FEParamDouble   m_tau1;     //!< relaxation time coeff. of 2nd term
    FEParamDouble   m_beta0;    //!< exponent
    FEParamDouble   m_beta1;    //!< coefficient of 2nd for beta
    FEParamDouble   m_alpha;    //!< exponent of 2nd term for tau and beta

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements Park's relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationParkDistUser : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationParkDistUser(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFunction1D*   m_tau;      //!< relaxation time
    FEFunction1D*   m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a power-law relaxation with constant relaxation time

class FEBondRelaxationPower : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPower(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau;      //!< relaxation time
    FEParamDouble   m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a power-law relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationPowerDistortion : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPowerDistortion(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau0;     //!< relaxation time at zero strain
    FEParamDouble   m_beta0;    //!< exponent of relaxation power law
    FEParamDouble   m_tau1;     //!< relaxation time coeff. of 2nd term
    FEParamDouble   m_beta1;    //!< coefficient of 2nd for beta
    FEParamDouble   m_alpha;    //!< exponent of 2nd term

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a power-law relaxation with a relaxation
// time that is a function of the distortional strain

class FEBondRelaxationPowerDistUser : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPowerDistUser(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFunction1D*   m_tau;      //!< relaxation time
    FEFunction1D*   m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a relaxation that produces Carreau fluid response
// under stead-state conditions

class FEBondRelaxationCarreau : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationCarreau(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau0;     //!< characteristic time constant
    FEParamDouble   m_lam;      //!< time constant
    FEParamDouble   m_n;        //!< power-law index
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements Prony series exponential relaxation with constant relaxation time

class FEBondRelaxationProny : public FEBondRelaxation
{
public:
    enum { MAX_TERMS = 6 };

public:
    //! constructor
    FEBondRelaxationProny(FEModel* pfem);
    
    //! data initialization and checking
    bool Validate() override;
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    double  m_g[MAX_TERMS];     //!< viscoelastic coefficients
    double  m_t[MAX_TERMS];     //!< relaxation times
    double  m_sg;               //!< sum of viscoelastic coefficients

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a Malkin relaxation with constant relaxation time

class FEBondRelaxationMalkin : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationMalkin(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble  m_tau1;     //!< lower relaxation time
    FEParamDouble  m_tau2;     //!< upper relaxation time
    FEParamDouble  m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a Malkin relaxation with adjustable relaxation time

class FEBondRelaxationMalkinDist : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationMalkinDist(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble  m_t1c0;     //!< constant coefficient of lower relaxation time
    FEParamDouble  m_t1c1;     //!< coefficient of exponential for lower relaxation time
    FEParamDouble  m_t1s0;     //!< time constant of exponential for lower relaxation time
    FEParamDouble  m_t2c0;     //!< constant coefficient of upper relaxation time
    FEParamDouble  m_t2c1;     //!< coefficient of exponential for upper relaxation time
    FEParamDouble  m_t2s0;     //!< time constant of exponential for upper relaxation time
    FEParamDouble  m_beta;     //!< exponent

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a Malkin relaxation with user-adjustable relaxation time

class FEBondRelaxationMalkinDistUser : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationMalkinDistUser(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFunction1D*   m_tau1;     //!< lower relaxation time
    FEFunction1D*   m_tau2;     //!< upper relaxation time
    FEFunction1D*   m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class implements a continuous spectrum exponential relaxation with constant relaxation time

class FEBondRelaxationCSexp : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationCSexp(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
public:
    FEParamDouble   m_tau;      //!< relaxation time
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
// This class implements a continuous spectrum exponential relaxation with adjustable relaxation time

class FEBondRelaxationCSexpDistUser : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationCSexpDistUser(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) override;
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFunction1D*   m_tau;      //!< relaxation time
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

