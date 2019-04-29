/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for bond relaxation of reactive viscoelastic materials.
//! These materials need to define a relaxation function.
//!
class FEBondRelaxation : public FEMaterial
{
public:
	FEBondRelaxation(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEBondRelaxation() {}
    
	//! relaxation
	virtual double Relaxation(FEMaterialPoint& pt, const double t, const mat3ds D) = 0;
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
	double	m_tau;      //!< relaxation time
    
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
    double	m_tau0;     //!< relaxation time
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_alpha;    //!< exponent of 2nd term for tau
    
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
    double	m_tau1;     //!< lower relaxation time
    double  m_tau2;     //!< upper relaxation time
    
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
    double	m_tau;      //!< relaxation time
    double  m_beta;     //!< exponent
    
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
    double	m_tau0;      //!< relaxation time
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_beta0;     //!< exponent
    double  m_beta1;    //!< coefficient of 2nd for beta
    double  m_alpha;    //!< exponent of 2nd term for tau and beta
    
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
    double	m_tau;      //!< relaxation time
    double  m_beta;     //!< exponent
    
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
    double	m_tau0;      //!< relaxation time at zero strain
    double  m_beta0;     //!< exponent of relaxation power law
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_beta1;    //!< coefficient of 2nd for beta
    double  m_alpha;    //!< exponent of 2nd term
    
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
    double	m_tau0;		//!< characteristic time constant
    double  m_lam;      //!< time constant
    double  m_n;        //!< power-law index
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
