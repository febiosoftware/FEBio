//
//  FEBondRelaxation.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEBondRelaxation__
#define __FEBioMech__FEBondRelaxation__

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
	virtual double Relaxation(FEMaterialPoint& pt, const double t) = 0;
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation with constant relaxation time

class FEBondRelaxationExponential : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationExponential(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
	double	m_tau;      //!< relaxation time
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
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
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau0;     //!< relaxation time
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_alpha;    //!< exponent of 2nd term for tau
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class implements Fung's relaxation

class FEBondRelaxationFung : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationFung(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau1;     //!< lower relaxation time
    double  m_tau2;     //!< upper relaxation time
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class implements Park's relaxation with constant relaxation time

class FEBondRelaxationPark : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPark(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau;      //!< relaxation time
    double  m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
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
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau0;      //!< relaxation time
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_beta0;     //!< exponent
    double  m_beta1;    //!< coefficient of 2nd for beta
    double  m_alpha;    //!< exponent of 2nd term for tau and beta
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class implements a power-law relaxation with constant relaxation time

class FEBondRelaxationPower : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationPower(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau;      //!< relaxation time
    double  m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
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
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_tau0;      //!< relaxation time at zero strain
    double  m_beta0;     //!< exponent of relaxation power law
    double	m_tau1;     //!< relaxation time coeff. of 2nd term
    double  m_beta1;    //!< coefficient of 2nd for beta
    double  m_alpha;    //!< exponent of 2nd term
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEBondRelaxation__) */
