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
    
	void Init();
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation triggered by any strain

class FEBondRelaxationExponential : public FEBondRelaxation
{
public:
	enum { MAX_TERMS = 6 };
    
public:
    //! constructor
    FEBondRelaxationExponential(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
    int     m_nt;           //!< number of non-zero relaxation times
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class implements exponential relaxation triggered by distortional strain

class FEBondRelaxationExpDistortion : public FEBondRelaxation
{
public:
    enum { MAX_TERMS = 6 };
    
public:
    //! constructor
    FEBondRelaxationExpDistortion(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
    double	m_t[MAX_TERMS];	//!< relaxation times
    int     m_nt;           //!< number of non-zero relaxation times
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class implements Fung's relaxation triggered by any strain

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
// This class implements Fung's relaxation triggered by distortional strain

class FEBondRelaxationFungDistortion : public FEBondRelaxation
{
public:
    //! constructor
    FEBondRelaxationFungDistortion(FEModel* pfem);
    
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
// This class implements Park's relaxation triggered by any strain

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
// This class implements Park's relaxation triggered by distortional strain

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
    double	m_tau;      //!< relaxation time
    double  m_beta;     //!< exponent
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEBondRelaxation__) */
