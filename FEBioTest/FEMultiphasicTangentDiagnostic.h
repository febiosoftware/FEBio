//
//  FEMultiphasicTangentDiagnostic.h
//  FEBio2
//
//  Created by Gerard Ateshian on 8/21/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBio2__FEMultiphasicTangentDiagnostic__
#define __FEBio2__FEMultiphasicTangentDiagnostic__

#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FEMultiphasicScenario : public FEDiagnosticScenario
{
public:
	FEMultiphasicScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }

public:
	double	m_dt;
};

//-----------------------------------------------------------------------------
class FEMultiphasicTangentUniaxial : public FEMultiphasicScenario
{
public:
	FEMultiphasicTangentUniaxial(FEDiagnostic* pdia);

	bool Init() override;

private:
    double		m_strain;
    double      m_pressure;
    double      m_concentration;

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! The FEMultiphasicTangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEMultiphasicTangentDiagnostic : public FEDiagnostic
{
public:
    FEMultiphasicTangentDiagnostic(FEModel& fem);
    virtual ~FEMultiphasicTangentDiagnostic(){}
    
    bool Init();
    
    bool Run();

	FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void deriv_residual(matrix& ke);
    
    void print_matrix(matrix& m);
    
public:
    FEMultiphasicScenario*	m_pscn;
};

#endif /* defined(__FEBio2__FEMultiphasicTangentDiagnostic__) */
