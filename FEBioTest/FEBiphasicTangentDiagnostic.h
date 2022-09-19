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



#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FEBiphasicScenario : public FEDiagnosticScenario
{
public:
	FEBiphasicScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }

public:
	double	m_dt;
};

//-----------------------------------------------------------------------------
class FEBiphasicTangentUniaxial : public FEBiphasicScenario
{
public:
	FEBiphasicTangentUniaxial(FEDiagnostic* pdia);

	bool Init() override;

private:
    double		m_strain;
    double      m_pressure;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! The FEBiphasicTangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEBiphasicTangentDiagnostic : public FEDiagnostic
{
public:
    FEBiphasicTangentDiagnostic(FEModel* fem);
    virtual ~FEBiphasicTangentDiagnostic(){}
    
    bool Init();
    
    bool Run();

	FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void deriv_residual(matrix& ke);
    
    void print_matrix(matrix& m);
    
public:
    FEBiphasicScenario*	m_pscn;
};
