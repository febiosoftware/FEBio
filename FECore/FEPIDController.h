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
#include "FELoadController.h"
#include "MathObject.h"

class FEPIDController : public FELoadController
{
public:
	FEPIDController(FEModel* fem);

	bool Init() override;

	double GetParameterValue() const { return m_paramVal; }
	double GetError() const { return m_error; }

	void Serialize(DumpStream& ar);

protected:
	double GetValue(double time) override;

private:
	std::string		m_paramName;	// the parameter to target
	double			m_trg;			// the target value

	double	m_Kp;
	double	m_Kd;
	double	m_Ki;

	FEParamValue	m_param;		//!< parameter that is tracked
	double			m_paramVal;		//!< current parameter value
	double			m_error;		//!< current error
	double			m_prev;			//!< error at previous time step
	double			m_prevTime;		//!< previous time step
	double			m_int;			//!< integral value

	DECLARE_FECORE_CLASS();
};
