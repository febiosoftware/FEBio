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
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEIdealGasPressure : public FESurfaceLoad
{
public:
	//! constructor
	FEIdealGasPressure(FEModel* pfem);

	//! initialization
	bool Init() override;

	void Activate() override;

	void Update() override;

public:
	double GetCurrentPressure() const { return m_currentPressure; }

public:
	//! calculate residual
	void LoadVector(FEGlobalVector& R) override;

	//! calculate stiffness
	void StiffnessMatrix(FELinearSystem& LS) override;

protected:
	double	m_initialPressure;	//!< initial pressure value
	bool	m_bsymm;			//!< use symmetric formulation
	bool	m_bshellb;			//!< flag for prescribing pressure on shell bottom

private:
	double	m_initialVolume;
	double	m_currentVolume;
	double	m_currentPressure;

	DECLARE_FECORE_CLASS();
};
