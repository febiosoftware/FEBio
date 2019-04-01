/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! This boundary condition sustains a fluid flux on a surface
//!
class FEBIOMIX_API FEFluidFlux : public FESurfaceLoad
{
public:
	//! constructor
	FEFluidFlux(FEModel* pfem);

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetMixture(bool bmix) { m_bmixture = bmix; }

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	//! calculate flux stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;

	//! unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, double dt, bool mixture);

	//! calculate stiffness for an element, for steady-state analysis
	void FluxStiffnessSS(FESurfaceElement& el, matrix& ke, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux, for steady-state analysis
	bool FlowRateSS(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower), for steady-state analysis
	bool LinearFlowRateSS(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture);

protected:
	FEParamDouble	m_flux;			//!< fluid flux
	bool	m_bmixture;		//!< mixture velocity or relative fluid flux
    bool    m_bshellb;      //!< flag for prescribing flux on shell bottom
	bool	m_blinear;		//!< type (linear or nonlinear)

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofP;
	int	m_dofVX;
	int	m_dofVY;
	int	m_dofVZ;
    int	m_dofSX;
    int	m_dofSY;
    int	m_dofSZ;
    int	m_dofQ;
    int	m_dofSVX;
    int	m_dofSVY;
    int	m_dofSVZ;

	DECLARE_FECORE_CLASS();
};
