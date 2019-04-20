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
#include <FECore/FESurfaceMap.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! This boundary condition applies a poro-elastic normal traction on a surface
//!
class FEBIOMIX_API FEPoroNormalTraction : public FESurfaceLoad
{
public:
	//! constructor
	FEPoroNormalTraction(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetEffective(bool beff) { m_beffective = beff; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate residual
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

protected:
	//! calculate stiffness for an element
	void TractionStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn, bool effective, bool bsymm);

	//! Calculates external pressure forces
	bool TractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearTractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	double	m_traction;		//!< traction value
	bool	m_blinear;		//!< linear or not (true is non-follower, false is follower)
    bool    m_bshellb;      //!< flag for prescribing traction on shell bottom
	bool	m_beffective;	//!< effective or total normal traction

	// pressure boundary data
	FESurfaceMap	m_PC;		//!< pressure boundary cards

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofP;
    int	m_dofSX;
    int	m_dofSY;
    int	m_dofSZ;
    int	m_dofQ;

	DECLARE_FECORE_CLASS();
};
