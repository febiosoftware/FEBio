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
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! This boundary condition applies a poro-elastic normal traction on a surface
//!
class FEBIOMIX_API FEPoroNormalTraction : public FESurfaceLoad
{
public:
	//! constructor
	FEPoroNormalTraction(FEModel* pfem);

	bool Init() override;

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetEffective(bool beff) { m_beffective = beff; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! calculate residual
	void LoadVector(FEGlobalVector& R) override;

private:
	double Traction(FESurfaceMaterialPoint& mp);

protected:
	FEParamDouble	m_traction;		//!< traction value
	bool	m_blinear;		//!< linear or not (true is non-follower, false is follower)
    bool    m_bshellb;      //!< flag for prescribing traction on shell bottom
	bool	m_beffective;	//!< effective or total normal traction

	DECLARE_FECORE_CLASS();
};
