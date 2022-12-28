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
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! 2D transversely isotropic Mooney-Rivlin

//! This class describes a transversely isotropic matrix where the base material
//! is Mooney-Rivlin. The difference between this material and the FETransIsoMooneyRivlin
//! material is that in this material the fibers lie in the plane that is perpendicular
//! to the transverse axis. 

class FE2DTransIsoMooneyRivlin : public FEUncoupledMaterial
{
	enum { NSTEPS = 12 };	// nr of integration steps

public:
	// material parameters
	double	m_c1;	//!< Mooney-Rivlin parameter c1
	double	m_c2;	//!< Mooney-Rivlin parameter c2

	// fiber parameters
	double	m_c3;
	double	m_c4;
	double	m_c5;
	double	m_lam1;
	double	m_w[2];
	double	m_epsf;
	FEVec3dValuator*	m_fiber;

	//--- active contraction stuff ---
	double	m_a[2];
	double	m_ac;
	// -------------------------------

public:
	//! constructor
	FE2DTransIsoMooneyRivlin(FEModel* pfem);
	
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	// declare parameter list
	DECLARE_FECORE_CLASS();

protected:
	static double	m_cth[NSTEPS];
	static double	m_sth[NSTEPS];
};
