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

#include "FEDiscreteElasticMaterial.h"
#include <FECore/FEFunction1D.h>

class FEDiscreteContractileMaterial : public FEDiscreteElasticMaterial
{
public:
	FEDiscreteContractileMaterial(FEModel* fem);

	// evaluate the force at a discrete element
	vec3d Force(FEDiscreteMaterialPoint& mp) override;

	// evaluate the stiffness at a discrete element (= dF / dr)
	mat3d Stiffness(FEDiscreteMaterialPoint& mp) override;

private:
	double force(FEDiscreteMaterialPoint& mp);
	double force_deriv_L(FEDiscreteMaterialPoint& mp);
	double force_deriv_V(FEDiscreteMaterialPoint& mp);
	double passive_force(double L, double V);
	double active_force(double L, double V);
	double passive_force_deriv_L(double L, double V);
	double passive_force_deriv_V(double L, double V);
	double active_force_deriv_L(double L, double V);
	double active_force_deriv_V(double L, double V);

private:
	double			m_Vmax;		// maximum shortening velocity
	double			m_ac;		// activation level
	double			m_Fmax;		// max force
	double			m_Ksh;		// shape parameter that determines the rise of exponential in passive element
	double			m_Lmax;		// relative length at which Fmax occurs
	double			m_L0;		// initial reference length

	FEFunction1D*	m_Sv;		// max velocity scale
	FEFunction1D*	m_Ftl;		// normalized tension-length curve for active element
	FEFunction1D*	m_Ftv;		// normalized tension-velocity curve for active element

	DECLARE_FECORE_CLASS();
};
