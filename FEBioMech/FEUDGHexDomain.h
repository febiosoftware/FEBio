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
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! domain class for uniform-deformation-gradient hex elements (UDG)
class FEBIOMECH_API FEUDGHexDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEUDGHexDomain(FEModel* pfem);

	//! Set the hourglass parameter
	void SetHourGlassParameter(double hg);

	//! allocate elements
	bool Create(int nelems, FE_Element_Spec spec) override;

public:
	//! calculates the residual
	void InternalForces(FEGlobalVector& R) override;

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;

protected: // element residual contributions
	//! Calculates the internal stress vector for enhanced strain hex elements
	void UDGInternalForces(FESolidElement& el, vector<double>& fe);

	//! calculates hourglass forces for the UDG element
	void UDGHourglassForces(FESolidElement& el, vector<double>& fe);

protected: // element stiffness contributions
	//! hourglass stiffness for UDG hex elements
	void UDGHourglassStiffness(FEModel& fem, FESolidElement& el, matrix& ke);

	//! geometrical stiffness for UDG hex elements
	void UDGGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness for UDG hex elements
	void UDGMaterialStiffness(FESolidElement& el, matrix& ke);

protected:
	void AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int state = 0);
	void AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8]);
	double HexVolume(FESolidElement& el, int state = 0);

public:
	double	m_hg;	//!< hourglass parameter

	DECLARE_FECORE_CLASS();
};
