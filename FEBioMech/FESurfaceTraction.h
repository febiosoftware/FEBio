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

class FESurfaceTraction : public FESurfaceLoad
{
public:
	// Constructor
	FESurfaceTraction(FEModel* fem);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! helper function for getting nodal coordinates
	void GetNodalCoordinates(FESurfaceElement& el, vec3d* re);

	//! helper function for getting nodal coordinates
	void GetReferenceNodalCoordinates(FESurfaceElement& el, vec3d* re);

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;

	//! calculate stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! Is the traction linear?
	bool IsLinear() const { return m_blinear; }

	//! Set the traction load linear or not
	void SetLinear(bool b) { m_blinear = b; }

protected:
	//! element contribution to residual
	virtual void ElementResidual(FESurfaceElement& el, std::vector<double>& fe);

	//! element contribution to stiffness
	virtual void ElementStiffness(FESurfaceElement& el, matrix& ke) {}

	//! traction at material point
	virtual vec3d Traction(const FESurfaceMaterialPoint& mp) { return vec3d(0, 0, 0);  }

protected:
	// degrees of freedom
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofSX;
	int	m_dofSY;
	int	m_dofSZ;

	bool	m_blinear;	//!< is the load linear (i.e. it will be calculated in the reference frame and assummed deformation independent)
	bool	m_bshellb;	//!< flag for prescribing pressure on shell bottom

	DECLARE_FECORE_CLASS();
};
