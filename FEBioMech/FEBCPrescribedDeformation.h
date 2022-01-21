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
#include <FECore/FEPrescribedBC.h>
#include <FECore/tens3d.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEBCPrescribedDeformation : public FEPrescribedNodeSet
{
public:
	FEBCPrescribedDeformation(FEModel* pfem);

	void SetDeformationGradient(const mat3d& F);

	void CopyFrom(FEBoundaryCondition* pbc) override;

protected:
	void GetNodalValues(int nodelid, std::vector<double>& val) override;

protected:
	double	m_scale;
	mat3d	m_F;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEBCPrescribedDeformation2O : public FEPrescribedNodeSet
{
public:
	FEBCPrescribedDeformation2O(FEModel* pfem);

	void SetScale(double s, int lc = -1);

	void SetReferenceNode(int n);

	bool Init() override;

	void SetDeformationGradient(const mat3d& F);
	void SetDeformationHessian(const tens3drs& G);

	void CopyFrom(FEBoundaryCondition* pbc) override;

protected:
	void GetNodalValues(int nodelist, std::vector<double>& val) override;

protected:
	double	m_scale;
	mat3d	m_F;
	tens3drs m_G;
	int	m_refNode;

	DECLARE_FECORE_CLASS();
};
