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
#include <FECore/FEPrescribedBC.h>
#include <FECore/tens3d.h>

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation(FEModel* pfem);

public:
	void AddNode(int n);
	void AddNodes(const FENodeSet& set) override;

	void Activate() override;

	void Deactivate() override;

	void PrepStep(std::vector<double>& ui, bool brel) override;

	void Update() override;

	int Items() const { return (int) m_node.size(); }

	int NodeID(int index) const { return m_node[index]; }

	void SetDeformationGradient(const mat3d& F);

	void CopyFrom(FEPrescribedBC* pbc) override;

protected:
	vec3d NodeValue(const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	vector<int>	m_node;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation2O : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation2O(FEModel* pfem);

public:
	void SetReferenceNode(int n);
	void AddNode(int n);
	void AddNodes(const FENodeSet& set) override;

	bool Init() override;

	void Activate() override;

	void Deactivate() override;

	void PrepStep(std::vector<double>& ui, bool brel) override;

	void Update() override;

	int Items() const { return (int)m_node.size(); }

	int NodeID(int index) const { return m_node[index]; }

	void SetDeformationGradient(const mat3d& F);
	void SetDeformationHessian(const tens3drs& G);

	void CopyFrom(FEPrescribedBC* pbc) override;

public:
	void SetScale(double s, int lc = -1);

protected:
	vec3d NodeValue(const vec3d& X1, const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	tens3drs m_G;
	vector<int>	m_node;
	int	m_refNode;

	DECLARE_FECORE_CLASS();
};
