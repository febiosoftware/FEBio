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
#include "FEPrescribedBC.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//! prescribed boundary condition data
//! \todo Should I make a derived class for the relative prescribed BC's?
class FECORE_API FEPrescribedDOF : public FEPrescribedBC
{
	struct ITEM
	{
		int		nid;	// nodal ID
		double	ref;	// reference value (for relative BC's)

		void Serialize(DumpStream& ar);
	};

public:
	FEPrescribedDOF(FEModel* pfem);
	FEPrescribedDOF(FEModel* pfem, const FEPrescribedDOF& bc);

	void AddNode(int node, double scale = 1.0);

	using FEPrescribedBC::AddNodes;
	void AddNodes(const FENodeSet& s, double scale);
	void AddNodes(const FENodeSet& s) override { AddNodes(s, 1.0); }

	int NodeID(int i) { return m_item[i].nid; }

	size_t Items() const { return m_item.size(); }

	const ITEM& GetItem(int i) const { return m_item[i]; }

public:
	void Serialize(DumpStream& ar) override;

	void Activate() override;

	void Deactivate() override;

	bool Init() override;

	void Update() override;

	void PrepStep(std::vector<double>& ui, bool brel = true) override;

	void CopyFrom(FEPrescribedBC* pbc) override;

public:
	FEPrescribedDOF& SetScale(double s, int lc = -1);
	FEPrescribedDOF& SetDOF(int dof) { m_dof = dof; return *this; }
	FEPrescribedDOF& SetRelativeFlag(bool br) { m_br = br; return *this; }

	int GetDOF() const { return m_dof; }

	double NodeValue(int n);

private:
	int			m_dof;		//!< dof
	FEParamDouble	m_scale;	//!< overall scale factor
	bool		m_br;		//!< flag for relative bc

	vector<ITEM>	m_item;		//!< item list

	DECLARE_FECORE_CLASS();
};
