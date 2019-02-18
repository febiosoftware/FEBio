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
