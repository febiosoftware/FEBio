#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FECORE_API FEDiscreteDomain : public FEDomain
{
public:
	FEDiscreteDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_DISCRETE, pm) {}

	void Create(int nelems, int elemType) override;
	int Elements() const override { return (int) m_Elem.size(); }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEDiscreteElement& Element(int n) { return m_Elem[n]; }

	bool Init() override;

public:
	void AddElement(int eid, int n[2]);

protected:
	vector<FEDiscreteElement>	m_Elem;
};
