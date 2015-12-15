#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for truss elements
class FETrussDomain : public FEDomain
{
public:
	FETrussDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_TRUSS, pm){}

	void create(int n) { m_Elem.resize(n); }
	int Elements() { return m_Elem.size(); }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

	//! Calculate the truss normal
	vec3d TrussNormal(FETrussElement& el);

protected:
	vector<FETrussElement>	m_Elem;
};
