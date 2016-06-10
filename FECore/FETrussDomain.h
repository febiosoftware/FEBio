#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for truss elements
class FETrussDomain : public FEDomain
{
public:
	FETrussDomain(FEMesh* pm);

public:
	void Create(int nsize, int elemType);

	int Elements() const { return (int)m_Elem.size(); }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) { return m_Elem[n]; }

public:
	//! Calculate the truss normal
	vec3d TrussNormal(FETrussElement& el);

protected:
	vector<FETrussElement>	m_Elem;
};
