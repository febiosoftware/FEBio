#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FEShellDomain : public FEDomain
{
public:
	//! constructor
	FEShellDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_SHELL, pm) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	//! Initialize elements
	void InitElements();

	//! Reset element data
	void Reset();

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElement& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

	//! Serialize domain data to archive
	void Serialize(DumpStream& ar);

protected:
	vector<FEShellElement>	m_Elem;	//!< array of elements
};
