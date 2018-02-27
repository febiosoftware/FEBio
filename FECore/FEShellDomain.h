#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell element domains
class FEShellDomain: public FEDomain
{
public:
	//! constructor
	FEShellDomain(FEMesh* pm);

	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! Reset element data
	void Reset();

	// get a shell element
	virtual FEShellElement& Element(int i) = 0;

	// get the element type (TODO: Move to FEDomain class?)
	int GetElementType() { return ElementRef(0).Type(); };
};

//-----------------------------------------------------------------------------
// Old director-based shell formulation
class FEShellDomainOld : public FEShellDomain
{
public:
	FEShellDomainOld(FEMesh* pm) : FEShellDomain(pm) {}

	//! create storage for elements
	void Create(int nsize, int elemType);

	//! Serialize domain data to archive
	void Serialize(DumpStream& ar);

public:
	//! return nr of elements
	int Elements() const { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	FEShellElementOld& ShellElement(int i) { return m_Elem[i]; }

protected:
	vector<FEShellElementOld>	m_Elem;	//!< array of elements
};

//-----------------------------------------------------------------------------
// New shell formulation
class FEShellDomainNew : public FEShellDomain
{
public:
	FEShellDomainNew(FEMesh* pm) : FEShellDomain(pm) {}

	//! create storage for elements
	void Create(int nsize, int elemType);

	//! Serialize domain data to archive
	void Serialize(DumpStream& ar);

public:
	//! return nr of elements
	int Elements() const { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

protected:
	vector<FEShellElement>	m_Elem;	//!< array of elements
};
