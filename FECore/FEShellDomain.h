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

public:
	// evaluate volume of element
	virtual double Volume(FEShellElement& el) { return 0.0; }

	// Initialize shell data (Called from FEMesh::InitShells)
	virtual void InitShells();
};

//-----------------------------------------------------------------------------
// Old director-based shell formulation
class FEShellDomainOld : public FEShellDomain
{
public:
	FEShellDomainOld(FEMesh* pm);

	//! create storage for elements
	void Create(int nsize, int elemType) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }

	FEShellElementOld& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

	void InitShells() override;

protected:
	vector<FEShellElementOld>	m_Elem;	//!< array of elements
};

//-----------------------------------------------------------------------------
// New shell formulation
class FEShellDomainNew : public FEShellDomain
{
public:
	FEShellDomainNew(FEMesh* pm);

	//! create storage for elements
	void Create(int nsize, int elemType) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }

	FEShellElementNew& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

protected:
	vector<FEShellElementNew>	m_Elem;	//!< array of elements
};
