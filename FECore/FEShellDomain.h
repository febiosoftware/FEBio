#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell element domains
class FECORE_API FEShellDomain : public FEDomain
{
public:
	//! constructor
	FEShellDomain(FEModel* fem);

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
class FECORE_API FEShellDomainOld : public FEShellDomain
{
public:
	FEShellDomainOld(FEModel* fem);

	//! create storage for elements
	void Create(int nsize, int elemType) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEShellElementOld& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

	void InitShells() override;

protected:
	vector<FEShellElementOld>	m_Elem;	//!< array of elements
};

//-----------------------------------------------------------------------------
// New shell formulation
class FECORE_API FEShellDomainNew : public FEShellDomain
{
public:
	FEShellDomainNew(FEModel* fem);

	//! create storage for elements
	void Create(int nsize, int elemType) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEShellElementNew& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

protected:
	vector<FEShellElementNew>	m_Elem;	//!< array of elements
};
