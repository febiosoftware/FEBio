/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include "FEDomain.h"
#include "FEShellElement.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell element domains
class FECORE_API FEShellDomain : public FEDomain
{
	FECORE_SUPER_CLASS(FESHELLDOMAIN_ID)
	FECORE_BASE_CLASS(FEShellDomain)

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
	// evaluate volume of element in reference frame
	virtual double Volume(FEShellElement& el) { return 0.0; }

	// evaluate volume of element in current frame
	virtual double CurrentVolume(FEShellElement& el) { return 0.0; }

	// Initialize shell data (Called from FEMesh::InitShells)
	virtual void InitShells();

public:
    //! get the current nodal coordinates
    void GetCurrentNodalCoordinates(const FEShellElement& el, vec3d* rt, const bool back = false);
    void GetCurrentNodalCoordinates(const FEShellElement& el, vec3d* rt, double alpha, const bool back = false);
    
    //! get the reference nodal coordinates
    void GetReferenceNodalCoordinates(const FEShellElement& el, vec3d* r0, const bool back = false);
    
    //! get the nodal coordinates at previous state
    void GetPreviousNodalCoordinates(const FEShellElement& el, vec3d* rp, const bool back = false);
    
public:
	void ForEachShellElement(std::function<void(FEShellElement& el)> f);
};

//-----------------------------------------------------------------------------
// Old director-based shell formulation
class FECORE_API FEShellDomainOld : public FEShellDomain
{
public:
	FEShellDomainOld(FEModel* fem);

	//! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

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
	bool Create(int nsize, FE_Element_Spec espec) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEShellElementNew& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

	double DefaultShellThickness() const { return m_h0; }

	void AssignDefaultShellThickness();

protected:
	double	m_h0;

protected:
	vector<FEShellElementNew>	m_Elem;	//!< array of elements

	DECLARE_FECORE_CLASS();
};
