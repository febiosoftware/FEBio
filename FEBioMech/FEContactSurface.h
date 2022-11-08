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

#include <FECore/FESurface.h>
#include <FECore/vec2d.h>
#include "FEContactInterface.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
// Stores material point data for contact interfaces
class FEBIOMECH_API FEContactMaterialPoint : public FESurfaceMaterialPoint
{
public:
	FEContactMaterialPoint()
	{
		m_gap = 0.0;
		m_Ln  = 0.0;
		m_pme = nullptr;
		m_pmep = nullptr;
	}

	void Init() override
	{
		FESurfaceMaterialPoint::Init();
		m_gap = 0.0;
		m_Ln = 0.0;
		m_pme = nullptr;
		m_pmep = nullptr;
	}

	void Serialize(DumpStream& ar) override
	{
		FESurfaceMaterialPoint::Serialize(ar);
		ar & m_gap & m_Ln;
	}

public:
	double	m_gap;	//!< gap function at integration points
	double	m_Ln;	//!< net contact pressure

	FESurfaceElement*	m_pme;	//!< target element
	FESurfaceElement*	m_pmep;	//!< previous target element
};

//-----------------------------------------------------------------------------
//! This class describes a contact surface

//!	this class is used in contact analyses to describe a contacting
//! surface in a contact interface.

class FEBIOMECH_API FEContactSurface : public FESurface
{
public:
	//! constructor
	FEContactSurface(FEModel* pfem);

	//! destructor
	~FEContactSurface();

	// initialization
	bool Init() override;

	// serialization
	void Serialize(DumpStream& ar) override;

	//! Set the sibling of this contact surface
	void SetSibling(FEContactSurface* ps);

    //! Set the parent of this contact surface
    void SetContactInterface(FEContactInterface* ps);
    
    //! Get the parent of this contact surface
    FEContactInterface* GetContactInterface() { return m_pContactInterface; }
    
	//! Unpack surface element data
	virtual void UnpackLM(FEElement& el, vector<int>& lm);

public:
    virtual void GetVectorGap      (int nface, vec3d& pg);
    virtual void GetContactTraction(int nface, vec3d& pt);
    
    virtual void GetNodalVectorGap      (int nface, vec3d* pg);
	virtual void GetNodalContactPressure(int nface, double* pg);
	virtual void GetNodalContactTraction(int nface, vec3d* pt);

    virtual void GetStickStatus(int nface, double& pt);
    
    void GetSurfaceTraction(int nface, vec3d& pt);
    void GetNodalSurfaceTraction(int nface, vec3d* pt);
    void GetGPSurfaceTraction(int nface, vec3d* pt);

	virtual vec3d GetContactForce();
    virtual double GetContactArea();

	FEModel* GetFEModel() { return m_pfem; }

protected:
	FEContactSurface* m_pSibling;
    FEContactInterface* m_pContactInterface;
	FEModel*	m_pfem;

	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
};
