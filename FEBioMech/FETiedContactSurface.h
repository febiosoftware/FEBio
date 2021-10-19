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
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class describes a contact surface used for tied contact.

//!	this class is used in contact analyses to describe a contacting
//! surface in a tied contact interface.

class FETiedContactSurface : public FEContactSurface
{
public:
	class Data : public FEContactMaterialPoint
	{
	public:
		Data();
		void Serialize(DumpStream& ar);

	public:
		vec3d				m_vgap;	//!< gap function at nodes
		vec2d				m_rs;	//!< natural coordinates of projection on secondary surface element
		vec3d				m_Lm;	//!< Lagrange multipliers
		vec3d				m_Tc;	//!< contact forces
		double				m_off;	//!< offset values (used for shells)
	};

public:
	//! constructor
	FETiedContactSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! data serialization
	void Serialize(DumpStream& ar);

	//! offset shell surfaces (must be called before Init())
	void SetShellOffset(bool b) { m_boffset = b; }

public:
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactPressure(int nface, double* pn);
	void GetNodalContactTraction(int nface, vec3d* tn);

public:
	vector<Data>	m_data;	//!< integration point data

protected:
	bool	m_boffset;		//!< offset shells
};
