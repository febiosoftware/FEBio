/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
//! This class describes a contact slave or master surface used for 
//! tied contact

//!	this class is used in contact analyses to describe a contacting
//! surface in a tied contact interface.

class FETiedContactSurface : public FEContactSurface
{
public:
	//! constructor
	FETiedContactSurface(FEModel* pfem) : FEContactSurface(pfem) { m_boffset = false; }

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
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers
	vector<vec3d>				m_Tc;	//!< contact forces
	vector<double>				m_off;	//!< offset values (used for shells)

protected:
	bool	m_boffset;		//!< offset shells
};
