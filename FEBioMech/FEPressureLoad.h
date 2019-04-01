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
#include "FESurfaceTraction.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceTraction
{
public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	//! serialize data
	void Serialize(DumpStream& ar) override;

protected:
	//! calculate traction
	vec3d Traction(const FESurfaceMaterialPoint& mp) override;

	//! calculate stiffness for an element
	void ElementStiffness(FESurfaceElement& el, matrix& ke) override;

	void SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);
	void UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);

protected:
	FEParamDouble	m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation

	DECLARE_FECORE_CLASS();
};
