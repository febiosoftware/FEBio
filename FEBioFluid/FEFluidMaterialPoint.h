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
#include <FECore/FEMaterial.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Fluid material point class.
//
class FEBIOFLUID_API FEFluidMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEFluidMaterialPoint(FEMaterialPointData* pt = 0);

    //! create a shallow copy
	FEMaterialPointData* Copy();

    //! data serialization
    void Serialize(DumpStream& ar);

    //! Data initialization
    void Init();

public:
    mat3ds RateOfDeformation() const { return m_Lf.sym(); }
    mat3da Spin() { return m_Lf.skew(); }
    vec3d  Vorticity() const { return vec3d(m_Lf(2,1)-m_Lf(1,2), m_Lf(0,2)-m_Lf(2,0), m_Lf(1,0)-m_Lf(0,1)); }
    
public:
    // fluid data
    vec3d       m_r0;       //!< material position
    vec3d       m_vft;      //!< fluid velocity
    vec3d       m_aft;      //!< fluid acceleration
    mat3d       m_Lf;       //!< fluid velocity gradient
    double      m_ef;       //!< fluid dilatation
    double      m_efdot;    //!< material time derivative of ef
    vec3d       m_gradef;   //!< gradient of ef
    double      m_pf;       //!< elastic fluid pressure
    mat3ds      m_sf;       //!< fluid stress
};
