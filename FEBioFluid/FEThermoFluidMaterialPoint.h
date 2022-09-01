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
class FEBIOFLUID_API FEThermoFluidMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
	FEThermoFluidMaterialPoint(FEMaterialPointData* pt);

    //! create a shallow copy
	FEMaterialPointData* Copy();

    //! data serialization
    void Serialize(DumpStream& ar);

    //! Data initialization
    void Init();

public:
    // fluid data
    double      m_T;        //!< temperature (relative to reference temperature)
    double      m_Tdot;     //!< material time derivative of temperature
    vec3d       m_gradT;    //!< temperature gradient
    double      m_k;        //!< bulk modulus
    double      m_K;        //!< thermal conductivity
    double      m_dKJ;      //!< derivative of thermal conductivity with respect to fluid volume ratio
    double      m_dKT;      //!< derivative of thermal conductivity with respect to temperature
    double      m_cv;       //!< isochoric specific heat capacity
    double      m_dcvJ;     //!< derivative of isochoric specific heat capacity w.r.t. fluid volume ratio
    double      m_dcvT;     //!< derivative of isochoric specific heat capacity w.r.t. temperature
    double      m_cp;       //!< isobaric specific heat capacity
    vec3d       m_q;        //!< heat flux
};

