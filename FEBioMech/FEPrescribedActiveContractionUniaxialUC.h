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
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of an uncoupled solid mixture material.
class FEPrescribedActiveContractionUniaxialUC : public FEUncoupledMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionUniaxialUC(FEModel* pfem);
    
    //! Validation
    bool Validate() override;

	//! serialization
	void Serialize(DumpStream& ar) override;
    
    //! stress
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! tangent
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
public:
    FEParamDouble	m_T0;       // prescribed active stress

private:
	vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of an uncoupled solid mixture material.
class FEPrescribedActiveContractionFiberUC : public FEUncoupledMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionFiberUC(FEModel* pfem);

    //! stress
    mat3ds DevStress(FEMaterialPoint& pt) override;

    //! tangent
    tens4ds DevTangent(FEMaterialPoint& pt) override;

private:
    FEParamDouble	m_T0;       // prescribed active stress
    FEParamVec3     m_n0;		// unit vector along fiber direction (local coordinate system)

    DECLARE_FECORE_CLASS();
};
