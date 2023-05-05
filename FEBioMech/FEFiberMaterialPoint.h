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
#include "FECore/FEMaterial.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
// Define a material point that stores the fiber pre-stretch
class FEBIOMECH_API FEFiberMaterialPoint : public FEMaterialPointData
{
public:
	FEFiberMaterialPoint(FEMaterialPointData* pt);
    
	FEMaterialPointData* Copy();
    
    void Init();
    
    void Serialize(DumpStream& ar);
    
public:
    // Set or clear pre-stretch, as needed in multigenerational materials (e.g., reactive viscoelasticity)
    void SetPreStretch(const mat3ds Us) { m_Us = Us; m_bUs = true; }
    void ResetPreStretch() { m_bUs = false; }
    vec3d FiberPreStretch(const vec3d a0);

public:
    mat3ds  m_Us;   //!< pre-stretch tensor for fiber
    bool    m_bUs;  //!< flag for pre-stretch
};
