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
#include <FEBioMech/FEElasticMaterial.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEFSIMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEFSIMaterialPoint(FEMaterialPointData* pt);
    
    //! create a shallow copy
	FEMaterialPointData* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // FSI material data
    vec3d       m_w;      //!< fluid flux relative to solid
    vec3d       m_aw;     //!< material time derivative of m_wt
    double      m_Jdot;   //!< time derivative of solid volume ratio
    mat3ds      m_ss;     //!< solid stress
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEFluidFSI : public FEMaterial
{
public:
    FEFluidFSI(FEModel* pfem);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
    
    // Get the elastic component (overridden from FEMaterial)
    FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFluid* Fluid() { return m_pFluid; }
    FEElasticMaterial* Solid() { return m_pSolid; }
    
protected: // material properties
    FEElasticMaterial*			m_pSolid;	//!< pointer to elastic solid material
    FEFluid*                    m_pFluid;	//!< pointer to fluid material
    
    DECLARE_FECORE_CLASS();
};
