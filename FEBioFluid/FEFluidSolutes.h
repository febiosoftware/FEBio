/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include <FEBioMix/FESolute.h>
#include <FEBioMix/FESoluteInterface.h>

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEFluidSolutesMaterialPoint : public FEMaterialPoint
{
public:
    //! constructor
    FEFluidSolutesMaterialPoint(FEMaterialPoint* pt);
    
    //! create a shallow copy
    FEMaterialPoint* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // solutes material data
    int                 m_nsol;     //!< number of solutes
    vector<double>      m_c;        //!< solute concentration
    vector<vec3d>       m_gradc;    //!< spatial gradient of solute concentration
    vector<vec3d>       m_j;        //!< solute molar flux
    vector<double>      m_cdot;     //!< material time derivative of solute concentration following fluid
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEFluidSolutes : public FEMaterial, public FESoluteInterface
{
public:
    FEFluidSolutes(FEModel* pfem);
    
    // returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFluid* Fluid() { return m_pFluid; }
    
    //! calculate solute molar flux
    vec3d SoluteFlux(FEMaterialPoint& pt, const int sol);
    
    //! actual concentration (as opposed to effective concentration)
    double Concentration(FEMaterialPoint& pt, const int sol);
    
    // solute interface
public:
    int Solutes() override { return (int)m_pSolute.size(); }
    FESolute* GetSolute(int i) override { return m_pSolute[i]; }
    
private: // material properties
    FEFluid*                m_pFluid;       //!< pointer to fluid material
    std::vector<FESolute*>  m_pSolute;      //!< pointer to solute materials

    DECLARE_FECORE_CLASS();
};
