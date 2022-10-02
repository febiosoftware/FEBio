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
#include "FEFluidFSI.h"
#include <FEBioMix/FEHydraulicPermeability.h>
#include <FEBioMix/FEBiphasic.h>
#include <FEBioMech/FEBodyForce.h>

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEBiphasicFSIMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEBiphasicFSIMaterialPoint(FEMaterialPointData* pt);
    
    //! create a shallow copy
	FEMaterialPointData* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // Biphasic FSI material data
    mat3d       m_Lw;       //!< grad of m_wt
    vec3d       m_gradJ;    //!< gradient of J
    double      m_phi0;     //!< solid volume fraction in reference configuration
    mat3ds      m_ss;       //!< solid stress
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEBiphasicFSI : public FEFluidFSI, public FEBiphasicInterface
{
public:
    FEBiphasicFSI(FEModel* pfem);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
    
    //! performs initialization
    bool Init() override;
    
public:
    //! calculate inner stress at material point
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! return the permeability tensor as a matrix
    void Permeability(double k[3][3], FEMaterialPoint& pt);
    
    //! return the permeability as a tensor
    mat3ds Permeability(FEMaterialPoint& pt);
    
    //! return the inverse permeability as a tensor
    mat3ds InvPermeability(FEMaterialPoint& pt);
    
    //! return the tangent permeability tensor
    tens4dmm Permeability_Tangent(FEMaterialPoint& pt);
    
    //! return the permeability property
    FEHydraulicPermeability* GetPermeability() { return m_pPerm; }
    
    //! porosity
    double Porosity(FEMaterialPoint& pt);
    
    //! Solid Volume
    double SolidVolumeFrac(FEMaterialPoint& pt);
    
    //! porosity gradient
    vec3d gradPorosity(FEMaterialPoint& pt);
    
    //! porosity gradient
    vec3d gradPhifPhis(FEMaterialPoint& pt);
    
    //! solid density
    double TrueSolidDensity(FEMaterialPoint& mp) { return Solid()->Density(mp); }
    
    //! true fluid density
    double TrueFluidDensity(FEMaterialPoint& mp) { return Fluid()->Density(mp); }
    
    //! solid density
    double SolidDensity(FEMaterialPoint& mp);
    
    //! fluid density
    double FluidDensity(FEMaterialPoint& mp);
    
public: // overridden from FEBiphasicInterface

    double GetReferentialSolidVolumeFraction(const FEMaterialPoint& mp) override {
        const FEBiphasicFSIMaterialPoint* pt = (mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        return pt->m_phi0;
    }

    //! solid referential apparent density
    double SolidReferentialApparentDensity(FEMaterialPoint& pt) override;

    //! solid referential volume fraction
    double SolidReferentialVolumeFraction(FEMaterialPoint& pt) override;

public: // material parameters
    double      m_rhoTw; //!< true fluid density
    FEParamDouble      m_phi0;  //!< solid volume fraction in reference configuration

    vector<FEBodyForce*>    m_bf;       //!< body forces acting on this biphasic material
    
protected: // material properties
    FEHydraulicPermeability*    m_pPerm;    //!< pointer to permeability material
    
    DECLARE_FECORE_CLASS();
};
