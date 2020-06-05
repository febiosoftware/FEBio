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
#include <FEBioMix/FEHydraulicPermeability.h>
#include <FEBioMech/FEBodyForce.h>

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEBiphasicFSIMaterialPoint : public FEMaterialPoint
{
public:
    //! constructor
    FEBiphasicFSIMaterialPoint(FEMaterialPoint* pt);
    
    //! create a shallow copy
    FEMaterialPoint* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // Biphasic FSI material data
    vec3d       m_w;      //!< fluid flux relative to solid
    vec3d       m_aw;     //!< material time derivative of m_wt
    mat3d       m_Lw;     //!< grad of m_wt
    double      m_Jdot;   //!< time derivative of solid volume ratio
    double      m_phis;   //!< solid volume fraction
    double      m_phif;   //!< fluid volume fraction
    vec3d       m_gradphif;   //!< gradient of fluid volume fraction
    vec3d       m_gradJ;      //!< gradient of J
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEBiphasicFSI : public FEMaterial
{
public:
    FEBiphasicFSI(FEModel* pfem);
    
    // returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
    // Get the elastic component (overridden from FEMaterial)
    FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }
    
    //! performs initialization
    bool Init() override;
    
public:
    //! calculate inner stress at material point
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! calculate inner tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt);
    
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
    
    //! solid density
    double TrueSolidDensity(FEMaterialPoint& mp) { return m_pSolid->Density(mp); }
    
    //! true fluid density
    double TrueFluidDensity(FEMaterialPoint& mp) { return m_pFluid->Density(mp); }
    
    //! solid density
    double SolidDensity(FEMaterialPoint& mp);
    
    //! fluid density
    double FluidDensity(FEMaterialPoint& mp);
    
    FEFluid* Fluid() { return m_pFluid; }
    FEElasticMaterial* Solid() { return m_pSolid; }
    
public: // material parameters
    double      m_rhoTw; //!< true fluid density
    double      m_phi0;  //!< solid volume fraction in reference configuration

    vector<FEBodyForce*>    m_bf;       //!< body forces acting on this biphasic material
    
private: // material properties
    FEElasticMaterial*            m_pSolid;    //!< pointer to elastic solid material
    FEFluid*                    m_pFluid;    //!< pointer to fluid material
    FEHydraulicPermeability*    m_pPerm;    //!< pointer to permeability material
    
    DECLARE_FECORE_CLASS();
};
