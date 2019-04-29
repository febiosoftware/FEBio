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
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FEBioMech/FEBodyForce.h>
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
//! Fluid material point class.
//
class FEBIOFLUID_API FEFluidMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
    FEFluidMaterialPoint(FEMaterialPoint* pt = 0);

	//! create a shallow copy
	FEMaterialPoint* Copy();

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
    double      m_Jf;       //!< determinant of fluid deformation gradient
    double      m_Jfdot;    //!< material time derivative of Jf
    vec3d       m_gradJf;   //!< gradient of Jf
	double		m_pf;		//!< elastic fluid pressure
    mat3ds		m_sf;		//!< fluid stress
};

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEBIOFLUID_API FEFluid : public FEMaterial
{
public:
	FEFluid(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;
	
public:
    //! initialization
    bool Init() override {
        m_Tr = GetFEModel()->GetGlobalConstant("T");
        return true;
    }
    
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp);
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp)  { return m_pViscous->Tangent_RateOfDeformation(mp); }
    
    //! elastic pressure
    double Pressure(FEMaterialPoint& mp);
    virtual double Pressure(const double e);

    //! tangent of elastic pressure with respect to strain J
    virtual double Tangent_Pressure_Strain(FEMaterialPoint& mp) { return -m_k; }
    
    //! 2nd tangent of elastic pressure with respect to strain J
    virtual double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) { return 0; }
    
	//! referential fluid density
	double ReferentialDensity() { return m_rhor; }

    //! calculate current fluid density
    double Density(FEMaterialPoint& pt);
    
    //! return viscous part
    FEViscousFluid* GetViscous() { return m_pViscous; }
    
    //! kinematic viscosity
    double KinematicViscosity(FEMaterialPoint& mp);
    
    //! acoustic speed
    double AcousticSpeed(FEMaterialPoint& mp);
    
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp);
    
    //! strain energy density
    virtual double StrainEnergyDensity(FEMaterialPoint& mp);
    
    //! kinetic energy density
    double KineticEnergyDensity(FEMaterialPoint& mp);
    
    //! strain + kinetic energy density
    double EnergyDensity(FEMaterialPoint& mp);
    
    //! invert pressure-dilatation relation
    virtual double Dilatation(const double p);
    
    //! evaluate temperature
    virtual double Temperature(FEMaterialPoint& mp) { return m_Tr; }
    
private: // material properties
    FEViscousFluid*		m_pViscous; //!< pointer to viscous part of fluid material
	
public:
    double						m_rhor;     //!< referential fluid density
    double                      m_k;        //!< bulk modulus at J=1
    double                      m_Tr;       //!< ambient temperature
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
