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
#include "FEHydraulicPermeability.h"
#include "FESolventSupply.h"
#include "FEActiveMomentumSupply.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Biphasic material point class.
//
class FEBIOMIX_API FEBiphasicMaterialPoint : public FEMaterialPointData
{
public:
	//! constructor
	FEBiphasicMaterialPoint(FEMaterialPointData* ppt);

	//! create a shallow copy
	FEMaterialPointData* Copy() override;

	//! data serialization
	void Serialize(DumpStream& ar) override;

	//! Data initialization
	void Init() override;

public:
	// poro-elastic material data
	// The actual fluid pressure is the same as the effective fluid pressure
	// in a poroelastic material without solute(s).  The actual fluid pressure
	// is included here so that models that include both poroelastic and
	// solute-poroelastic domains produce plotfiles with consistent fluid
	// pressure fields.
	double		m_p;		//!< fluid pressure
	vec3d		m_gradp;	//!< spatial gradient of p
    vec3d		m_gradpp;	//!< gradp at previous time
	vec3d		m_w;		//!< fluid flux
	double		m_pa;		//!< actual fluid pressure
    double      m_phi0;     //!< referential solid volume fraction at initial time
	double		m_phi0t;	//!< referential solid volume fraction at current time
	double		m_phi0p;	//!< referential solid volume fraction at previous time
	double		m_phi0hat;	//!< referential solid volume fraction supply at current time
    double      m_Jp;       //!< determinant of solid deformation gradient at previous time
    mat3ds      m_ss;       //!< solid (elastic or effective) stress
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FEBiphasicInterface
{
public:
	FEBiphasicInterface() {}
	virtual ~FEBiphasicInterface() {}

public:
	virtual double GetReferentialSolidVolumeFraction(const FEMaterialPoint& mp) { return 0.0; }

	// TODO: These are only used by multiphasic materials. Perhapse move to separate interface class? 
	//! solid referential apparent density
	virtual double SolidReferentialApparentDensity(FEMaterialPoint& pt) { return 0.0; }

	//! solid referential volume fraction
	virtual double SolidReferentialVolumeFraction(FEMaterialPoint& pt) { return 0.0; };

	// TODO: This is a bit of a hack to get the fluid pressure. 
	virtual double GetActualFluidPressure(const FEMaterialPoint& pt) { return 0.0; }
};

//-----------------------------------------------------------------------------
//! Base class for biphasic materials.

class FEBIOMIX_API FEBiphasic : public FEMaterial, public FEBiphasicInterface
{
public:
	FEBiphasic(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }
	
public:
    //! initialize
    bool Init() override;
    
    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! return the permeability tensor as a matrix
	void Permeability(double k[3][3], FEMaterialPoint& pt);

	//! return the permeability as a tensor
	mat3ds Permeability(FEMaterialPoint& pt);

	//! return the permeability property
	FEHydraulicPermeability* GetPermeability() { return m_pPerm; }
	
	//! calculate actual fluid pressure
	double Pressure(FEMaterialPoint& pt);

	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
    //! solid density
    double SolidDensity(FEMaterialPoint& mp) { return m_pSolid->Density(mp); }
    
	//! fluid density
	double FluidDensity() { return m_rhoTw; }

	//! get the solvent supply
	double SolventSupply(FEMaterialPoint& mp) { return (m_pSupp? m_pSupp->Supply(mp) : 0); }

	//! get the solvent supply property
	FESolventSupply* GetSolventSupply() { return m_pSupp; }

	//! Get the active momentum supply
	FEActiveMomentumSupply* GetActiveMomentumSupply() { return m_pAmom; }

public: // overridden from FEBiphasicInterface

	double GetReferentialSolidVolumeFraction(const FEMaterialPoint& mp) override {
		const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
		return pt->m_phi0t;
	}

	double GetActualFluidPressure(const FEMaterialPoint& mp) override { 
		const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
		return pt->m_pa;
	}

    //! evaluate and return solid referential volume fraction
    double SolidReferentialVolumeFraction(FEMaterialPoint& mp) override {
        double phisr = m_phi0(mp);
        FEBiphasicMaterialPoint* bp = (mp.ExtractData<FEBiphasicMaterialPoint>());
        bp->m_phi0 = bp->m_phi0t = phisr;
        return phisr;
    };
    
public: // material parameters
	double						m_rhoTw;	//!< true fluid density
	FEParamDouble               m_phi0;		//!< solid volume fraction in reference configuration
    double                      m_tau;      //!< characteristic time constant for stabilization
    vector<FEBodyForce*>        m_bf;       //!< body forces acting on this biphasic material

private: // material properties
	FEElasticMaterial*			m_pSolid;	//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;	//!< pointer to permeability material
	FESolventSupply*			m_pSupp;	//!< pointer to solvent supply
	FEActiveMomentumSupply*		m_pAmom;	//!< pointer to active momentum supply
	
	DECLARE_FECORE_CLASS();
};
