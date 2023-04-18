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
#include "FEElasticMaterial.h"
#include <FECore/FEMaterialPoint.h>

//-----------------------------------------------------------------------------
// Define a material point that stores the Donnan equilibrium variables.
class FEDonnanEquilibriumMaterialPoint : public FEMaterialPointData
{
public:
    FEDonnanEquilibriumMaterialPoint(FEMaterialPointData*pt) : FEMaterialPointData(pt) {}
    
	FEMaterialPointData* Copy();
    
    void Init();
    
    void Serialize(DumpStream& ar);
    
public:
    double      m_cFr;      //!< referential fixed-charge density
    double      m_cF;       //!< fixed-charge density
    double      m_osm;      //!< osmolarity
    double      m_p;        //!< osmotic pressure
    double      m_bpi;      //!< osmotic pressure tangent
};

//-----------------------------------------------------------------------------
//! Material class that implements Donnan equilibrium. 
//! When used on its own (not in a solid mixture), this materials
//! is intrinsically unstable
class FEDonnanEquilibrium : public FEElasticMaterial
{
public:
	// constructor
	FEDonnanEquilibrium(FEModel* pfem);
	
	//! Initialization routine
	bool Init() override;

	//! Returns the Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! Returs the spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

    //! Return the fixed-charge density
    double FixedChargeDensity(FEMaterialPoint& mp);
    
    //! Return the osmotic pressure
    double OsmoticPressure(FEMaterialPoint& mp);
    
    //! Return the osmotic pressure tangent
    double OsmoticPressureTangent(FEMaterialPoint& mp);
    
    //! Return the osmolarity
    double Osmolarity(FEMaterialPoint& mp);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override
    {
        return new FEDonnanEquilibriumMaterialPoint(new FEElasticMaterialPoint);
    }
    
    // update fatigue material point at each iteration
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
    double	m_phisr;	//!< referential solid volume fraction (may evolve with time)
	FEParamDouble	m_cFr;		//!< fixed charge density in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
    double  m_Phi;      //!< osmotic coefficient
    bool    m_bnew;     //!< flag for using old or new method
    bool    m_binit;    //!< initialization flag
};
