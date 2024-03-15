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
#include <FECore/FEMaterialPoint.h>
#include <FECore/tens4d.h>
#include <FECore/FEMaterialPointProperty.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FEBIOMECH_API FEElasticMaterialPoint : public FEMaterialPointData
{
public:
	//! constructor
	FEElasticMaterialPoint(FEMaterialPointData* mp = nullptr);

	//! Initialize material point data
	void Init() override;

	//! create a shallow copy
	FEMaterialPointData* Copy() override;

	//! serialize material point data
	void Serialize(DumpStream& ar) override;

public:
	mat3ds Strain() const;
	mat3ds SmallStrain() const;

	mat3ds RightCauchyGreen() const;
	mat3ds LeftCauchyGreen () const;

	mat3ds DevRightCauchyGreen() const;
	mat3ds DevLeftCauchyGreen () const;
    
    mat3ds RightStretch() const;
    mat3ds LeftStretch () const;
    
    mat3ds RightStretchInverse() const;
    mat3ds LeftStretchInverse () const;
    
    mat3ds RightHencky() const;
    mat3ds LeftHencky () const;
    
    mat3d Rotation() const;
    
    mat3ds RateOfDeformation() const { return m_L.sym(); }

	mat3ds pull_back(const mat3ds& A) const;
	mat3ds push_forward(const mat3ds& A) const;

	tens4ds pull_back(const tens4ds& C) const;
	tens4ds push_forward(const tens4ds& C) const;

public:
    bool    m_buncoupled;   //!< set to true if this material point was created by an uncoupled material
    
	// deformation data at intermediate time
	mat3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
	double	m_J_star;	//!< desired J
    vec3d   m_gradJ;  //!< gradient of J
    vec3d   m_v;    //!< velocity
    vec3d   m_a;    //!< acceleration
    mat3d   m_L;    //!< spatial velocity gradient

	// solid material data
	mat3ds		m_s;		//!< Cauchy stress
    
    // uncoupled pressure
    double      m_p;        //!< only for uncoupled materials
    
    // current time data
    double      m_Wt;       //!< strain energy density at current time
    
    // previous time data
    double      m_Wp;       //!< strain energy density
    
    // augmented Lagrangian variable
    double      m_Lk;       //!< strain energy density
};
