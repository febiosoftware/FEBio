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

//-----------------------------------------------------------------------------
//! Material point data for visco-elastic materials
class FEViscoElasticMaterialPoint : public FEMaterialPointData
{
public:
	enum { MAX_TERMS = 6 };

public:
	//! constructor
	FEViscoElasticMaterialPoint(FEMaterialPointData* mp = nullptr);

	//! copy material point data
	FEMaterialPointData* Copy();

	//! Initialize material point data
	void Init();

	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	//! Serialize data to archive
	void Serialize(DumpStream& ar);

public:
	mat3ds	m_Se;	//!< elastic PK2 stress
	mat3ds	m_Sep;	//!< elastic 2nd PK stress at previous time

	mat3ds	m_H[MAX_TERMS];		//!< internal variables
	mat3ds	m_Hp[MAX_TERMS];	//!< internal variables at previous timestep

    double  m_alpha[MAX_TERMS];     //!< exponent of right-stretch tensor in series spring
    double  m_alphap[MAX_TERMS];    //!< alpha at previous time step
    
	double	m_sed;	//!< elastic strain energy density
	double	m_sedp;	//!< elastic strain energy density at previous time
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FEViscoElasticMaterial :	public FEElasticMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };

public:
	//! default constructor
	FEViscoElasticMaterial(FEModel* pfem);

	//! get the elastic base material
	FEElasticMaterial* GetBaseMaterial();

	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase);

public:
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! strain energy density
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! calculate exponent of right-stretch tensor in series spring
    bool SeriesStretchExponent(FEMaterialPoint& pt);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

public: 
	// material parameters
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	double	m_t[MAX_TERMS];	//!< relaxation times

private:
	FEElasticMaterial*	m_Base;	//!< pointer to elastic solid material

public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
