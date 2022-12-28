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
#include "FEBioMech/FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEPeriodicBoundary1O.h"
#include "FECore/FECallBack.h"
#include "FERVEModel.h"

class FERVEProbe;

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint : public FEElasticMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint();

	//! Initialize material point data
	void Init();

	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	//! create a shallow copy
	FEMaterialPointData* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	mat3ds		m_S;				// 2nd Piola-Kirchhoff stress
	mat3d		m_F_prev;			// deformation gradient from last time step

	double     m_macro_energy;	// Macroscopic strain energy
	double	   m_micro_energy;	// Volume-average of strain energy throughout the RVE solution
	double	   m_energy_diff;	// Difference between macro energy and volume averaged energy of RVE (should be zero) 

	double	   m_macro_energy_inc;	// Macroscopic strain energy increment
	double	   m_micro_energy_inc;	// Microscopic strain energy increment

	FERVEModel	m_rve;				// Local copy of the parent rve
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial :	public FEElasticMaterial
{
public:
	FEMicroMaterial(FEModel* pfem);
	~FEMicroMaterial(void);

public:
	std::string	m_szrve;	//!< filename for RVE file
	std::string	m_szbc;		//!< name of nodeset defining boundary
	int			m_bctype;		//!< periodic bc flag
	double		m_scale;		//!< RVE scale factor
	FERVEModel	m_mrve;			//!< the parent RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! data initialization
	bool Init() override;

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() override;

	// calculate the average PK1 stress
	mat3d AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp);

	// calculate the average PK2 stress
	mat3ds AveragedStressPK2(FEModel& rve, FEMaterialPoint &mp);

	// average RVE energy
	double micro_energy(FEModel& rve);

public:
	int Probes() { return (int) m_probe.size(); }
	FERVEProbe& Probe(int i) { return *m_probe[i]; }

protected:
	std::vector<FERVEProbe*>	m_probe;

public:
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
