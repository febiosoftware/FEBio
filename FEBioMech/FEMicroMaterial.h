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
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEPeriodicBoundary1O.h"
#include "FECore/FECallBack.h"
#include "FERVEModel.h"

//-----------------------------------------------------------------------------
class FEBioPlotFile;

//-----------------------------------------------------------------------------
class FERVEProbe : public FECallBack
{
public:
	// The first FEModel (fem) is the macro-problem, i.e. the model that will generate the callbacks
	// The second FEModel (rve) is the micro-problem that needs to be tracked.
	FERVEProbe(FEModel& fem, FEModel& rve, const char* szfile);

	bool Execute(FEModel& fem, int nwhen);

	void Save();

	void SetDebugFlag(bool b) { m_bdebug = b; }
	bool GetDebugFlag() const { return m_bdebug; }

private:
	FEModel&			m_rve;		//!< The RVE model to keep track of
	FEBioPlotFile*		m_xplt;		//!< the actual plot file
	std::string			m_file;		//!< file name
	bool				m_bdebug;
};

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint(FEMaterialPoint* mp);

	//! Initialize material point data
	void Init();

	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	//! create a shallow copy
	FEMaterialPoint* Copy();

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

	FERVEModel	m_rve;				// Local copy of the master rve
};

//-----------------------------------------------------------------------------
// The FEMicroProbe class is not really a material, but we abuse the framework
// here in order to read in the probe information. 
class FEMicroProbe : public FEMaterial
{
	enum { MAX_FILE = 128 };

public:
	FEMicroProbe(FEModel* pfem);
	~FEMicroProbe();

public:
	int			m_neid;			//!< element Id
	int			m_ngp;			//!< Gauss-point (one-based!)
	std::string	m_szfile;		//!< file name
	bool		m_bdebug;		//!< debug flag
	FERVEProbe*	m_probe;

	DECLARE_FECORE_CLASS();
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
	FERVEModel	m_mrve;			//!< the master RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! data initialization
	bool Init() override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() override;

	// calculate the average PK1 stress
	mat3d AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp);

	// calculate the average PK2 stress
	mat3ds AveragedStressPK2(FEModel& rve, FEMaterialPoint &mp);

	// average RVE energy
	double micro_energy(FEModel& rve);

public:
	int Probes() { return (int) m_probe.size(); }
	FEMicroProbe& Probe(int i) { return *m_probe[i]; }

protected:
	std::vector<FEMicroProbe*>	m_probe;

public:
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
