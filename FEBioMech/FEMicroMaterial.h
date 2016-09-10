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
	mat3ds		m_S;				// LTE - 2nd Piola-Kirchhoff stress
	mat3d		m_F_prev;			// deformation gradient from last time step

	double     m_macro_energy;	// LTE - Macroscopic strain energy
	double	   m_micro_energy;	// LTE - Volume-average of strain energy throughout the RVE solution
	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero) 

	double	   m_macro_energy_inc;	// LTE - Macroscopic strain energy increment
	double	   m_micro_energy_inc;	// LTE - Microscopic strain energy increment

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
	int		m_neid;					//!< element Id
	int		m_ngp;					//!< Gauss-point (one-based!)
	char	m_szfile[MAX_FILE];		//!< file name
	bool	m_bdebug;				//!< debug flag
	FERVEProbe*	m_probe;

	DECLARE_PARAMETER_LIST();
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
	char		m_szrve[256];	//!< filename for RVE file
	char		m_szbc[256];	//!< name of nodeset defining boundary
	int			m_bctype;		//!< periodic bc flag
	FERVEModel	m_mrve;			//!< the master RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

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
	FEVecPropertyT<FEMicroProbe>	m_probe;

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
