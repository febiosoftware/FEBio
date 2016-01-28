#pragma once
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEPeriodicBoundary1O.h"
#include "FECore/FECallBack.h"

//-----------------------------------------------------------------------------
class FEBioPlotFile;

//-----------------------------------------------------------------------------
class FERVEProbe : public FECallBack
{
public:
	// The first FEModel (fem) is the macro-problem, i.e. the model that will generate the callbacks
	// The second FEModel (rve) is the micro-problem that needs to be tracked.
	FERVEProbe(FEModel& fem, FEModel& rve, const char* szfile);

	void Execute(FEModel& fem, int nwhen);

private:
	FEModel&			m_rve;		//!< The RVE model to keep track of
	FEBioPlotFile*		m_xplt;		//!< the actual plot file
	std::string			m_file;		//!< file name
};

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint(FEMaterialPoint* mp);

	//! Initialize material point data
	void Init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	mat3d      m_PK1;			// LTE - 1st Piola-Kirchhoff stress
	mat3ds     m_S;				// LTE - 2nd Piola-Kirchhoff stress
	mat3d      m_PK1_prev;

	mat3ds     m_inf_str;		// LTE - infinitesimal strain
	mat3ds     m_E;				// LTE - Green-Lagrange strain
	mat3ds     m_e;				// LTE - Euler-Almansi strain
	
	double     m_macro_energy;	// LTE - Macroscopic strain energy
	double	   m_micro_energy;	// LTE - Volume-average of strain energy throughout the RVE solution
	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero) 

	double	   m_macro_energy_inc;	// LTE - Macroscopic strain energy increment
	double	   m_micro_energy_inc;	// LTE - Microscopic strain energy increment

	bool    m_rve_init;			// LTE - Flag indicating that the rve has been initialized
	FEModel m_rve;				// LTE - Current copy of the rve		
	FEModel m_rve_prev;			// LTE - Previous converged state of the rve

	tens4ds	m_Ka;				// LTE - Averaged rank 4 material stiffness
};

//-----------------------------------------------------------------------------
// The FEMicroProbe class is not really a material, but we abuse the framework
// here in order to read in the probe information. 
class FEMicroProbe : public FEMaterial
{
	enum { MAX_FILE = 128 };

public:
	FEMicroProbe(FEModel* pfem);

public:
	int		m_neid;					//!< element Id
	int		m_ngp;					//!< Gauss-point (one-based!)
	char	m_szfile[MAX_FILE];		//!< file name

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
	char	m_szrve[256];	//!< filename for RVE file
	char	m_szbc[256];	//!< name of nodeset defining boundary
	bool	m_bperiodic;	//!< periodic bc flag
	FEModel	m_mrve;			//!< the master RVE (Representive Volume Element)

protected:
	double	m_V0;			//!< initial volume of RVE
	vector<int> m_BN;		//!< boundary node flags

	double m_bb_x; double m_bb_y; double m_bb_z;  // LTE - RVE bounding box
	int m_num_ext_node;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	bool PrepRVE();
	bool PrepDisplacementBC();
	bool PrepPeriodicBC();
	void FindBoundaryNodes();

	double InitVolume();

	void UpdateBC(FEModel& rve, mat3d& F);
	
	mat3ds AveragedStress(FEModel& rve, FEMaterialPoint& pt);
	tens4ds AveragedStiffness(FEModel& rve, FEMaterialPoint& pt);

	mat3d AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp);
	mat3ds AveragedStressPK2(FEModel& rve, FEMaterialPoint &mp);

	void calc_energy_diff(FEMaterialPoint& pt);

public:
	int Probes() { return (int) m_probe.size(); }
	FEMicroProbe& Probe(int i) { return *m_probe[i]; }

protected:
	FEVecPropertyT<FEMicroProbe>	m_probe;

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
