#pragma once
#include "FEElasticMaterial2O.h"
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>
#include "FE2OMicroConstraint.h"
#include "FEMicroMaterial.h"
#include "FERVEModel2O.h"

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint2O : public FEMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint2O(FEMaterialPoint* mp);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	FEMicroModel2O m_rve;				//!< local copy of the rve		
	int		m_elem_id;		//!< element ID
	int		m_gpt_id;		//!< Gauss point index (0-based)
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial2O : public FEElasticMaterial2O
{
public:
	FEMicroMaterial2O(FEModel* pfem);
	~FEMicroMaterial2O(void);

public:
	char			m_szrve[256];	//!< filename for RVE file
	char			m_szbc[256];	//!< name of nodeset defining boundary
	bool			m_bperiodic;	//!< periodic bc flag
	double			m_scale;		//!< geometry scale factor
	FERVEModel2O	m_mrve;			//!< the master RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	void Stress(FEMaterialPoint& mp, mat3d& P, tens3drs& Q);

	//! calculate tangent stiffness at material point
	void Tangent(FEMaterialPoint &mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J);
	
	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

public:
	int Probes() { return (int) m_probe.size(); }
	FEMicroProbe& Probe(int i) { return *m_probe[i]; }

protected:
	FEVecPropertyT<FEMicroProbe>	m_probe;

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
