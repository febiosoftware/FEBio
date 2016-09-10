#pragma once
#include "FECore/FEModel.h"

class FEBCPrescribedDeformation;

//-----------------------------------------------------------------------------
// Class describing the RVE model.
// This is used by the homogenization code.
class FERVEModel : public FEModel
{
public:
	enum RVE_TYPE
	{
		DISPLACEMENT,		// prescribed displacement
		PERIODIC_AL,		// periodic, augmented Lagrangian
		PERIODIC_LC			// periodic, linear constraints
	};

public:
	FERVEModel();
	~FERVEModel();

	//! one time initialization
	bool InitRVE(int rveType, const char* szbc, const char* szforce = 0);

	//! Return the initial volume (calculated in Init)
	double InitialVolume() const { return m_V0; }

	//! see if node is boundary node
	bool IsBoundaryNode(int i) const { return (m_BN[i]==1); }

	//! Update the RVE (before it is solved)
	void Update(const mat3d& F);

	// copy from the master RVE
	void CopyFrom(FERVEModel& rve);

protected:
	//! Calculate the initial volume
	void EvalInitialVolume();

	//! find the list of boundary nodes
	void FindBoundaryNodes(vector<int>& BN);

	//! Center the RVE
	void CenterRVE();

	bool PrepDisplacementBC();
	bool PrepPeriodicBC(const char* szbc);

private:
	double			m_V0;				//!< initial volume
	int				m_bctype;			//!< RVE type
	FEBoundingBox	m_bb;				//!< bounding box of mesh
	vector<int>		m_BN;				//!< boundary node flags

public:
	vector<int>		m_FN;				//!< list of force nodes
};
