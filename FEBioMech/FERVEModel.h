#pragma once
#include "FECore/FEModel.h"

class FEBCPrescribedDeformation;

//-----------------------------------------------------------------------------
// Class describing the RVE model.
// This is used by the homogenization code.
class FERVEModel : public FEModel
{
public:
	FERVEModel();
	~FERVEModel();

	//! one time initialization
	bool InitRVE(bool bperiodic, const char* szbc);

public:
	//! Return the initial volume (calculated in Init)
	double InitialVolume() const { return m_V0; }

	//! see if node is boundary node
	bool IsBoundaryNode(int i) const { return (m_BN[i]==1); }

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
	double			m_V0;			//!< initial volume
	bool			m_bperiodic;	//!< periodic BCs flag
	FEBoundingBox	m_bb;		//!< bounding box of mesh
	vector<int>		m_BN;			//!< boundary node flags
	FEBCPrescribedDeformation*	m_PD;	//!< prescribed deformation
};
