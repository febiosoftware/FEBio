#pragma once

#include "FECore/FEBoundaryCondition.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/FENLSolver.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! The FEAnalysisStep contains the data that describes a complete, single analysis.

class FEAnalysisStep : public FEAnalysis
{
public:
	//! constructor
	FEAnalysisStep(FEModel& fem);

	//! descructor
	virtual ~FEAnalysisStep(void);

	//! initialize step
	bool Init();

	//! Solve the analysis step
	bool Solve();

	//! wrap it up
	void Finish();

	//! Serialize data from and to a binary archive
	void Serialize(DumpFile& ar);

	//! add a boundary condition to the analysis
	void AddBoundaryCondition(FEBoundaryCondition* pbc) { m_BC.push_back(pbc); }

protected:
	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

public:
	// --- Control Data ---
	//{
		bool	m_baugment;		//!< use Lagrangian augmentation
	//}

	// --- Boundary conditions data ---
	//{
		vector<FEBoundaryCondition*>	m_BC;	//!< array of boundary conditions
	//}
};
