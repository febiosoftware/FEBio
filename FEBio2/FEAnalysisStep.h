#pragma once

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
	FEAnalysisStep(FEModel& fem, int ntype = 0);

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

protected:
	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

public:
	// --- Control Data ---
	bool	m_baugment;		//!< use Lagrangian augmentation (TODO: move to solver class?)
};

//-----------------------------------------------------------------------------
//! Analysis class for structural mechanics problems
class FESolidAnalysis : public FEAnalysisStep
{
public:
	FESolidAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_SOLID) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for biphasic problems
class FEBiphasicAnalysis : public FEAnalysisStep
{
public:
	FEBiphasicAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_BIPHASIC) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for biphasic-solute problems
class FEBiphasicSoluteAnalysis : public FEAnalysisStep
{
public:
	FEBiphasicSoluteAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_POROSOLUTE) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for triphasic problems
class FETriphasicAnalysis : public FEAnalysisStep
{
public:
	FETriphasicAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_TRIPHASIC) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for triphasic problems
class FEHeatTransferAnalysis : public FEAnalysisStep
{
public:
	FEHeatTransferAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_HEAT) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for linear elastic problems
class FELinearSolidAnalysis : public FEAnalysisStep
{
public:
	FELinearSolidAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_LINEAR_SOLID) {}
};

//-----------------------------------------------------------------------------
//! Analysis class for linear thermo-elastic problems
class FEThermoElasticAnalysis : public FEAnalysisStep
{
public:
	FEThermoElasticAnalysis(FEModel& fem) : FEAnalysisStep(fem, FE_HEAT_SOLID) {}
};
