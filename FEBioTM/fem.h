#pragma once
#include "FECore/FEModel.h"

class FEM : public FEModel
{
public:
	// Initialization
	virtual bool Init() { return false; }

	//! Resets data structures
	virtual bool Reset() { return false; }

	// solve the model
	virtual bool Solve() { return false; }

	virtual void PushState() {}
	virtual void PopState () {}

	//! Evaluate parameter list
	virtual void EvaluateParameterList(FEParameterList& pl) {}
	virtual void EvaluateMaterialParameters(FEMaterial* pm) {}

	//! return a pointer to the named variable
	virtual double* FindParameter(const char* szname) { return 0; }

	//! find a boundary condition from the ID
	virtual FEBoundaryCondition* FindBC(int nid) { return 0; }

	//! check for user interruption
	virtual void CheckInterruption() {}

	// input from file
	virtual bool Input(const char* szfile) { return false; }

	//! write to plot file
	virtual void Write() {}

	//! write data to log file
	virtual void WriteData() {}

	//! write data to dump file
	virtual void DumpData() {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return false; }

	//! Add data record
	virtual void AddDataRecord(DataRecord* pd) {}
};
