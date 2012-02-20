#pragma once
#include "FEBioLib/FEBioModel.h"
#include "FEBioPlot/PlotFile.h"
#include "FEBioLib/Timer.h"
#include "FEBioLib/DataStore.h"
#include <FL/Fl_Progress.H>
#include <FL/Fl.H>

//-----------------------------------------------------------------------------
class FETMProgress : public Progress
{
public:
	FETMProgress(Fl_Progress* pw) : m_pw(pw) { pw->maximum(100.f); pw->minimum(0.f); pw->value(0.f); }
	void SetProgress(double f) { m_pw->value((float) f); Fl::flush(); }
protected:
	Fl_Progress*	m_pw;
};

//-----------------------------------------------------------------------------
// The FE model class
class FEM : public FEBioModel
{
public:
	// constructor
	FEM();

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

	//! write to plot file
	virtual void Write() {}

	//! write data to log file
	virtual void WriteData() {}

	//! write data to dump file
	virtual void DumpData() {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return false; }
};
