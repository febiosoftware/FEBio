/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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



#include "stdafx.h"
#include "FEBioModel.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioXML/FERestartImport.h"
#include <FECore/NodeDataRecord.h>
#include <FECore/FaceDataRecord.h>
#include <FECore/ElementDataRecord.h>
#include <FEBioMech/ObjectDataRecord.h>
#include <FECore/NLConstraintDataRecord.h>
#include <FEBioMech/FERigidConnector.h>
#include <FEBioMech/FEGenericRigidJoint.h>
#include <FEBioMech/FERigidSphericalJoint.h>
#include <FEBioMech/FERigidPrismaticJoint.h>
#include <FEBioMech/FERigidRevoluteJoint.h>
#include <FEBioMech/FERigidCylindricalJoint.h>
#include <FEBioMech/FERigidPlanarJoint.h>
#include <FEBioMech/FERigidDamper.h>
#include <FEBioMech/FERigidSpring.h>
#include <FEBioMech/FERigidAngularDamper.h>
#include <FEBioMech/FERigidContractileForce.h>
#include "FEBioModelBuilder.h"
#include "FECore/log.h"
#include "FECore/FECoreKernel.h"
#include "FECore/DumpFile.h"
#include "FECore/DOFS.h"
#include <FECore/FEAnalysis.h>
#include <NumCore/MatrixTools.h>
#include <FECore/LinearSolver.h>
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/FETimeStepController.h>
#include "febio.h"
#include "version.h"
#include <iostream>
#include <sstream>
#include <fstream>

#ifdef WIN32
size_t FEBIOLIB_API GetPeakMemory();	// in memory.cpp
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEBioModel, FEMechModel)
	ADD_PARAMETER(m_title   , "title"    );
	ADD_PARAMETER(m_logLevel, "log_level");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// echo the input data to the log file
extern void echo_input(FEBioModel& fem);

bool FEBioModel::handleCB(FEModel* fem, int unsigned nwhen, void* pd)
{
	FEBioModel* febioModel = (FEBioModel*)pd;
	return febioModel->processEvent(nwhen);
}

//-----------------------------------------------------------------------------
bool FEBioModel::processEvent(int nevent)
{
	// write output files (but not while serializing)
	if ((nevent == CB_SERIALIZE_LOAD) || (nevent == CB_SERIALIZE_SAVE)) return true;
	Write(nevent);

	// process event handlers
	switch (nevent)
	{
	case CB_STEP_SOLVED: on_cb_stepSolved(); break;
	case CB_SOLVED     : on_cb_solved(); break;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Constructor of FEBioModel class.
FEBioModel::FEBioModel()
{
	m_logLevel = 1;

	m_dumpLevel = FE_DUMP_NEVER;
	m_dumpStride = 1;

	// --- I/O-Data ---
	m_ndebug = 0;
	m_becho = true;
	m_plot = nullptr;
	m_writeMesh = false;

	m_stats.ntimeSteps = 0;
	m_stats.ntotalIters = 0;
	m_stats.ntotalRHS = 0;
	m_stats.ntotalReforms = 0;

	m_pltAppendOnRestart = true;

	m_lastUpdate = -1;

	m_bshowErrors = true;

	// Add the output callback
	// We call this function always since we want to flush the logfile for each event.
	AddCallback(handleCB, CB_ALWAYS, this);
}

//-----------------------------------------------------------------------------
FEBioModel::~FEBioModel()
{
	// close the plot file
	if (m_plot) { delete m_plot; m_plot = 0; }
	m_log.close();
}

//-----------------------------------------------------------------------------
void FEBioModel::ShowWarningsAndErrors(bool b)
{
	m_bshowErrors = b;
}

//-----------------------------------------------------------------------------
bool FEBioModel::ShowWarningsAndErrors() const
{
	return m_bshowErrors;
}

//-----------------------------------------------------------------------------
Timer& FEBioModel::GetSolveTimer()
{
	return *GetTimer(Timer_ModelSolve);
}

//-----------------------------------------------------------------------------
//! return number of seconds of time spent in linear solver
int FEBioModel::GetLinearSolverTime()
{
	Timer* t = GetTimer(TimerID::Timer_LinSolve);
	return (int)t->peek();
}

//-----------------------------------------------------------------------------
//! set the debug level
void FEBioModel::SetDebugLevel(int debugLvl) { m_ndebug = debugLvl; }

//! get the debug level
int FEBioModel::GetDebugLevel() { return m_ndebug; }

//! set the dump level (for cold restarts)
void FEBioModel::SetDumpLevel(int dumpLevel) { m_dumpLevel = dumpLevel; }

//! get the dump level
int FEBioModel::GetDumpLevel() const { return m_dumpLevel; }

//! Set the dump stride
void FEBioModel::SetDumpStride(int n) { m_dumpStride = n; }

//! get the dump stride
int FEBioModel::GetDumpStride() const { return m_dumpStride; }

//! Set the log level
void FEBioModel::SetLogLevel(int logLevel) { m_logLevel = logLevel; }

//! Get the stats 
ModelStats FEBioModel::GetModelStats() const
{
	return m_stats;
}

//-----------------------------------------------------------------------------
//! Set the title of the model
void FEBioModel::SetTitle(const char* sz)
{
	m_title = sz;
}

//-----------------------------------------------------------------------------
//! Return the title of the model
const std::string& FEBioModel::GetTitle() const
{
	return m_title;
}

//-----------------------------------------------------------------------------
double FEBioModel::GetEndTime() const
{
	return GetCurrentStep()->m_tend;
}

//=============================================================================
//
//		FEBioModel: I-O Functions
//
//=============================================================================

//-----------------------------------------------------------------------------
//! Add a data record to the data store
void FEBioModel::AddDataRecord(DataRecord* pd)
{
	DataStore& dataStore = GetDataStore();
	dataStore.AddRecord(pd); 
}

//-----------------------------------------------------------------------------
//! Get the plot file
PlotFile* FEBioModel::GetPlotFile()
{
	return m_plot;
}

//-----------------------------------------------------------------------------
//! Sets the name of the FEBio input file
void FEBioModel::SetInputFilename(const std::string& sfile)
{ 
	m_sfile = sfile;
	size_t npos = sfile.rfind('/');
	if (npos != string::npos) m_sfile_title = sfile.substr(npos+1, std::string::npos);
	if (m_sfile_title.empty()) 
	{
		npos = sfile.rfind('\\');
		if (npos != string::npos) m_sfile_title = sfile.substr(npos, string::npos);
		if (m_sfile_title.empty()) m_sfile_title = m_sfile;
	}
}

//-----------------------------------------------------------------------------
//! Set the name of the log file
void FEBioModel::SetLogFilename(const std::string& sfile)
{ 
	m_slog = sfile; 
}

//-----------------------------------------------------------------------------
//! Set the name of the plot file
void FEBioModel::SetPlotFilename(const std::string& sfile)
{ 
	m_splot = sfile;
}

//-----------------------------------------------------------------------------
//! Set the name of the restart archive (i.e. the dump file)
void FEBioModel::SetDumpFilename(const std::string& sfile)
{ 
	m_sdump = sfile;
}

//-----------------------------------------------------------------------------
//! Return the name of the input file
const std::string& FEBioModel::GetInputFileName()
{ 
	return m_sfile; 
}

//-----------------------------------------------------------------------------
//! Return the name of the log file
const std::string& FEBioModel::GetLogfileName()
{ 
	return m_slog;  
}

//-----------------------------------------------------------------------------
//! Return the name of the plot file
const std::string& FEBioModel::GetPlotFileName()
{ 
	return m_splot; 
}

//-----------------------------------------------------------------------------
//! Return the dump file name.
const std::string& FEBioModel::GetDumpFileName()
{
	return	m_sdump;
}

//-----------------------------------------------------------------------------
//! get the file title (i.e. name of input file without the path)
const std::string& FEBioModel::GetFileTitle()
{ 
	return m_sfile_title; 
}

//-----------------------------------------------------------------------------
// set append-on-restart flag
void FEBioModel::SetAppendOnRestart(bool b)
{
	m_pltAppendOnRestart = b;
}

//-----------------------------------------------------------------------------
bool FEBioModel::AppendOnRestart() const
{
	return m_pltAppendOnRestart;
}

//=============================================================================
//    I N P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FEBioModel::Input(const char* szfile)
{
	// start the timer
	TimerTracker t(&m_InputTime);

	// create file reader
	FEBioImport fim;

	// override the default model builder
	fim.SetModelBuilder(new FEBioModelBuilder(*this));

	feLog("Reading file %s ...", szfile);

	// Load the file
	if (fim.Load(*this, szfile) == false)
	{
		feLog("FAILED!\n");
		char szerr[256];
		fim.GetErrorMessage(szerr);
		feLogError(szerr);

		return false;
	}
	else feLog("SUCCESS!\n");

	// set the input file name
	SetInputFilename(szfile);

	// see if user redefined output filenames
	if (fim.m_szdmp[0]) SetDumpFilename(fim.m_szdmp);
	if (fim.m_szlog[0]) SetLogFilename (fim.m_szlog);
	if (fim.m_szplt[0]) SetPlotFilename(fim.m_szplt);

	// add the data records
	int ND = (int)fim.m_data.size();
	for (int i=0; i<ND; ++i) AddDataRecord(fim.m_data[i]);

	// we're done reading
	return true;
}

//-----------------------------------------------------------------------------
//! This function finds all the domains that have a certain material
void FEBioModel::DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom)
{
	FEMesh& mesh = GetMesh();

	// make sure the list is empty
	if (ldom.empty() == false) ldom.clear();

	// loop over all domains
	int ND = mesh.Domains();
	int NM = (int)lmat.size();
	for (int i = 0; i<ND; ++i)
	{
		FEDomain& di = mesh.Domain(i);
		int dmat = di.GetMaterial()->GetID();
		for (int j = 0; j<NM; ++j)
		{
			if (dmat == lmat[j])
			{
				ldom.push_back(i);
				break;
			}
		}
	}
}

//=============================================================================
//    O U T P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! Export state to plot file.
void FEBioModel::Write(unsigned int nevent)
{
	TimerTracker t(&m_IOTimer);

	// get the current step
	FEAnalysis* pstep = GetCurrentStep();

	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	if ((nevent == CB_INIT) && m_becho)
	{
		Logfile::MODE old_mode = m_log.GetMode();

		// don't output when no output is requested
		if (old_mode != Logfile::LOG_NEVER)
		{
			// we only output this data to the log file and not the screen
			m_log.SetMode(Logfile::LOG_FILE);

			// write output
			echo_input();

			// reset log mode
			m_log.SetMode(old_mode);
		}
	}

	// update plot file
	WritePlot(nevent);

	// Dump converged state to the archive
	DumpData(nevent);

	// write the output data
	WriteData(nevent);
}

//-----------------------------------------------------------------------------
void FEBioModel::WritePlot(unsigned int nevent)
{
	// get the current step
	FEAnalysis* pstep = GetCurrentStep();

	// get the plot level
	int nplt = pstep->GetPlotLevel();

	// if we don't want to plot anything we return
	if (nplt != FE_PLOT_NEVER)
	{
		// try to open the plot file
		if ((nevent == CB_INIT) || (nevent == CB_STEP_ACTIVE))
		{
			// If the first step did not request output, m_plot can still be null
			if (m_plot == 0) InitPlotFile();

			if (m_plot->IsValid() == false)
			{
				// Add the plot objects
				UpdatePlotObjects();

				if (m_plot->Open(m_splot.c_str()) == false)
				{
					feLog("ERROR : Failed creating PLOT database\n");
					delete m_plot;
					m_plot = 0;
				}

				// Since it is assumed that for the first timestep
				// there are no loads or initial displacements, the case n=0 is skipped.
				// Therefor we can output those results here.
				// TODO: Offcourse we should actually check if this is indeed
				//       the case, otherwise we should also solve for t=0
				// Only output the initial state if requested
				if (m_plot)
				{
					bool bout = true;

					// if we're using the fixed time stepper, we check the plot range and zero state flag
					if (pstep->m_timeController == nullptr) bout = (pstep->m_nplotRange[0] == 0) || (pstep->m_bplotZero);

					// for multi-step analyses, some plot variables can't be plot until the first step
					// is activated. So, in that case, we'll wait. 
					if ((nevent == CB_INIT) && (Steps() > 1)) bout = false;

					// store initial time step (i.e. time step zero)
					if (bout)
					{
						double time = GetTime().currentTime;
						m_plot->Write((float)time);
					}
				}
			}
			else
			{
				// for multi-step analyses, we did not write the initial time step during CB_INIT
				// so we'll do it during the activation of the first step. 
				if ((nevent == CB_STEP_ACTIVE) && (Steps() > 1) && (GetCurrentStepIndex() == 0))
				{
					double time = GetTime().currentTime;
					m_plot->Write((float)time);
				}
			}
		}
		else
		{
			// assume we won't be writing anything
			bool bout = false;

			// see if we need to output something
			int ndebug = GetDebugLevel();

			// write a new mesh section if needed
			if (nevent == CB_REMESH)
			{
				m_writeMesh = true;
				m_lastUpdate = -1;
			}

			if (ndebug == 1)
			{
				if ((nevent == CB_INIT) || (nevent == CB_MODEL_UPDATE) || (nevent == CB_MINOR_ITERS) || (nevent == CB_SOLVED) || (nevent == CB_REMESH))
				{
					bout = true;
				}

				if (nevent == CB_MAJOR_ITERS)
				{
					bout = true;
					m_lastUpdate = -1;
				}
			}
			else
			{
				int currentStep = pstep->m_ntimesteps;
				int lastStep = pstep->m_ntime;
				int nmin = pstep->m_nplotRange[0]; if (nmin < 0) nmin = lastStep + nmin + 1;
				int nmax = pstep->m_nplotRange[1]; if (nmax < -1) nmax = lastStep + nmax + 1;

				bool inRange = true;
				bool isStride = true;
				if (pstep->m_timeController == nullptr)
				{
					inRange = false;
					if ((currentStep >= nmin) && ((currentStep <= nmax) || (nmax == -1))) inRange = true;

				}
				isStride = ((pstep->m_ntimesteps - nmin) % pstep->m_nplot_stride) == 0;

				switch (nevent)
				{
				case CB_MINOR_ITERS:
				{
					if (nplt == FE_PLOT_MINOR_ITRS) bout = true;
					if ((ndebug == 2) && (NegativeJacobian::IsThrown()))
					{
						bout = true;
						NegativeJacobian::clearFlag();
					}
				}
				break;
				case CB_MAJOR_ITERS:
					if ((nplt == FE_PLOT_MAJOR_ITRS) && inRange && isStride) bout = true;
					if ((nplt == FE_PLOT_MUST_POINTS) && (pstep->m_timeController) && (pstep->m_timeController->m_nmust >= 0)) bout = true;
					if (nplt == FE_PLOT_AUGMENTATIONS) bout = true;
					break;
				case CB_AUGMENT:
					if (nplt == FE_PLOT_AUGMENTATIONS)
					{
						// Note that this is called before the augmentations.
						// The reason we store the state prior to the augmentations
						// is because the augmentations are going to change things such that
						// the system no longer in equilibrium. Since the model has to be converged
						// before we do augmentations, storing the model now will store an actual converged state.
						bout = true;
					}
					break;
				case CB_SOLVED:
					if (nplt == FE_PLOT_FINAL) bout = true;
					if (nplt == FE_PLOT_MAJOR_ITRS)
					{
						// we want to force storing the final time step, but we have to make sure 
						// it hasn't been stored already during the CB_MAJOR_ITERS callback.
						if ((inRange == false) || (isStride == false)) bout = true;
					}
					break;
				case CB_STEP_SOLVED: if (nplt == FE_PLOT_STEP_FINAL) bout = true;  break;
				case CB_USER1: if ((nplt == FE_PLOT_USER1) && inRange && isStride) bout = true; break;
				}
			}

			// output the state if requested
			if (bout && (m_lastUpdate != UpdateCounter()))
			{
				m_lastUpdate = UpdateCounter();

				// update the plot objects
				UpdatePlotObjects();

				// see if we need to write a new mesh section
				if (m_writeMesh) {
					FEBioPlotFile* plt = dynamic_cast<FEBioPlotFile*>(m_plot);
					feLogDebug("writing mesh section to plot file");
					plt->WriteMeshSection(*this);
				}

				// set the status flag
				int statusFlag = 0;
				if (m_writeMesh) statusFlag = 1;
				else if (nevent != CB_MAJOR_ITERS)
				{
					statusFlag = 2;
				}

				// write the state section
				double time = GetTime().currentTime;
				if (m_plot)
				{
					feLogDebug("writing to plot file; time = %lg; flag = %d", time, statusFlag);
					m_plot->Write((float)time, statusFlag);
				}

				// make sure to reset write mesh flag
				m_writeMesh = false;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Write user data to the logfile
void FEBioModel::WriteData(unsigned int nevent)
{
	// get the current step
	FEAnalysis* pstep = GetCurrentStep();
	int nout = pstep->GetOutputLevel();
	if (nout == FE_OUTPUT_NEVER) return;

	// see if we need to output
	bool bout = false;
	switch (nevent)
	{
	case CB_INIT: if (nout == FE_OUTPUT_MAJOR_ITRS) bout = true; break;
	case CB_MINOR_ITERS: if (nout == FE_OUTPUT_MINOR_ITRS) bout = true; break;
	case CB_MAJOR_ITERS:
		if (nout == FE_OUTPUT_MAJOR_ITRS) bout = true;
		if ((nout == FE_OUTPUT_MUST_POINTS) && (pstep->m_timeController) && (pstep->m_timeController->m_nmust >= 0)) bout = true;
		break;
	case CB_SOLVED:
		if (nout == FE_OUTPUT_FINAL) bout = true;
		break;
	}

	// output data
	if (bout)
	{
		DataStore& dataStore = GetDataStore();
		dataStore.Write();
	}
}

//-----------------------------------------------------------------------------
//! Dump state to archive for restarts
void FEBioModel::DumpData(int nevent)
{
	// get the current step
	FEAnalysis* pstep = GetCurrentStep();
	int ndump = GetDumpLevel();
	int stride = GetDumpStride();
	if (ndump == FE_DUMP_NEVER) return;

	bool bdump = false;
	switch (nevent)
	{
	case CB_MAJOR_ITERS:
		if (ndump == FE_DUMP_MAJOR_ITRS)
		{
			if (stride <= 1) bdump = true;
			else
			{
				int niter = pstep->m_ntimesteps;
				bdump = ((niter % stride) == 0);
			}
		}
		if ((ndump == FE_DUMP_MUST_POINTS) && (pstep->m_timeController) && (pstep->m_timeController->m_nmust >= 0)) bdump = true;
		break;
	case CB_STEP_SOLVED: if (ndump == FE_DUMP_STEP) bdump = true; break;
	}
	
	if (bdump)
	{
		DumpFile ar(*this);
		if (ar.Create(m_sdump.c_str()) == false)
		{
			feLogWarning("Failed creating restart file (%s).\n", m_sdump.c_str());
		}
		else
		{
			Serialize(ar);
			feLogInfo("\nRestart point created. Archive name is %s.", m_sdump.c_str());
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::Log(int ntag, const char* szmsg)
{
	if      (ntag == 0) m_log.printf(szmsg);
	else if ((ntag == 1) && m_bshowErrors) m_log.printbox("WARNING", szmsg);
	else if ((ntag == 2) && m_bshowErrors) m_log.printbox("ERROR", szmsg);
	else if (ntag == 3) m_log.printbox(nullptr, szmsg);
	else if (ntag == 4)
	{
		if (GetDebugLevel() > 0)
			m_log.printf("debug>%s\n", szmsg);
	}

	// Flushing the logfile each time we get here might be a bit overkill.
	// For now, I'm flushing the log file in the output_cb method.
//	m_log.flush();
}

//-----------------------------------------------------------------------------

class FEPlotRigidBodyPosition : public FEPlotObjectData
{
public:
	FEPlotRigidBodyPosition(FEModel* fem, FERigidBody* prb) : FEPlotObjectData(fem), m_rb(prb) {}

	bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
	{
		assert(m_rb);
		ar << m_rb->m_rt;
		return true;
	}

private:
	FERigidBody* m_rb;
};

class FEPlotRigidBodyRotation : public FEPlotObjectData
{
public:
	FEPlotRigidBodyRotation(FEModel* fem, FERigidBody* prb) : FEPlotObjectData(fem), m_rb(prb) {}

	bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
	{
		assert(m_rb);
		quatd q = m_rb->GetRotation();
		vec3d e;
		q.GetEuler(e.x, e.y, e.z);
		e *= RAD2DEG;
		ar << e;
		return true;
	}

private:
	FERigidBody* m_rb;
};

class FEPlotRigidBodyForce : public FEPlotObjectData
{
public:
	FEPlotRigidBodyForce(FEModel* fem, FERigidBody* prb) : FEPlotObjectData(fem), m_rb(prb) {}

	bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
	{
		assert(m_rb);
		ar << m_rb->m_Fr;
		return true;
	}

private:
	FERigidBody* m_rb;
};

class FEPlotRigidBodyMoment : public FEPlotObjectData
{
public:
	FEPlotRigidBodyMoment(FEModel* fem, FERigidBody* prb) : FEPlotObjectData(fem), m_rb(prb) {}

	bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
	{
		assert(m_rb);
		ar << m_rb->m_Mr;
		return true;
	}

private:
	FERigidBody* m_rb;
};


//-----------------------------------------------------------------------------
class FEPlotRigidConnectorTranslationLCS : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorTranslationLCS(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->RelativeTranslation(false);
        return true;
    }

private:
    FERigidConnector* m_rc;
};

class FEPlotRigidConnectorRotationLCS : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorRotationLCS(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->RelativeRotation(false);
        return true;
    }

private:
    FERigidConnector* m_rc;
};

class FEPlotRigidConnectorTranslationGCS : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorTranslationGCS(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->RelativeTranslation(true);
        return true;
    }

private:
    FERigidConnector* m_rc;
};

class FEPlotRigidConnectorRotationGCS : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorRotationGCS(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->RelativeRotation(true);
        return true;
    }

private:
    FERigidConnector* m_rc;
};

class FEPlotRigidConnectorForce : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorForce(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->m_F;
        return true;
    }

private:
    FERigidConnector* m_rc;
};

class FEPlotRigidConnectorMoment : public FEPlotObjectData
{
public:
    FEPlotRigidConnectorMoment(FEModel* fem, FERigidConnector* prb) : FEPlotObjectData(fem), m_rc(prb) {}

    bool Save(FEBioPlotFile::PlotObject* po, FEDataStream& ar)
    {
        assert(m_rc);
        ar << m_rc->m_M;
        return true;
    }

private:
    FERigidConnector* m_rc;
};

//-----------------------------------------------------------------------------
void FEBioModel::UpdatePlotObjects()
{
	FEBioPlotFile* plt = dynamic_cast<FEBioPlotFile*>(m_plot);
	if (plt == nullptr) return;

	int nrb = RigidBodies();
	if (nrb == 0) return;

	FEModel& fem = *GetFEModel();

	if (plt->PointObjects() == 0)
	{
		int nid = 1;
		for (int i = 0; i < nrb; ++i)
		{
			FERigidBody* rb = GetRigidBody(i);
			string name = rb->GetName();
			if (name.empty())
			{
				stringstream ss;
				ss << "Object" << nid;
				name = ss.str();
			}

			FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
			po->m_tag = 1;
			po->m_pos = rb->m_r0;
			po->m_rot = quatd(0, vec3d(1,0,0));

			po->AddData("Position", PLT_VEC3F, new FEPlotRigidBodyPosition(this, rb));
			po->AddData("Euler angles", PLT_VEC3F, new FEPlotRigidBodyRotation(this, rb));
			po->AddData("Force" , PLT_VEC3F, new FEPlotRigidBodyForce(this, rb));
			po->AddData("Moment", PLT_VEC3F, new FEPlotRigidBodyMoment(this, rb));

			nid++;
		}

		// check rigid connectors
		for (int i = 0; i < fem.NonlinearConstraints(); ++i)
		{
			FENLConstraint* pc = fem.NonlinearConstraint(i);

			string name = pc->GetName();
			if (name.empty())
			{
				stringstream ss;
				ss << "Object" << nid;
				name = ss.str();
			}

			FEGenericRigidJoint* rj = dynamic_cast<FEGenericRigidJoint*>(pc);
			if (rj)
			{
				FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
				po->m_tag = 2;
				po->m_pos = rj->InitialPosition();
				po->m_rot = quatd(0, vec3d(1, 0, 0));
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rj));
			}

            FERigidSphericalJoint* rsj = dynamic_cast<FERigidSphericalJoint*>(pc);
            if (rsj)
            {
                FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
                po->m_tag = 3;
                po->m_pos = rsj->InitialPosition();
                po->m_rot = quatd(0, vec3d(1, 0, 0));
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rsj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rsj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rsj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rsj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rsj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rsj));
            }

			FERigidPrismaticJoint* rpj = dynamic_cast<FERigidPrismaticJoint*>(pc);
			if (rpj)
			{
				FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
				po->m_tag = 4;
				po->m_pos = rpj->InitialPosition();
				po->m_rot = rpj->Orientation();
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rpj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rpj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rpj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rpj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rpj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rpj));
			}

			FERigidRevoluteJoint* rrj = dynamic_cast<FERigidRevoluteJoint*>(pc);
			if (rrj)
			{
				FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
				po->m_tag = 5;
				po->m_pos = rrj->InitialPosition();
				po->m_rot = rrj->Orientation();
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rrj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rrj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rrj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rrj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rrj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rrj));
			}

			FERigidCylindricalJoint* rcj = dynamic_cast<FERigidCylindricalJoint*>(pc);
			if (rcj)
			{
				FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
				po->m_tag = 6;
				po->m_pos = rcj->InitialPosition();
				po->m_rot = rcj->Orientation();
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rcj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rcj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rcj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rcj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rcj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rcj));
			}
            
            FERigidPlanarJoint* rlj = dynamic_cast<FERigidPlanarJoint*>(pc);
            if (rlj)
            {
                FEBioPlotFile::PointObject* po = plt->AddPointObject(name);
                po->m_tag = 7;
                po->m_pos = rlj->InitialPosition();
                po->m_rot = rlj->Orientation();
                po->AddData("Relative translation (LCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationLCS(this, rlj));
                po->AddData("Relative rotation (LCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationLCS(this, rlj));
                po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rlj));
                po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rlj));
                po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rlj));
                po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rlj));
            }
			nid++;
		}
	}
	else
	{
		for (int i = 0; i < nrb; ++i)
		{
			FERigidBody* rb = GetRigidBody(i);
			FEBioPlotFile::PointObject* po = plt->GetPointObject(i);
			po->m_pos = rb->m_rt;
			po->m_rot = rb->GetRotation();
		}

		// check rigid connectors
		int n = nrb;
		for (int i = 0; i < fem.NonlinearConstraints(); ++i)
		{
			FENLConstraint* pc = fem.NonlinearConstraint(i);

			FEGenericRigidJoint* rj = dynamic_cast<FEGenericRigidJoint*>(pc);
			if (rj)
			{
				FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
				po->m_pos = rj->Position();
				po->m_rot = quatd(0, vec3d(1, 0, 0));
			}

            FERigidSphericalJoint* rsj = dynamic_cast<FERigidSphericalJoint*>(pc);
            if (rsj)
            {
                FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
                po->m_pos = rsj->Position();
                po->m_rot = quatd(0, vec3d(1, 0, 0));
            }

			FERigidPrismaticJoint* rpj = dynamic_cast<FERigidPrismaticJoint*>(pc);
			if (rpj)
			{
				FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
				po->m_pos = rpj->Position();
				po->m_rot = rpj->Orientation();
			}

			FERigidRevoluteJoint* rrj = dynamic_cast<FERigidRevoluteJoint*>(pc);
			if (rrj)
			{
				FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
				po->m_pos = rrj->Position();
				po->m_rot = rrj->Orientation();
			}

			FERigidCylindricalJoint* rcj = dynamic_cast<FERigidCylindricalJoint*>(pc);
			if (rcj)
			{
				FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
				po->m_pos = rcj->Position();
				po->m_rot = rcj->Orientation();
			}
            
            FERigidPlanarJoint* rlj = dynamic_cast<FERigidPlanarJoint*>(pc);
            if (rlj)
            {
                FEBioPlotFile::PointObject* po = plt->GetPointObject(n++);
                po->m_pos = rlj->Position();
                po->m_rot = rlj->Orientation();
            }
		}
	}

	if (plt->LineObjects() == 0)
	{
		int nid = 1;

		// check rigid connectors
		for (int i = 0; i < fem.NonlinearConstraints(); ++i)
		{
			FERigidConnector* prc = dynamic_cast<FERigidConnector*>(fem.NonlinearConstraint(i));
			if (prc)
			{
				string name = prc->GetName();
				if (name.empty())
				{
					stringstream ss;
					ss << "Object" << nid;
					name = ss.str();
				}

				vec3d ra = GetRigidBody(prc->m_nRBa)->m_r0;
				vec3d rb = GetRigidBody(prc->m_nRBb)->m_r0;

				FERigidSpring* rs = dynamic_cast<FERigidSpring*>(prc);
				if (rs)
				{
					FEBioPlotFile::LineObject* po = plt->AddLineObject(name);
					po->m_tag = 1;
					po->m_rot = quatd(0, vec3d(1, 0, 0));
					po->m_r1 = rs->m_at;
					po->m_r2 = rs->m_bt;
                    po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rs));
                    po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rs));
                    po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rs));
                    po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rs));
				}

				FERigidDamper* rd = dynamic_cast<FERigidDamper*>(prc);
				if (rd)
				{
					FEBioPlotFile::LineObject* po = plt->AddLineObject(name);
					po->m_tag = 2;
                    po->m_rot = quatd(0, vec3d(1, 0, 0));
					po->m_r1 = rd->m_at;
					po->m_r2 = rd->m_bt;
                    po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rd));
                    po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rd));
                    po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rd));
                    po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rd));
				}
                
                FERigidAngularDamper* rad = dynamic_cast<FERigidAngularDamper*>(prc);
                if (rad)
                {
                    FEBioPlotFile::LineObject* po = plt->AddLineObject(name);
                    po->m_tag = 3;
                    po->m_r1 = ra;
                    po->m_r2 = rb;
                    po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rad));
                    po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rad));
                    po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rad));
                    po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rad));
                }
                
                FERigidContractileForce* rcf = dynamic_cast<FERigidContractileForce*>(prc);
                if (rcf)
                {
                    FEBioPlotFile::LineObject* po = plt->AddLineObject(name);
                    po->m_tag = 4;
                    po->m_rot = quatd(0, vec3d(1, 0, 0));
                    po->m_r1 = rcf->m_at;
                    po->m_r2 = rcf->m_bt;
                    po->AddData("Relative translation (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorTranslationGCS(this, rcf));
                    po->AddData("Relative rotation (GCS)", PLT_VEC3F, new FEPlotRigidConnectorRotationGCS(this, rcf));
                    po->AddData("Reaction force (GCS)" , PLT_VEC3F, new FEPlotRigidConnectorForce(this, rcf));
                    po->AddData("Reaction moment (GCS)", PLT_VEC3F, new FEPlotRigidConnectorMoment(this, rcf));
                }
			}
		}
	}
	else
	{
		// check rigid connectors
		int n = 0;
		for (int i = 0; i < fem.NonlinearConstraints(); ++i)
		{
			FERigidConnector* prc = dynamic_cast<FERigidConnector*>(fem.NonlinearConstraint(i));
			if (prc)
			{
				vec3d ra = GetRigidBody(prc->m_nRBa)->m_rt;
				vec3d rb = GetRigidBody(prc->m_nRBb)->m_rt;

				FERigidSpring* rs = dynamic_cast<FERigidSpring*>(prc);
				if (rs)
				{
					FEBioPlotFile::LineObject* po = plt->GetLineObject(n++);
					po->m_r1 = rs->m_at;
					po->m_r2 = rs->m_bt;
				}

				FERigidDamper* rd = dynamic_cast<FERigidDamper*>(prc);
				if (rd)
				{
					FEBioPlotFile::LineObject* po = plt->GetLineObject(n++);
					po->m_r1 = ra;
					po->m_r2 = rb;
				}
                
                FERigidAngularDamper* rad = dynamic_cast<FERigidAngularDamper*>(prc);
                if (rad)
                {
                    FEBioPlotFile::LineObject* po = plt->GetLineObject(n++);
                    po->m_r1 = ra;
                    po->m_r2 = rb;
                }
                
                FERigidContractileForce* rcf = dynamic_cast<FERigidContractileForce*>(prc);
                if (rcf)
                {
                    FEBioPlotFile::LineObject* po = plt->GetLineObject(n++);
                    po->m_r1 = rcf->m_at;
                    po->m_r2 = rcf->m_bt;
                }
			}
		}
	}
}

//=============================================================================
//    R E S T A R T
//=============================================================================

class restart_exception : public std::runtime_error
{
public:
	restart_exception() : std::runtime_error("restart error") {}
	restart_exception(const char* msg) : std::runtime_error(msg) {}
};

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa DumpFile
void FEBioModel::Serialize(DumpStream& ar)
{
	// don't need to do anything for running restarts
	if (ar.IsShallow())
	{
		// serialize model data
		FEMechModel::Serialize(ar);
	}
	else
	{
		if (ar.IsSaving())
		{
			// --- version number ---
			ar << (int) RSTRTVERSION;
		}
		else
		{
			// --- version ---
			int nversion;
			ar >> nversion;

			// make sure it is the right version
			if (nversion != RSTRTVERSION) throw restart_exception("incorrect version number");
		}

		// serialize model data
		FEMechModel::Serialize(ar);

		// --- Save IO Data
		SerializeIOData(ar);
	}
}

//-----------------------------------------------------------------------------
//! Serialization of FEBioModel data
void FEBioModel::SerializeIOData(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		// file names
		ar << m_sfile << m_splot << m_slog << m_sdump;

		// plot file
		int npltfmt = 2;
		ar << npltfmt;

		SerializePlotData(ar);

		// data records
		SerializeDataStore(ar);
	}
	else
	{
		// file names
		string splot, slog, sdmp;
		ar >> m_sfile >> splot >> slog >> sdmp;

		// don't forget to call store the input file name so
		// that m_szfile_title gets initialized
		SetInputFilename(m_sfile);

		// If we append, use the original names
		// otherwise we use the names as was initialized by the command line parser
		if (m_pltAppendOnRestart)
		{
			m_splot = splot;
			m_slog = slog;
			m_sdump = sdmp;
		}

		// get the plot file format (should be 2)
		int npltfmt = 0;
		ar >> npltfmt;
		assert(npltfmt == 2);

		SerializePlotData(ar);

		// remove the plot file (if any)
		if (m_plot) { delete m_plot; m_plot = 0; }

		// create the plot file
		FEBioPlotFile* pplt = new FEBioPlotFile(this);
		m_plot = pplt;

		if (m_pltAppendOnRestart)
		{
			// Open for appending
			if (m_plot->Append(m_splot.c_str()) == false)
			{
				printf("FATAL ERROR: Failed reopening plot database %s\n", m_splot.c_str());
				throw "FATAL ERROR";
			}
		}
		else
		{
			// set the software string
			const char* szver = febio::getVersionString();
			char szbuf[256] = { 0 };
			sprintf(szbuf, "FEBio %s", szver);
			pplt->SetSoftwareString(szbuf);
		}

		// data records
		SerializeDataStore(ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::SerializePlotData(DumpStream& ar)
{
	GetPlotDataStore().Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEBioModel::SerializeDataStore(DumpStream& ar)
{
	DataStore& dataStore = GetDataStore();
	if (ar.IsSaving())
	{
		int N = dataStore.Size();
		ar << N;
		for (int i=0; i<N; ++i)
		{
			DataRecord* pd = dataStore.GetDataRecord(i);

			int ntype = pd->m_type;
			ar << ntype;
			pd->Serialize(ar);
		}
	}
	else
	{
		int N;
		dataStore.Clear();
		ar >> N;
		for (int i=0; i<N; ++i)
		{
			int ntype;
			ar >> ntype;

			DataRecord* pd = 0;
			switch(ntype)
			{
			case FE_DATA_NODE: pd = new NodeDataRecord        (this); break;
			case FE_DATA_FACE: pd = new FaceDataRecord        (this); break;
			case FE_DATA_ELEM: pd = new ElementDataRecord     (this); break;
			case FE_DATA_RB  : pd = new ObjectDataRecord      (this); break;
			case FE_DATA_NLC : pd = new NLConstraintDataRecord(this); break;
			}
			assert(pd);
			pd->Serialize(ar);
			dataStore.AddRecord(pd);
		}
	}
}

//=============================================================================
//    I N I T I A L I Z A T I O N
//=============================================================================

//-----------------------------------------------------------------------------
// Initialize plot file
bool FEBioModel::InitPlotFile()
{
	FEBioPlotFile* pplt = new FEBioPlotFile(this);
	m_plot = pplt;

	// set the software string
	const char* szver = febio::getVersionString();
	char szbuf[256] = { 0 };
	sprintf(szbuf, "FEBio %s", szver);
	pplt->SetSoftwareString(szbuf);
	
	// see if a valid plot file name is defined.
	const std::string& splt = GetPlotFileName();
	if (splt.empty())
	{
		// if not, we take the input file name and set the extension to .xplt
		char sz[1024] = { 0 };
		strcpy(sz, GetInputFileName().c_str());
		char* ch = strrchr(sz, '.');
		if (ch) *ch = 0;
		strcat(sz, ".xplt");
		SetPlotFilename(sz);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEBioModel::Init()
{
	TimerTracker t(&m_InitTime);

	// Open the logfile
	if (m_logLevel != 0)
	{
		if (InitLogFile() == false) return false;
	}

	FEBioPlotFile* pplt = nullptr;
	m_lastUpdate = -1;

	// open plot database file
	FEAnalysis* step = GetCurrentStep();
	if (step->GetPlotLevel() != FE_PLOT_NEVER)
	{
		if (m_plot == 0) InitPlotFile();
	}

	// see if a valid dump file name is defined.
	const std::string& sdmp = GetDumpFileName();
	if (sdmp.empty())
	{
		// if not, we take the input file name and set the extension to .dmp
		char sz[1024] = { 0 };
		strcpy(sz, GetInputFileName().c_str());
		char* ch = strrchr(sz, '.');
		if (ch) *ch = 0;
		strcat(sz, ".dmp");
		SetDumpFilename(sz);
	}

	// initialize data records
	DataStore& dataStore = GetDataStore();
	for (int i = 0; i < dataStore.Size(); ++i)
	{
		if (dataStore.GetDataRecord(i)->Initialize() == false) return false;
	}

	// initialize model data
	if (FEMechModel::Init() == false) 
	{
		feLogError("Model initialization failed");
		return false;
	}

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
// Opens the log file
bool FEBioModel::InitLogFile()
{
	// Only do this if the log file is not valid
	if (!m_log.is_valid())
	{
		// see if a valid log file name is defined.
		const std::string& slog = GetLogfileName();
		if (slog.empty())
		{
			// if not, we take the input file name and set the extension to .log
			char sz[1024] = {0};
			strcpy(sz, GetInputFileName().c_str());
			char *ch = strrchr(sz, '.');
			if (ch) *ch = 0;
			strcat(sz, ".log");
			SetLogFilename(sz);
		}

		// create a log stream
		LogFileStream* fp = new LogFileStream;
		m_log.SetFileStream(fp);
		if (fp->open(m_slog.c_str()) == false)
		{
			feLogError("Failed creating log file");
			return false;
		}

		// make sure we have a step
		FEAnalysis* step = GetCurrentStep();
		if (step == 0)
		{
			feLogError("No step defined.");
			return false;
		}

		// print welcome message to file
		Logfile::MODE m = m_log.SetMode(Logfile::LOG_FILE);
		febio::Hello(*fp);
		m_log.SetMode(m);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEBioModel::Reset()
{
	// Reset model data
	FEMechModel::Reset();

	// re-initialize the log file
	if (m_logLevel != 0)
	{
		if (InitLogFile() == false) return false;
	}

	// open plot database file
	FEAnalysis* step =  GetCurrentStep();
	if (step->GetPlotLevel() != FE_PLOT_NEVER)
	{
		int hint = step->GetPlotHint();
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(this);
			hint = 0;
		}

		if (hint != FE_PLOT_APPEND)
		{
			if (m_plot->Open(m_splot.c_str()) == false)
			{
				feLogError("Failed creating PLOT database.");
				return false;
			}
		}
	}

	m_stats.ntimeSteps = 0;
	m_stats.ntotalIters = 0;
	m_stats.ntotalRHS = 0;
	m_stats.ntotalReforms = 0;

	// do the callback
	DoCallback(CB_INIT);

	// All data is reset successfully
	return true;
}

//=============================================================================
//                               S O L V E
//=============================================================================

//-----------------------------------------------------------------------------
void FEBioModel::on_cb_solved()
{
	FEAnalysis* step = GetCurrentStep();
	if (step == nullptr) return;

	// for multistep analysis we'll print a grand total
	if (Steps() > 1)
	{
		feLog("\n\n N O N L I N E A R   I T E R A T I O N   S U M M A R Y\n\n");
		feLog("\tNumber of time steps completed .................... : %d\n\n", m_stats.ntimeSteps);
		feLog("\tTotal number of equilibrium iterations ............ : %d\n\n", m_stats.ntotalIters);
		feLog("\tTotal number of right hand evaluations ............ : %d\n\n", m_stats.ntotalRHS);
		feLog("\tTotal number of stiffness reformations ............ : %d\n\n", m_stats.ntotalReforms);
	}

	// get and print elapsed time
	char sztime[64];
	Timer* solveTimer = GetTimer(TimerID::Timer_LinSolve);
	solveTimer->time_str(sztime);
	feLog("\tTime in linear solver: %s\n\n", sztime);

	// always flush the log
	m_log.flush();

	// get peak memory usage
#ifdef WIN32
	size_t memsize = GetPeakMemory();
	if (memsize != 0)
	{
		double mb = (double)memsize / 1048576.0;
		feLog(" Peak memory  : %.1lf MB\n", mb);
	}
#endif

	// print the elapsed time
	GetSolveTimer().time_str(sztime);
	feLog("\n Elapsed time : %s\n\n", sztime);

	// print additional stats to the log file only
	if (m_log.GetMode() & Logfile::LOG_FILE)
	{
		// print more detailed timing info to the log file
		Logfile::MODE old_mode = m_log.SetMode(Logfile::LOG_FILE);

		// sum up all the times spend in the linear solvers
		double total_time   = 0.0;
		double input_time   = m_InputTime.GetTime(); total_time += input_time;
		double init_time    = m_InitTime .GetTime(); total_time += init_time;
		double solve_time   = GetTimer(TimerID::Timer_ModelSolve)->GetTime(); total_time += solve_time;
		double io_time      = m_IOTimer.GetTime();
		double total_linsol = GetTimer(TimerID::Timer_LinSolve )->GetTime();
		double total_reform = GetTimer(TimerID::Timer_Reform   )->GetTime();
		double total_stiff  = GetTimer(TimerID::Timer_Stiffness)->GetTime();
		double total_rhs    = GetTimer(TimerID::Timer_Residual )->GetTime();
		double total_update = GetTimer(TimerID::Timer_Update   )->GetTime();
		double total_qn     = GetTimer(TimerID::Timer_QNUpdate )->GetTime();

		feLog(" T I M I N G   I N F O R M A T I O N\n\n");
		Timer::time_str(input_time  , sztime); feLog("\tInput time ...................... : %s (%lg sec)\n\n", sztime, input_time);
		Timer::time_str(init_time   , sztime); feLog("\tInitialization time ............. : %s (%lg sec)\n\n", sztime, init_time);
		Timer::time_str(solve_time  , sztime); feLog("\tSolve time ...................... : %s (%lg sec)\n\n", sztime, solve_time);
		Timer::time_str(io_time     , sztime); feLog("\t   IO-time (plot, dmp, data) .... : %s (%lg sec)\n\n", sztime, io_time);
		Timer::time_str(total_reform, sztime); feLog("\t   reforming stiffness .......... : %s (%lg sec)\n\n", sztime, total_reform);
		Timer::time_str(total_stiff , sztime); feLog("\t   evaluating stiffness ......... : %s (%lg sec)\n\n", sztime, total_stiff);
		Timer::time_str(total_rhs   , sztime); feLog("\t   evaluating residual .......... : %s (%lg sec)\n\n", sztime, total_rhs);
		Timer::time_str(total_update, sztime); feLog("\t   model update ................. : %s (%lg sec)\n\n", sztime, total_update);
		Timer::time_str(total_qn    , sztime); feLog("\t   QN updates ................... : %s (%lg sec)\n\n", sztime, total_qn);
		Timer::time_str(total_linsol, sztime); feLog("\t   time in linear solver ........ : %s (%lg sec)\n\n", sztime, total_linsol);
		Timer::time_str(total_time  , sztime); feLog("\tTotal elapsed time .............. : %s (%lg sec)\n\n", sztime, total_time);

		m_log.SetMode(old_mode);

		bool bconv = IsSolved();
		if (bconv)
		{
			feLog("\n N O R M A L   T E R M I N A T I O N\n\n");
		}
		else
		{
			feLog("\n E R R O R   T E R M I N A T I O N\n\n");
		}

		// flush the log file
		m_log.flush();
	}

	// close the plot file
	int hint = GetStep(Steps() - 1)->GetPlotHint();
	if (hint != FE_PLOT_APPEND)
		if (m_plot) m_plot->Close();
}

//-----------------------------------------------------------------------------
void FEBioModel::on_cb_stepSolved()
{
	FEAnalysis* step = GetCurrentStep();
	if (step == nullptr) return;

	// output report
	feLog("\n\n N O N L I N E A R   I T E R A T I O N   I N F O R M A T I O N\n\n");
	feLog("\tNumber of time steps completed .................... : %d\n\n", step->m_ntimesteps);
	feLog("\tTotal number of equilibrium iterations ............ : %d\n\n", step->m_ntotiter);
	feLog("\tAverage number of equilibrium iterations .......... : %lg\n\n", (step->m_ntimesteps != 0 ? (double)step->m_ntotiter / (double)step->m_ntimesteps : 0));
	feLog("\tTotal number of right hand evaluations ............ : %d\n\n", step->m_ntotrhs);
	feLog("\tTotal number of stiffness reformations ............ : %d\n\n", step->m_ntotref);

	// print linear solver stats
	FESolver* ps = step->GetFESolver();
	if (ps)
	{
		LinearSolver* ls = step->GetFESolver()->GetLinearSolver();
		if (ls)
		{
			LinearSolverStats stats = ls->GetStats();
			int nsolves = stats.backsolves;
			int niters = stats.iterations;
			double avgiters = (nsolves != 0 ? (double)niters / (double)nsolves : (double)niters);
			feLog("\n L I N E A R   S O L V E R   S T A T S\n\n");
			feLog("\tTotal calls to linear solver ........ : %d\n\n", nsolves);
			feLog("\tAvg iterations per solve ............ : %lg\n\n", avgiters);
		}
	}

	// add to stats
	m_stats.ntimeSteps    += step->m_ntimesteps;
	m_stats.ntotalIters   += step->m_ntotiter;
	m_stats.ntotalRHS     += step->m_ntotrhs;
	m_stats.ntotalReforms += step->m_ntotref;
}
