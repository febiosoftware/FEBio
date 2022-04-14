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
#include "FEBioControlSection3.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FENewtonSolver.h"

class FEObsoleteStepParamHandler : public FEObsoleteParamHandler
{
public:
	FEObsoleteStepParamHandler(XMLTag& tag, FEAnalysis* step) : m_step(step), FEObsoleteParamHandler(tag, step)
	{
		AddParam("solver.max_ups"           , "solver.qn_method.max_ups"        , FE_PARAM_INT);
		AddParam("solver.qn_max_buffer_size", "solver.qn_method.max_buffer_size", FE_PARAM_INT);
		AddParam("solver.qn_cycle_buffer"   , "solver.qn_method.cycle_buffer"   , FE_PARAM_BOOL);
		AddParam("solver.cmax"              , "solver.qn_method.cmax"           , FE_PARAM_DOUBLE);
	}

	bool ProcessTag(XMLTag& tag) override
	{
		if (tag == "qnmethod")
		{
			const char* szv = tag.szvalue();
			int l = strlen(szv);
			if      ((strnicmp(szv, "BFGS"   , l) == 0) || (strnicmp(szv, "0", l) == 0)) m_qnmethod = QN_BFGS;
			else if ((strnicmp(szv, "BROYDEN", l) == 0) || (strnicmp(szv, "1", l) == 0)) m_qnmethod = QN_BROYDEN;
			else if ((strnicmp(szv, "JFNK"   , l) == 0) || (strnicmp(szv, "2", l) == 0)) m_qnmethod = QN_JFNK;
			else return false;

			return true;
		}
		else return FEObsoleteParamHandler::ProcessTag(tag);
	}

	void MapParameters() override
	{
		FEModel* fem = m_step->GetFEModel();

		// first, make sure that the QN method is allocated
		if (m_qnmethod != -1)
		{
			FENewtonSolver& solver = dynamic_cast<FENewtonSolver&>(*m_step->GetFESolver());
			FEProperty& qn = *solver.FindProperty("qn_method");
			switch (m_qnmethod)
			{
			case QN_BFGS   : solver.SetSolutionStrategy(fecore_new<FENewtonStrategy>("BFGS"   , fem)); break;
			case QN_BROYDEN: solver.SetSolutionStrategy(fecore_new<FENewtonStrategy>("Broyden", fem)); break;
			case QN_JFNK   : solver.SetSolutionStrategy(fecore_new<FENewtonStrategy>("JFNK"   , fem)); break;
			default:
				assert(false);
			}
		}

		// now, process the rest
		FEObsoleteParamHandler::MapParameters();
	}

private:
	int	m_qnmethod = -1;
	FEAnalysis* m_step;
};

//-----------------------------------------------------------------------------
FEBioControlSection3::FEBioControlSection3(FEFileImport* pim) : FEFileSection(pim)
{
}

//-----------------------------------------------------------------------------
void FEBioControlSection3::Parse(XMLTag& tag)
{
	// get the step
	FEAnalysis* pstep = GetBuilder()->GetStep();
	if (pstep == 0)
	{
		throw XMLReader::InvalidTag(tag);
	}

	// Get the solver
	FESolver* psolver = pstep->GetFESolver();
	if (psolver == 0) 
	{
		string m = GetBuilder()->GetModuleName();
		throw FEBioImport::FailedAllocatingSolver(m.c_str());
	}

	// prepare obsolete parameter mapping.
	FEObsoleteStepParamHandler stepParamHandler(tag, pstep);

	// read the step parameters
	SetInvalidTagHandler(&stepParamHandler);
	ReadParameterList(tag, pstep);

	stepParamHandler.MapParameters();
}
