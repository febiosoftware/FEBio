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
#include "FEBioParamRun.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEShellDomain.h>
#include <FECore/log.h>
#include <XML/XMLReader.h>

//! class constructor
FEBioParamRun::FEBioParamRun(FEModel* pfem) : FECoreTask(pfem)
{
	m_febioOutput = false;
}

//! initialization
bool FEBioParamRun::Init(const char* szfile)
{
	feLog("P A R A M   R U N   M O D U L E\n\n");

	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	// read the input file
	if (Input(szfile) == false) return false;

	// don't plot anything
	if (m_febioOutput == false)
	{
		// NOTE: I need to call GetParameterList to ensure that the parameters are allocated.
		FEParameterList& pl = fem->GetParameterList();
		fem->SetParameter("log_level", 0);
		for (int i = 0; i < fem->Steps(); ++i)
		{
			fem->GetStep(i)->SetPlotLevel(FE_PLOT_NEVER);
			fem->GetStep(i)->SetOutputLevel(FE_OUTPUT_NEVER);
		}
	}

	// do the initialization of the task
	if (m_febioOutput == false) GetFEModel()->BlockLog();
	if (fem->Init() == false) return false;
	if (m_febioOutput == false) GetFEModel()->UnBlockLog();

	// initialize all parameters
	for (FEModelParameter* v : m_inVar)
	{
		if (v->Init() == false) return false;

		// set the initial value
		v->SetValue(v->InitValue());
	}
	for (FEDataParameter* v : m_outVar)
	{
		if (v->Init() == false) return false;
	}

	// since we can change shell thickness now, we need to reinitialize 
	// the shell elements with the new thickness
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEShellDomainNew* shellDomain = dynamic_cast<FEShellDomainNew*>(&mesh.Domain(i));
		if (shellDomain) shellDomain->AssignDefaultShellThickness();
	}
	fem->InitShells();

	return true;
}

//! read control file
bool FEBioParamRun::Input(const char* szfile)
{
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	// find the root tag
	XMLTag tag;
	if (xml.FindTag("febio_run", tag) == false) return false;

	// read the tags
	++tag;
	do
	{
		if (tag == "Parameters")
		{
			++tag;
			do
			{
				if (tag == "param")
				{
					// read parameter
					FEModelParameter* var = new FEModelParameter(fem);

					// get the variable name
					const char* sz = tag.AttributeValue("name");
					var->SetName(sz);

					// set the value
					double val = 0.0;
					tag.value(val);
					var->InitValue() = val;

					m_inVar.push_back(var);
				}
				else return false;
				++tag;
			} while (!tag.isend());
		}
		else if (tag == "Output")
		{
			++tag;
			do
			{
				if (tag == "file")
				{
					tag.value(m_outFile);
				}
				else if (tag == "generate_febio_output")
				{
					tag.value(m_febioOutput);
				}
				else if (tag == "param")
				{
					// read parameter
					FEDataParameter* var = new FEDataParameter(fem);

					// get the variable name
					const char* sz = tag.AttributeValue("name");
					var->SetParameterName(sz);

					m_outVar.push_back(var);
				}
				else return false;
				++tag;
			}
			while (!tag.isend());
		}
		else return false;
		++tag;
	} while (!tag.isend());

	// cleanup
	xml.Close();
	return true;
}

//! Run the task
bool FEBioParamRun::Run()
{
	// get the model
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	// solve the model
	if (m_febioOutput == false) GetFEModel()->BlockLog();
	if (fem->Solve() == false) return false;
	if (m_febioOutput == false) GetFEModel()->UnBlockLog();

	// output the values
	FILE* fp = fopen(m_outFile.c_str(), "wt");
	for (FEDataParameter* p : m_outVar)
	{
		double v = p->value();
		fprintf(fp, "%lg\n", v);
	}
	fclose(fp);

	return true;
}
