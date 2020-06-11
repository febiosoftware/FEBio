/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEOptimizeData.h"
#include "FELMOptimizeMethod.h"
#include "FEPowellOptimizeMethod.h"
#include "FEScanOptimizeMethod.h"
#include "FEConstrainedLMOptimizeMethod.h"
#include "FEOptimizeInput.h"
#include <FECore/log.h>
#include <FEBioXML/XMLReader.h>

//=============================================================================
// FEOptimizeInput
//=============================================================================

//-----------------------------------------------------------------------------
//! Read the data from the xml input file
//!
bool FEOptimizeInput::Input(const char* szfile, FEOptimizeData* pOpt)
{
	// try to open the file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "\nFATAL ERROR: Failed to load file %s\n", szfile);
		return false;
	}

	// find the root element
	XMLTag tag;
	if (xml.FindTag("febio_optimize", tag) == false) return false;

	// check for the version attribute
	const char* szversion = tag.AttributeValue("version", true);
	if ((szversion == 0) || (strcmp(szversion, "2.0") != 0))
	{
		fprintf(stderr, "\nFATAL ERROR: Invalid version number for febio_optimize!\n\n");
		return false;
	}

	m_opt = pOpt;

	// build the file section map
	m_map["Options"    ] = new FEOptionsSection(pOpt, this);
	m_map["Parameters" ] = new FEParametersSection(pOpt, this);
	m_map["Constraints"] = new FEConstraintsSection(pOpt, this);
	m_map["Objective"  ] = new FEObjectiveSection(pOpt, this);
	m_map["Task"       ] = new FETaskSection(pOpt, this);

	// parse the file
	bool ret = ParseFile(tag);

	// all done
	xml.Close();

	return ret;
}

//=================================================================================================
void FEOptionsSection::Parse(XMLTag& tag)
{
	FEOptimizeMethod* popt = 0;
	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) popt = new FELMOptimizeMethod;
	else
	{
		if (strcmp(szt, "levmar") == 0) popt = new FELMOptimizeMethod;
		else if (strcmp(szt, "powell") == 0) popt = new FEPowellOptimizeMethod;
		else if (strcmp(szt, "scan") == 0) popt = new FEScanOptimizeMethod;
#ifdef HAVE_LEVMAR
		else if (strcmp(szt, "constrained levmar") == 0) popt = new FEConstrainedLMOptimizeMethod;
#endif
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}

	// get the parameter list
	FEParameterList& pl = popt->GetParameterList();

	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (ReadParameter(tag, pl) == false)
			{
				if (tag == "log_level")
				{
					char szval[256];
					tag.value(szval);
					if (strcmp(szval, "LOG_DEFAULT") == 0) {} // don't change the plot level
					else if (strcmp(szval, "LOG_NEVER") == 0) popt->m_loglevel = LogLevel::LOG_NEVER;
					else if (strcmp(szval, "LOG_FILE_ONLY") == 0) popt->m_loglevel = LogLevel::LOG_FILE;
					else if (strcmp(szval, "LOG_SCREEN_ONLY") == 0) popt->m_loglevel = LogLevel::LOG_SCREEN;
					else if (strcmp(szval, "LOG_FILE_AND_SCREEN") == 0) popt->m_loglevel = LogLevel::LOG_FILE_AND_SCREEN;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "print_level")
				{
					char szval[256];
					tag.value(szval);
					if      (strcmp(szval, "PRINT_ITERATIONS") == 0) popt->m_print_level = PRINT_ITERATIONS;
					else if (strcmp(szval, "PRINT_VERBOSE"   ) == 0) popt->m_print_level = PRINT_VERBOSE;
					else
					{
						int print_level = atoi(szval);
						if ((print_level == PRINT_ITERATIONS) || (print_level == PRINT_VERBOSE))
						{
							popt->m_print_level = print_level;
						}
						else throw XMLReader::InvalidValue(tag);
					}
				}
				else throw XMLReader::InvalidTag(tag);
			}
			++tag;
		} while (!tag.isend());
	}

	m_opt->SetSolver(popt);
}

//=================================================================================================
void FETaskSection::Parse(XMLTag& tag)
{
	m_opt->m_pTask = fecore_new<FECoreTask>(tag.szvalue(), m_opt->GetFEModel());
	if (m_opt->m_pTask == nullptr) throw XMLReader::InvalidValue(tag);
}

//=================================================================================================
void FEObjectiveSection::Parse(XMLTag& tag)
{
	FEModel& fem = *m_opt->GetFEModel();

	// get the type attribute
	const char* sztype = tag.AttributeValue("type");

	if (strcmp(sztype, "data-fit") == 0)
	{
		FEDataFitObjective* obj = new FEDataFitObjective(&fem);
		m_opt->SetObjective(obj);

		++tag;
		do
		{
			if (tag == "fnc")
			{
				FEDataSource* src = ParseDataSource(tag, *m_opt);
				obj->SetDataSource(src);
			}
			else if (tag == "data")
			{
				vector<pair<double, double> > data;

				// see if the user wants to read the data from a text file
				const char* szf = tag.AttributeValue("import", true);
				if (szf)
				{
					// make sure this tag is a leaf
					if ((tag.isempty() == false) || (tag.isleaf() == false)) throw XMLReader::InvalidValue(tag);

					// read the data form a text file
					FILE* fp = fopen(szf, "rt");
					if (fp == 0) throw XMLReader::InvalidAttributeValue(tag, "import", szf);
					char szline[256] = { 0 };
					do
					{
						fgets(szline, 255, fp);
						double t, v;
						int n = sscanf(szline, "%lg%lg", &t, &v);
						if (n == 2)
						{
							pair<double, double> pt(t, v);
							data.push_back(pt);
						}
						else break;
					} while ((feof(fp) == 0) && (ferror(fp) == 0));

					fclose(fp);
				}
				else
				{
					double v[2] = { 0 };
					++tag;
					do
					{
						tag.value(v, 2);
						data.push_back(pair<double, double>(v[0], v[1]));
						++tag;
					} while (!tag.isend());
				}

				obj->SetMeasurements(data);
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		} while (!tag.isend());
	}
	else if (strcmp(sztype, "target") == 0)
	{
		FEMinimizeObjective* obj = new FEMinimizeObjective(&fem);
		m_opt->SetObjective(obj);

		++tag;
		do
		{
			if (tag == "var")
			{
				const char* szname = tag.AttributeValue("name");

				double d[2] = { 0 };
				tag.value(d, 2);

				if (obj->AddFunction(szname, d[0]) == false) throw XMLReader::InvalidAttributeValue(tag, "name", szname);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}
	else if (strcmp(sztype, "element-data") == 0)
	{
		FEElementDataTable* obj = new FEElementDataTable(&fem);
		m_opt->SetObjective(obj);

		++tag;
		do
		{
			if (tag == "var")
			{
				const char* sztype = tag.AttributeValue("type");

				// try to allocate the element data record
				FELogElemData* var = fecore_new<FELogElemData>(sztype, &fem);
				if (var == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

				obj->SetVariable(var);
			}
			else if (tag == "data")
			{
				++tag;
				do
				{
					if (tag == "elem")
					{
						const char* szid = tag.AttributeValue("id");
						int nid = atoi(szid);
						double v = 0.0;
						tag.value(v);

						obj->AddValue(nid, v);
					}
					else throw XMLReader::InvalidTag(tag);

					++tag;
				}
				while (!tag.isend());
			}
			++tag;
		}
		while (!tag.isend());
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	FEOptimizeMethod* solver = m_opt->GetSolver();
	if (solver)
	{
		FEObjectiveFunction& obj = m_opt->GetObjective();
		if (solver->m_print_level == PRINT_ITERATIONS) obj.SetVerbose(false);
		else obj.SetVerbose(true);
	}
}

//-----------------------------------------------------------------------------
FEDataSource* FEObjectiveSection::ParseDataSource(XMLTag& tag, FEOptimizeData& opt)
{
	FEModel& fem = *opt.GetFEModel();

	const char* sztype = tag.AttributeValue("type");
	if (strcmp(sztype, "parameter") == 0)
	{
		FEDataParameter* src = new FEDataParameter(&fem);
		++tag;
		do
		{
			if (tag == "param")
			{
				const char* szname = tag.AttributeValue("name");
				src->SetParameterName(szname);
			}
			else if (tag == "ordinate")
			{
				const char* szname = tag.AttributeValue("name");
				src->SetOrdinateName(szname);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());

		return src;
	}
	else if (strcmp(sztype, "filter_positive_only") == 0)
	{
		FEDataFilterPositive* src = new FEDataFilterPositive(&fem);

		++tag;
		do
		{
			if (tag == "source")
			{
				FEDataSource* s = ParseDataSource(tag, opt);
				src->SetDataSource(s);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());

		return src;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// we shouldn't be here
	assert(false);
	return 0;
}

//=============================================================================
void FEParametersSection::Parse(XMLTag& tag)
{
	FEModel& fem = *m_opt->GetFEModel();

	// read the parameters
	++tag;
	do
	{
		if (tag == "param")
		{
			FEModelParameter* var = new FEModelParameter(&fem);

			// get the variable name
			const char* sz = tag.AttributeValue("name");
			var->SetName(sz);

			// set initial values and bounds
			double d[4] = { 0, 0, 0, 1 };
			tag.value(d, 4);
			var->InitValue() = d[0];
			var->MinValue() = d[1];
			var->MaxValue() = d[2];
			var->ScaleFactor() = d[3];

			// add the variable
			m_opt->AddInputParameter(var);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

//=============================================================================
void FEConstraintsSection::Parse(XMLTag& tag)
{
	int NP = m_opt->InputParameters();
	if ((NP > OPT_MAX_VAR) || (NP < 2)) throw XMLReader::InvalidTag(tag);

	double v[OPT_MAX_VAR + 1];
	++tag;
	do
	{
		if (tag == "constraint")
		{
			int m = tag.value(v, OPT_MAX_VAR + 1);
			if (m != NP + 1) throw XMLReader::InvalidValue(tag);

			OPT_LIN_CONSTRAINT con;
			for (int i = 0; i<NP; ++i) con.a[i] = v[i];
			con.b = v[NP];

			m_opt->AddLinearConstraint(con);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}
