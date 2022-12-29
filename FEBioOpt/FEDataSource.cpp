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
#include <cwctype>
#include "FEDataSource.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/NodeDataRecord.h>
#include <FECore/ElementDataRecord.h>
#include <FECore/SurfaceDataRecord.h>
#include <FECore/DomainDataRecord.h>

//=================================================================================================
FEDataSource::FEDataSource(FEModel* fem) : m_fem(*fem)
{
	
}

FEDataSource::~FEDataSource()
{
	
}

bool FEDataSource::Init()
{
	return true;
}

void FEDataSource::Reset()
{
	
}

//=================================================================================================
bool FEDataParameter::update(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEDataParameter& src = *((FEDataParameter*)pd);
	src.update();

	return true;
}

void FEDataParameter::update()
{
	// get the current time value
	double time = m_fem.GetTime().currentTime;

	// evaluate the current reaction force value
	double x = m_fx();
	double y = m_fy();

	// add the data pair to the loadcurve
	m_rf.Add(x, y);
}

FEDataParameter::FEDataParameter(FEModel* fem) : FEDataSource(fem)
{
	m_ord = "fem.time";
}

void FEDataParameter::SetParameterName(const std::string& name)
{
	m_param = name;
}

// set the ordinate name
void FEDataParameter::SetOrdinateName(const std::string& name)
{
	m_ord = name;
}

bool FEDataParameter::Init()
{
	FEModel* fem = &m_fem;

	// find all the parameters
	FEParamValue val = m_fem.GetParameterValue(ParamString(m_param.c_str()));
	if (val.isValid() == false) {

		// see if it's a data parameter
		if (strstr(m_param.c_str(), "fem.node_data"))
		{
			char buf[256] = { 0 };
			strcpy(buf, m_param.c_str());
			char* sz = buf + 14;
			char* c1 = strchr(sz, ',');
			*c1++ = 0;

			int nid = atoi(c1);

			c1 = strrchr(sz, '\'');
			if (sz[0] == '\'') sz++;
			*c1 = 0;

			FELogNodeData* pd = fecore_new<FELogNodeData>(sz, fem);
			if (pd == nullptr) { feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str()); return false; }

			FEMesh& mesh = fem->GetMesh();
			FENode* pn = mesh.FindNodeFromID(nid);
			if (pn == nullptr) { feLogErrorEx(fem, "Invalid node id"); return false; }

			m_fy = [=]() { return pd->value(*pn); };
		}
		else if (strstr(m_param.c_str(), "fem.element_data"))
		{
			char buf[256] = { 0 };
			strcpy(buf, m_param.c_str());
			char* sz = buf + 17;
			char* c1 = strchr(sz, ',');
			*c1++ = 0;

			int eid = atoi(c1);

			c1 = strrchr(sz, '\'');
			if (sz[0] == '\'') sz++;
			*c1 = 0;

			FELogElemData* pd = fecore_new<FELogElemData>(sz, fem);
			if (pd == nullptr) { feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str()); return false; }

			FEMesh& mesh = fem->GetMesh();
			FEElement* pe = mesh.FindElementFromID(eid);
			if (pe == nullptr) { feLogErrorEx(fem, "Invalid element id"); return false; }

			m_fy = [=]() { return pd->value(*pe); };
		}
		else if (strstr(m_param.c_str(), "fem.surface_data"))
		{
			char buf[256] = { 0 };
			strcpy(buf, m_param.c_str());
			char* sz = buf + 17;
			char* c1 = strchr(sz, ',');
			*c1++ = 0;

			char* szsurf = strchr(c1, '\'');
			if (szsurf == 0) { feLogErrorEx(fem, "Syntax error in parameter name %s", m_param.c_str()); return false; }
			szsurf++;

			c1 = strrchr(sz, '\'');
			if (sz[0] == '\'') sz++;
			*c1 = 0;

			c1 = strrchr(szsurf, '\'');
			if (c1) *c1 = 0;

			FELogSurfaceData* pd = fecore_new<FELogSurfaceData>(sz, fem);
			if (pd == nullptr) { feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str()); return false; }

			FEMesh& mesh = fem->GetMesh();
			FESurface* surf = mesh.FindSurface(szsurf);
			if (surf == nullptr) { feLogErrorEx(fem, "Invalid surface name %s", szsurf); return false; }

			m_fy = [=]() { return pd->value(*surf); };
		}
		else if (strstr(m_param.c_str(), "fem.domain_data"))
		{
			const char* c = m_param.c_str() + 15;
			char buf[256] = { 0 };
			int n = 0;
			bool skipws = true;
			while (c && *c)
			{
				if (*c == '\'')
				{
					if (skipws) skipws = false; else skipws = true;
				}

				if ((skipws == false) || (iswspace(*c) == 0)) buf[n++] = *c;
				c++;
			}
			buf[n] = 0;

			char* sz = buf;
			int l = strlen(sz);
			if (sz[0  ] != '(') { feLogErrorEx(fem, "Syntax error in parameter name %s", m_param.c_str()); return false; }
			if (sz[l-1] != ')') { feLogErrorEx(fem, "Syntax error in parameter name %s", m_param.c_str()); return false; }
			sz[l - 1] = 0;
			*sz++ = 0;

			// separate string in data variable and domain name
			char* c1 = strrchr(sz, ',');
			*c1++ = 0;

			// process part name
			char* szdom = strchr(c1, '\'');
			if (szdom == 0) { feLogErrorEx(fem, "Syntax error in parameter name %s", m_param.c_str()); return false; }
			szdom++;

			c1 = strrchr(szdom, '\'');
			if (c1) *c1 = 0; else { feLogErrorEx(fem, "Syntax error in parameter name %s", m_param.c_str()); return false; }

			FEMesh& mesh = fem->GetMesh();
			FEDomain* dom = mesh.FindDomain(szdom);
			if (dom == nullptr) { feLogErrorEx(fem, "Invalid domain name %s", szdom); return false; }

			// process data name
			FELogDomainData* pd = nullptr;
			if (sz[0] == '\'')
			{
				sz++;
				c1 = strrchr(sz, '\'');
				if (c1 == nullptr) { feLogErrorEx(fem, "Invalid domain name %s", szdom); return false; }
				*c1 = 0;
				pd = fecore_new<FELogDomainData>(sz, fem);
			}
			else
			{
				c1 = strchr(sz, '(');
				if (c1 == nullptr) { feLogErrorEx(fem, "Syntax error"); return false; }
				char* cvar = c1 + 1;
				*c1 = 0;
				if ((c1 = strrchr(cvar, ')')) == nullptr) { feLogErrorEx(fem, "Syntax error"); return false; }
				*c1 = 0;

				std::vector<std::string> varList;
				do
				{
					if (cvar[0] == '\'')
					{
						*cvar++ = 0;
						if ((c1 = strchr(cvar, '\'')) == nullptr) { feLogErrorEx(fem, "Syntax error"); return false; }
						*c1 = 0;
						varList.push_back(cvar);
						cvar = c1 + 1;
					}
					else varList.push_back(cvar);
					cvar = strchr(cvar, ',');
					if (cvar) cvar++;
				} while (cvar && *cvar);

				pd = fecore_new<FELogDomainData>(sz, fem);
				if (pd == nullptr) { feLogErrorEx(fem, "Syntax error"); return false; }
				if (pd->SetParameters(varList) == false) { feLogErrorEx(fem, "Syntax error"); return false; }
			}
			if (pd == nullptr) { feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str()); return false; }

			m_fy = [=]() { return pd->value(*dom); };
		}
		else
		{
			feLogErrorEx(fem, "Invalid parameter name %s", m_param.c_str());
			return false;
		}
	}
	else
	{
		if (val.type() != FE_PARAM_DOUBLE) {
			feLogErrorEx(fem, "Invalid type for parameter %s", m_param.c_str());
			return false;
		}
		m_fy = [=]() { return val.value<double>(); };
	}

	// find the ordinate
	val = m_fem.GetParameterValue(ParamString(m_ord.c_str()));
	if (val.isValid() == false) {

		// see if it's a data parameter
        if (strstr(m_ord.c_str(), "fem.node_data"))
        {
            char buf[256] = { 0 };
            strcpy(buf, m_ord.c_str());
            char* sz = buf + 14;
            char* c1 = strchr(sz, ',');
            *c1++ = 0;
            
            int nid = atoi(c1);
            
            c1 = strrchr(sz, '\'');
            if (sz[0] == '\'') sz++;
            *c1 = 0;
            
            FELogNodeData* pd = fecore_new<FELogNodeData>(sz, fem);
            if (pd == nullptr) { feLogErrorEx(fem, "Invalid ordinate name %s", m_ord.c_str()); return false; }
            
            FEMesh& mesh = fem->GetMesh();
            FENode* pn = mesh.FindNodeFromID(nid);
            if (pn == nullptr) { feLogErrorEx(fem, "Invalid node id"); return false; }
            
            m_fx = [=]() { return pd->value(*pn); };
        }
		else if (strstr(m_ord.c_str(), "fem.element_data"))
		{
			char buf[256] = { 0 };
			strcpy(buf, m_ord.c_str());
			char* sz = buf + 17;
			char* c1 = strchr(sz, ',');
			*c1++ = 0;

			int eid = atoi(c1);

			c1 = strrchr(sz, '\'');
			if (sz[0] == '\'') sz++;
			*c1 = 0;

			FELogElemData* pd = fecore_new<FELogElemData>(sz, fem);
			if (pd == nullptr) { feLogErrorEx(fem, "Invalid ordinate name %s", m_ord.c_str()); return false; }

			FEMesh& mesh = fem->GetMesh();
			FEElement* pe = mesh.FindElementFromID(eid);
			if (pe == nullptr) { feLogErrorEx(fem, "Invalid element id"); return false; }

			m_fx = [=]() { return pd->value(*pe); };
		}
		else
		{
			feLogErrorEx(fem, "Invalid ordinate name %s", m_ord.c_str());
			return false;
		}
	}
	else
	{
		if (val.type() != FE_PARAM_DOUBLE) {
			feLogErrorEx(fem, "Invalid type for ordinate %s", m_ord.c_str());
			return false;
		}
		m_fx = [=]() { return val.value<double>(); };
	}

	// register callback
	m_fem.AddCallback(update, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return FEDataSource::Init();
}

void FEDataParameter::Reset()
{
	// reset the reaction force load curve
	m_rf.Clear();
	FEDataSource::Reset();
}

double FEDataParameter::Evaluate(double x)
{
	return m_rf.value(x);
}

//=================================================================================================
FEDataFilterPositive::FEDataFilterPositive(FEModel* fem) : FEDataSource(fem)
{
	m_src = 0;
}

FEDataFilterPositive::~FEDataFilterPositive()
{
	if (m_src) delete m_src;
}

void FEDataFilterPositive::SetDataSource(FEDataSource* src)
{
	if (m_src) delete m_src;
	m_src = src;
}

bool FEDataFilterPositive::Init()
{
	if (m_src == 0) return false;

	return m_src->Init();
}

void FEDataFilterPositive::Reset()
{
	if (m_src) m_src->Reset();
}

double FEDataFilterPositive::Evaluate(double t)
{
	double v = m_src->Evaluate(t);
	return (v >= 0.0 ? v : -v);
}


//=================================================================================================
FEDataFilterSum::FEDataFilterSum(FEModel* fem) : FEDataSource(fem)
{
	m_data = nullptr;
	m_nodeSet = nullptr;
}

FEDataFilterSum::~FEDataFilterSum()
{
	delete m_data;
}

void FEDataFilterSum::SetData(FELogNodeData* data, FENodeSet* nodeSet)
{
	m_data = data;
	m_nodeSet = nodeSet;
}

// Initialize data
bool FEDataFilterSum::Init()
{
	if (m_data == nullptr) return false;
	if (m_nodeSet == nullptr) return false;

	if (m_data->Init() == false) return false;

	// register callback
	m_fem.AddCallback(update, CB_MAJOR_ITERS, (void*)this);

	return FEDataSource::Init();
}

// reset data
void FEDataFilterSum::Reset()
{
	m_rf.Clear();
	m_rf.Add(0, 0);
}

// evaluate data source at x
double FEDataFilterSum::Evaluate(double x)
{
	return m_rf.value(x);
}

bool FEDataFilterSum::update(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEDataFilterSum& src = *((FEDataFilterSum*)pd);
	src.update();

	return true;
}

void FEDataFilterSum::update()
{
	// get the current time value
	double time = m_fem.GetTime().currentTime;

	FEMesh* mesh = m_nodeSet->GetMesh();
	FENodeSet& ns = *m_nodeSet;
	double sum = 0.0;
	for (int i = 0; i < m_nodeSet->Size(); ++i)
	{
		double vi = m_data->value(*ns.Node(i));
		sum += vi;
	}

	// evaluate the current reaction force value
	double x = time;
	double y = sum;

	// add the data pair to the loadcurve
	m_rf.Add(x, y);
}
