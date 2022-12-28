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
#pragma once
#include <FECore/FECallBack.h>
#include "febiorve_api.h"

//-----------------------------------------------------------------------------
class FEBioPlotFile;
class FEMaterialPoint;

//-----------------------------------------------------------------------------
// Base class for RVE probes
class FERVEProbe : public FECallBack
{
	FECORE_BASE_CLASS(FERVEProbe);

public:
	// The first FEModel (fem) is the macro-problem, i.e. the model that will generate the callbacks
	// The second FEModel (rve) is the micro-problem that needs to be tracked.
	FERVEProbe(FEModel* fem);

	bool Init() override;

	bool Execute(FEModel& fem, int nwhen) override;

	void Save();

	void SetDebugFlag(bool b) { m_bdebug = b; }
	bool GetDebugFlag() const { return m_bdebug; }

	void SetRVEModel(FEModel* rve);

protected:
	std::string	m_file;			//!< file name
	bool		m_bdebug;		//!< debug flag

private:
	FEBioPlotFile* m_xplt;		//!< the actual plot file
	FEModel* m_rve;		//!< the RVE model
};

//-----------------------------------------------------------------------------
// Probe used by FEMicroMaterial
class FEMicroProbe : public FERVEProbe
{
public:
	FEMicroProbe(FEModel* fem);

	bool Init() override;

private:
	int			m_neid;			//!< element Id
	int			m_ngp;			//!< Gauss-point (one-based!)

	DECLARE_FECORE_CLASS();
};
