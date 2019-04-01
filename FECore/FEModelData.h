/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FECoreBase.h"
#include "ElementDataRecord.h"


//-------------------------------------------------------------------
// Class that can used to query model data
class FECORE_API FEModelData : public FECoreBase
{
	FECORE_SUPER_CLASS

	struct ELEMREF
	{
		int	ndom;
		int	nid;
	};

public:
	// constructor
	FEModelData(FEModel* fem, FELogElemData* eval, vector<int>& items);

	// update model data (i.e. evaluate it based on current state of model)
	void Update();

private:
	FELogElemData*		m_eval;		//!< class that evaluates data
	std::vector<int>	m_item;		//!< item list
	std::vector<double>	m_data;		//!< model data values

protected:
	void BuildELT();
	vector<ELEMREF>	m_ELT;
	int	m_offset;

	DECLARE_FECORE_CLASS();
};
