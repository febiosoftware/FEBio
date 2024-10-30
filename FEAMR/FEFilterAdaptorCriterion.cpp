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
#include "FEFilterAdaptorCriterion.h"

BEGIN_FECORE_CLASS(FEMinMaxFilterAdaptorCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_min, "min");
	ADD_PARAMETER(m_max, "max");
	ADD_PARAMETER(m_clamp, "clamp");

	ADD_PROPERTY(m_data, "data");
END_FECORE_CLASS();

FEMinMaxFilterAdaptorCriterion::FEMinMaxFilterAdaptorCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_min = -1.0e37;
	m_max =  1.0e37;
	m_clamp = true;
	m_data = nullptr;
}

bool FEMinMaxFilterAdaptorCriterion::GetElementValue(FEElement& el, double& value)
{
	if (m_data == nullptr) return false;
	bool b = m_data->GetElementValue(el, value);
	if (b)
	{
		if (m_clamp)
		{
			if (value < m_min) value = m_min;
			if (value > m_max) value = m_max;
		}
		else if ((value < m_min) || (value > m_max))
		{
			b = false;
		}
	}
	return b;
}

bool FEMinMaxFilterAdaptorCriterion::GetMaterialPointValue(FEMaterialPoint& mp, double& value)
{
	if (m_data == nullptr) return false;
	bool b = m_data->GetMaterialPointValue(mp, value);
	if (b)
	{
		if (m_clamp)
		{
			if (value < m_min) value = m_min;
			if (value > m_max) value = m_max;
		}
		else if ((value < m_min) || (value > m_max))
		{
			b = false;
		}
	}
	return b;
}
