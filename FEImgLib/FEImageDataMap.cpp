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
#include "FEImageDataMap.h"
#include <FECore/FEDomainMap.h>
#include <FECore/FEMesh.h>
#include "image_tools.h"

BEGIN_FECORE_CLASS(FEImageDataMap, FEElemDataGenerator)
	ADD_PARAMETER(m_r0, "range_min");
	ADD_PARAMETER(m_r1, "range_max");
	ADD_PARAMETER(m_blur, "blur");
	ADD_PROPERTY(m_imgSrc, "image");
END_FECORE_CLASS();

FEImageDataMap::FEImageDataMap(FEModel* fem) : FEElemDataGenerator(fem), m_map(m_im)
{
	m_imgSrc = nullptr;
	m_blur = 0.0;
	m_data = nullptr;
}

bool FEImageDataMap::Init()
{
	if (m_imgSrc == nullptr) return false;

	if (m_imgSrc->GetImage3D(m_im0) == false)
	{
		return false;
	}
	m_im = m_im0;
	m_map.SetRange(m_r0, m_r1);

	return FEElemDataGenerator::Init();
}

FEDataMap* FEImageDataMap::Generate()
{
	FEElementSet* elset = GetElementSet();
	if (elset == nullptr) return nullptr;

	// TODO: Can I use FMT_NODE?
	if (m_data == nullptr) m_data = new FEDomainMap(FEDataType::FE_DOUBLE, Storage_Fmt::FMT_MULT);
	m_data->Create(elset);

	GenerateData();
	return m_data;
}

void FEImageDataMap::GenerateData()
{
	assert(m_data);
	FEElementSet& set = *GetElementSet();
	FEMesh& mesh = *set.GetMesh();
	int N = set.Elements();
	for (int i = 0; i < N; ++i)
	{
		FEElement& el = set.Element(i);
		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j)
		{
			vec3d ri = mesh.Node(el.m_node[j]).m_r0;
			double d = m_map.value(m_map.map(ri));
			m_data->setValue(i, j, d);
		}
	}
}

void FEImageDataMap::Evaluate(double time)
{
	if (m_blur > 0)
	{
		if (m_im0.depth() == 1) blur_image_2d(m_im, m_im0, (float)m_blur);
		else blur_image(m_im, m_im0, (float)m_blur);
	}
	else m_im = m_im0;

	GenerateData();
}
