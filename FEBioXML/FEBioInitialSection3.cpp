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

#include "stdafx.h"
#include "FEBioInitialSection3.h"
#include <FECore/FEInitialCondition.h>

FEBioInitialSection3::FEBioInitialSection3(FEFileImport* pim) : FEFileSection(pim) 
{
}

void FEBioInitialSection3::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();

	++tag;
	do
	{
		if (tag == "ic")
		{
			// read the type attribute
			const char* sztype = tag.AttributeValue("type");

			// try to allocate the initial condition
			FEInitialCondition* pic = fecore_new<FEInitialCondition>(sztype, fem);
			if (pic == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// get the node set
			const char* szset = tag.AttributeValue("node_set", true);
			if (szset)
			{
				FENodeSet* nodeSet = mesh.FindNodeSet(szset);
				if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

				// add the nodes to the IC
				pic->AddNodes(*nodeSet);
			}

			// add it to the model
			GetBuilder()->AddInitialCondition(pic);

			// Read the parameter list
			ReadParameterList(tag, pic);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
