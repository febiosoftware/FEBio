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
#include "FEBioMaterialSection.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEMaterial.h"
#include <FECore/log.h>
#include <sstream>

//-----------------------------------------------------------------------------
//! This function creates a material by checking the type attribute against
//! registered materials. Also, if the tag defines attributes (other than
//! type and name), the material is offered a chance to process the attributes.
FEMaterial* FEBioMaterialSection::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);
	
	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// create a new material of this type
	FEMaterial* pmat = GetBuilder()->CreateMaterial(sztype);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	return pmat;
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		if (tag == "material")
		{
			// check that the ID attribute is defined and that it 
			// equals the number of materials + 1.
			int nid = -1;
			tag.AttributeValue("id", nid);
			int nmat = fem.Materials();
			if (nid != nmat+1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// make sure that the name is unique
			std::string name;
			const char* szname = tag.AttributeValue("name", true);
			if (szname == nullptr)
			{
				stringstream ss;
				ss << "Material" << nid;
				name = ss.str();

				feLogWarningEx((&fem), "Material %d has no name.\nIt was given the name %s.", nid, name.c_str());
			}
			else name = szname;

			FEMaterial* mat = fem.FindMaterial(name);
			if (mat)
			{
				throw XMLReader::InvalidAttributeValue(tag, "name");
			}

			// create a material from this tag
			FEMaterial* pmat = CreateMaterial(tag); assert(pmat);
			pmat->SetName(name);

			// set the material's ID
			++m_nmat;
			pmat->SetID(m_nmat);

			// parse the material parameters
			ReadParameterList(tag, pmat);

			// add the material
			GetBuilder()->AddMaterial(pmat);
		}
		else throw XMLReader::InvalidTag(tag);

		// read next tag
		++tag;
	}
	while (!tag.isend());
}

//===============================================================================

//-----------------------------------------------------------------------------
//! This function creates a material by checking the type attribute against
//! registered materials. Also, if the tag defines attributes (other than
//! type and name), the material is offered a chance to process the attributes.
FEMaterial* FEBioMaterialSection3::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);

	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// create a new material of this type
	FEMaterial* pmat = GetBuilder()->CreateMaterial(sztype);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	return pmat;
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
void FEBioMaterialSection3::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		if (tag == "material")
		{
			// check that the ID attribute is defined and that it 
			// equals the number of materials + 1.
			int nid = -1;
			tag.AttributeValue("id", nid);
			int nmat = fem.Materials();
			if (nid != nmat + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// make sure that the name is unique
			const char* szname = tag.AttributeValue("name");
			FEMaterial* mat = fem.FindMaterial(szname);
			if (mat)
			{
				throw XMLReader::InvalidAttributeValue(tag, "name");
			}

			// create a material from this tag
			FEMaterial* pmat = CreateMaterial(tag); assert(pmat);
			pmat->SetName(szname);

			// add the material
			fem.AddMaterial(pmat);
			++m_nmat;

			// set the material's ID
			pmat->SetID(m_nmat);

			// parse the material parameters
			ReadParameterList(tag, pmat);
		}
		else throw XMLReader::InvalidTag(tag);

		// read next tag
		++tag;
	} while (!tag.isend());
}
