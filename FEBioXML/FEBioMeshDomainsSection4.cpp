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
#include "FEBioMeshDomainsSection4.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEShellDomain.h>
#include <FECore/FETrussDomain.h>
#include <FECore/FEDomain2D.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FENodeNodeList.h>
#include <sstream>

FEBioMeshDomainsSection4::FEBioMeshDomainsSection4(FEBioImport* pim) : FEBioFileSection(pim) 
{
	m_noff = 0;
}

void FEBioMeshDomainsSection4::Parse(XMLTag& tag)
{
	// build the node ID lookup table
	BuildNLT();
	
	// read all sections
	if (tag.isleaf() == false)
	{
		++tag;
		do
		{
			if      (tag == "SolidDomain") ParseSolidDomainSection(tag);
			else if (tag == "ShellDomain") ParseShellDomainSection(tag);
			else if (tag == "BeamDomain" ) ParseBeamDomainSection(tag);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}

	// let's build the part
	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (feb.BuildPart(*GetFEModel(), *part, false) == false)
	{
		throw XMLReader::Error("Failed building parts.");
	}

	// tell the file reader to rebuild the node ID table
	GetBuilder()->BuildNodeList();

	// At this point the mesh is completely read in.
	// allocate material point data
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
	}

	// Now we can allocate the degrees of freedom.
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

void FEBioMeshDomainsSection4::BuildNLT()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) return;

	// build node-index lookup table
	int noff = -1, maxID = 0;
	int N0 = mesh.Nodes();
	int NN = part->Nodes();
	for (int i = 0; i < NN; ++i)
	{
		int nid = part->GetNode(i).id;
		if ((noff < 0) || (nid < noff)) noff = nid;
		if (nid > maxID) maxID = nid;
	}
	m_NLT.assign(maxID - noff + 1, -1);
	for (int i = 0; i < NN; ++i)
	{
		int nid = part->GetNode(i).id - noff;
		m_NLT[nid] = i + N0;
	}

	m_noff = noff;
}

void FEBioMeshDomainsSection4::ParseSolidDomainSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) throw XMLReader::InvalidTag(tag);

	// get the domain name
	const char* szname = tag.AttributeValue("name");
	FEBModel::Domain* partDomain = part->FindDomain(szname);
	if (partDomain == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name", szname);

	// get the material name
	const char* szmat = tag.AttributeValue("mat");
	FEMaterial* mat = fem.FindMaterial(szmat);
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "mat", szmat);

	// set the material name
	partDomain->SetMaterialName(szmat);

	// see if the element type is specified
	const char* szelem = tag.AttributeValue("elem_type", true);
	if (szelem)
	{
		FE_Element_Spec elemSpec = partDomain->ElementSpec();
		FE_Element_Spec newSpec = GetBuilder()->ElementSpec(szelem);

		// make sure it's valid
		if ((FEElementLibrary::IsValid(newSpec) == false) || (elemSpec.eshape != newSpec.eshape))
		{
			throw XMLReader::InvalidAttributeValue(tag, "elem_type", szelem);
		}

		partDomain->SetElementSpec(newSpec);
	}

	// get the (optional) type attribute
	const char* sztype = tag.AttributeValue("type", true);

	// --- build the domain --- 
	// we'll need the kernel for creating domains
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// element count
	int elems = partDomain->Elements();

	// get the element spect
	FE_Element_Spec spec = partDomain->ElementSpec();

	// create the domain
	FEDomain* dom = nullptr;
	if (sztype)
	{
		// if the type attribute is defined, try to allocate the domain class directly. 
		dom = febio.CreateDomainExplicit(FESOLIDDOMAIN_ID, sztype, &fem);
		if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, sztype);
		dom->SetMaterial(mat);
	}
	else
	{
		// if not, then use "old" logic, which tries to match the domain using the 
		// domain factories. 
		dom = febio.CreateDomain(spec, &mesh, mat);
		if (dom == 0) throw XMLReader::InvalidTag(tag);
	}

	// add it to the mesh
	mesh.AddDomain(dom);

	// Allocate elements
	if (dom->Create(elems, spec) == false)
	{
		throw XMLReader::InvalidTag(tag);
	}

	// assign the material
	dom->SetMatID(mat->GetID() - 1);

	// get the part name
	string partName = part->Name();
	if (partName.empty() == false) partName += ".";

	string domName = partName + partDomain->Name();
	dom->SetName(domName);

	// read additional parameters
	if (tag.isleaf() == false)
	{
		ReadParameterList(tag, dom);
	}

	// process element data
	for (int j = 0; j < elems; ++j)
	{
		const FEBModel::ELEMENT& domElement = partDomain->GetElement(j);

		FEElement& el = dom->ElementRef(j);
		el.SetID(domElement.id);

		int ne = el.Nodes();
		for (int n = 0; n < ne; ++n) el.m_node[n] = m_NLT[domElement.node[n] - m_noff];
	}
}

void FEBioMeshDomainsSection4::ParseShellDomainSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) throw XMLReader::InvalidTag(tag);

	// get the domain name
	const char* szname = tag.AttributeValue("name");
	FEBModel::Domain* partDomain = part->FindDomain(szname);
	if (partDomain == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name", szname);

	// get the material name
	const char* szmat = tag.AttributeValue("mat");
	FEMaterial* mat = fem.FindMaterial(szmat);
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "mat", szmat);

	// set the material name
	partDomain->SetMaterialName(szmat);

	// see if the element type is specified
	const char* szelem = tag.AttributeValue("elem_type", true);
	if (szelem)
	{
		FE_Element_Spec elemSpec = partDomain->ElementSpec();
		FE_Element_Spec newSpec = GetBuilder()->ElementSpec(szelem);

		// make sure it's valid
		if ((FEElementLibrary::IsValid(newSpec) == false) || (elemSpec.eshape != newSpec.eshape))
		{
			throw XMLReader::InvalidAttributeValue(tag, "elem_type", szelem);
		}

		partDomain->SetElementSpec(newSpec);
	}

	// get the (optional) type attribute
	const char* sztype = tag.AttributeValue("type", true);

	// --- build the domain --- 
	// we'll need the kernel for creating domains
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// element count
	int elems = partDomain->Elements();

	// get the element spect
	FE_Element_Spec spec = partDomain->ElementSpec();

	// create the domain
	FEDomain* dom = nullptr;
	if (sztype)
	{
		// if the type attribute is defined, try to allocate the domain class directly. 
		dom = febio.CreateDomainExplicit(FESHELLDOMAIN_ID, sztype, &fem);
		if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, sztype);
		dom->SetMaterial(mat);
	}
	else
	{
		// if not, then use "old" logic, which tries to match the domain using the 
		// domain factories. 
		dom = febio.CreateDomain(spec, &mesh, mat);
		if (dom == 0) throw XMLReader::InvalidTag(tag);
	}

	mesh.AddDomain(dom);

	if (dom->Create(elems, spec) == false)
	{
		throw XMLReader::InvalidTag(tag);
	}

	// assign the material
	dom->SetMatID(mat->GetID() - 1);

	// get the part name
	string partName = part->Name();
	if (partName.empty() == false) partName += ".";

	string domName = partName + partDomain->Name();
	dom->SetName(domName);

	// read additional parameters
	if (tag.isleaf() == false)
	{
		ReadParameterList(tag, dom);
	}

	// process element data
	FEShellDomainNew* shellDomain = dynamic_cast<FEShellDomainNew*>(dom);
	if (shellDomain)
	{
		for (int j = 0; j < elems; ++j)
		{
			const FEBModel::ELEMENT& domElement = partDomain->GetElement(j);

			FEShellElement& el = shellDomain->Element(j);
			el.SetID(domElement.id);

			int ne = el.Nodes();
			for (int n = 0; n < ne; ++n)
			{
				el.m_node[n] = m_NLT[domElement.node[n] - m_noff];
				el.m_h0[n] = 0.0;
			}
		}

		shellDomain->AssignDefaultShellThickness();
	}
}

void FEBioMeshDomainsSection4::ParseBeamDomainSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) throw XMLReader::InvalidTag(tag);

	// get the domain name
	const char* szname = tag.AttributeValue("name");
	FEBModel::Domain* partDomain = part->FindDomain(szname);
	if (partDomain == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name", szname);

	// get the material name
	const char* szmat = tag.AttributeValue("mat");
	FEMaterial* mat = fem.FindMaterial(szmat);
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "mat", szmat);

	// set the material name
	partDomain->SetMaterialName(szmat);

	// see if the element type is specified
	const char* szelem = tag.AttributeValue("elem_type", true);
	if (szelem)
	{
		FE_Element_Spec elemSpec = partDomain->ElementSpec();
		FE_Element_Spec newSpec = GetBuilder()->ElementSpec(szelem);

		// make sure it's valid
		if ((FEElementLibrary::IsValid(newSpec) == false) || (elemSpec.eshape != newSpec.eshape))
		{
			throw XMLReader::InvalidAttributeValue(tag, "elem_type", szelem);
		}

		partDomain->SetElementSpec(newSpec);
	}

	// get the (optional) type attribute
	const char* sztype = tag.AttributeValue("type", true);

	// --- build the domain --- 
	// we'll need the kernel for creating domains
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// element count
	int elems = partDomain->Elements();

	// get the element spect
	FE_Element_Spec spec = partDomain->ElementSpec();

	// create the domain
	FEDomain* dom = nullptr;
	if (sztype)
	{
		// if the type attribute is defined, try to allocate the domain class directly. 
		dom = febio.CreateDomainExplicit(FEBEAMDOMAIN_ID, sztype, &fem);
		if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, sztype);
		dom->SetMaterial(mat);
	}
	else
	{
		// if not, then use "old" logic, which tries to match the domain using the 
		// domain factories. 
		dom = febio.CreateDomain(spec, &mesh, mat);
		if (dom == 0) throw XMLReader::InvalidTag(tag);
	}

	// add it to the mesh
	mesh.AddDomain(dom);

	// Allocate elements
	if (dom->Create(elems, spec) == false)
	{
		throw XMLReader::InvalidTag(tag);
	}

	// assign the material
	dom->SetMatID(mat->GetID() - 1);

	// get the part name
	string partName = part->Name();
	if (partName.empty() == false) partName += ".";

	string domName = partName + partDomain->Name();
	dom->SetName(domName);

	// read additional parameters
	if (tag.isleaf() == false)
	{
		ReadParameterList(tag, dom);
	}

	// process element data
	for (int j = 0; j < elems; ++j)
	{
		const FEBModel::ELEMENT& domElement = partDomain->GetElement(j);

		FEElement& el = dom->ElementRef(j);
		el.SetID(domElement.id);

		// TODO: This assumes one-based indexing of all nodes!
		int ne = el.Nodes();
		for (int n = 0; n < ne; ++n) el.m_node[n] = m_NLT[domElement.node[n] - m_noff];
	}
}
