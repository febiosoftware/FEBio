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
#include "FEItemList.h"
#include "DumpStream.h"
#include "FEModel.h"

FEItemList::FEItemList(FEModel* fem)// : FECoreBase(fem)
{
	if (fem)
	{
		m_mesh = &fem->GetMesh();
	}
}

FEItemList::FEItemList(FEMesh* mesh)// : FECoreBase(fem)
{
	m_mesh = mesh;
}

FEItemList::~FEItemList() {}

void FEItemList::Serialize(DumpStream& ar)
{
	ar & m_name;
}

// get the mesh
FEMesh* FEItemList::GetMesh() const
{
	return m_mesh;
}

void FEItemList::SetMesh(FEMesh* mesh)
{
	m_mesh = mesh;
}

const std::string& FEItemList::GetName() const
{
	return m_name;
}

void FEItemList::SetName(const std::string& name)
{
	m_name = name;
}

FEItemList* FEItemList::LoadClass(DumpStream& ar, FEItemList* p)
{
	assert(false);
	return nullptr;
}

void FEItemList::SaveClass(DumpStream& ar, FEItemList* p)
{
	assert(false);
}
