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
#include "FEModelComponent.h"
#include "FEModel.h"
#include "DumpStream.h"
#include <string.h>

//-----------------------------------------------------------------------------
FEModelComponent::FEModelComponent(FEModel* fem) : FECoreBase(fem)
{
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{
	
}

//-----------------------------------------------------------------------------
void FEModelComponent::Update()
{

}

//-----------------------------------------------------------------------------
double FEModelComponent::CurrentTime() const
{
	return GetFEModel()->GetTime().currentTime;
}

//-----------------------------------------------------------------------------
double FEModelComponent::CurrentTimeIncrement() const
{
	return GetFEModel()->GetTime().timeIncrement;
}

//-----------------------------------------------------------------------------
double FEModelComponent::GetGlobalConstant(const char* sz) const
{
	return GetFEModel()->GetGlobalConstant(sz);
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetDOFIndex(const char* szvar, int n) const
{
	return GetFEModel()->GetDOFIndex(szvar, n);
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetDOFIndex(const char* szdof) const
{
	return GetFEModel()->GetDOFIndex(szdof);
}

//-----------------------------------------------------------------------------
//! Get the model's mesh
FEMesh& FEModelComponent::GetMesh()
{
	return GetFEModel()->GetMesh();
}

//-----------------------------------------------------------------------------
const FETimeInfo& FEModelComponent::GetTimeInfo() const
{
	return GetFEModel()->GetTime();
}

//-----------------------------------------------------------------------------
void FEModelComponent::AttachLoadController(const char* szparam, int lc)
{
	FEParam* p = GetParameter(szparam); assert(p);
	if (p)
	{
		FEModel* fem = GetFEModel();
		fem->AttachLoadController(p, lc);
	}
}

//-----------------------------------------------------------------------------
void FEModelComponent::AttachLoadController(void* pd, int lc)
{
	FEParam* p = FindParameterFromData(pd); assert(p);
	if (p)
	{
		FEModel* fem = GetFEModel();
		fem->AttachLoadController(p, lc);
	}
}
