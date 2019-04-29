/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include <string.h>
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEModelComponent::FEModelComponent(FEModel* fem) : FECoreBase(fem)
{
	// the ID can be used by derived class to define a identifier for derived classes
	// This value needs to be set in the constructor
	m_nID = 0;

	// initialize parameters
	m_bactive = true;
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{
	
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetID() const
{
	return m_nID;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetID(int n)
{
	m_nID = n;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::IsActive() const
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Deactivate()
{ 
	m_bactive = false; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Update()
{

}

//-----------------------------------------------------------------------------
void FEModelComponent::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_nID;
	ar & m_bactive;
}
