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
#include <vector>
#include "fecore_enum.h"
#include "fecore_api.h"

class FEElement;
class FEElementTraits;
class FEElementShape;

//-----------------------------------------------------------------------------
//! This class stores the different element traits classes

//! The purpose of this class is to store all the different element traits classes.
//! It are these traits classes that define the different element types. All different
//! traits must be registered before they can be assigned to elements.

class FECORE_API FEElementLibrary
{
public:
	//! destructor
	~FEElementLibrary();

	//! return the element library
	static FEElementLibrary* GetInstance();

	//! Assign a traits class to an element
	static void SetElementTraits(FEElement& el, int id);

	//! return element traits data
	static FEElementTraits* GetElementTraits(int ntype);

	//! return element shape class
	static FEElementShape* GetElementShapeClass(FE_Element_Shape eshape);

	//! return the element shape of a given element type
	static FE_Element_Shape GetElementShape(int ntype);

	//! return the element class of a given element type
	static FE_Element_Class GetElementClass(int ntype);

	//! checks if the element spec is valid
	static bool IsValid(const FE_Element_Spec& c);

	//! initialize library
	static void Initialize();

	//! get the element spec from the type
	static FE_Element_Spec GetElementSpecFromType(FE_Element_Type elemType);

private:
	//! constructor
	FEElementLibrary(){}
	FEElementLibrary(const FEElementLibrary&){}

	//! Function to register an element shape class
	int RegisterShape(FEElementShape* pshape);

	//! Function to register an element traits class
	int RegisterTraits(FEElementTraits* ptrait);

private:
	std::vector<FEElementTraits*>	m_Traits;	//!< pointer to registered element traits
	std::vector<FEElementShape*>	m_Shape;	//!< pointer to registered element shapes
	static FEElementLibrary*	m_pThis;
};
