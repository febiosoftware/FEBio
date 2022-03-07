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
#include "FEMeshPartition.h"

// forward declaration of material class
class FEMaterial;

// Base class for solid and shell parts. Domains can also have materials assigned.
class FECORE_API FEDomain : public FEMeshPartition
{
public:
	FEDomain(int nclass, FEModel* fem);

	//! get the material of this domain
	virtual FEMaterial* GetMaterial() { return 0; }

	// assign a material to this domain
	virtual void SetMaterial(FEMaterial* pm);

	//! set the material ID of all elements
	void SetMatID(int mid);

	//! Allocate material point data for the elements
	//! This is called after elements get read in from the input file.
	//! And must be called before material point data can be accessed.
	//! \todo Perhaps I can make this part of the "creation" routine
	void CreateMaterialPointData();

	// serialization
	void Serialize(DumpStream& ar) override;

	//! augmentation
	// NOTE: This is here so that the FESolver can do the augmentations
	// for the 3-field hex/shell domains.
	virtual bool Augment(int naug) { return true; }

	// create function
	virtual bool Create(int elements, FE_Element_Spec espec) = 0;

public:
	//! Get the list of dofs on this domain
	virtual const FEDofList& GetDOFList() const = 0;

	//! Unpack the LM data for an element of this domain
	virtual void UnpackLM(FEElement& el, vector<int>& lm);

	//! build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M);

	//! Activate the domain
	virtual void Activate();

protected:
	// helper function for activating dof lists
	void Activate(const FEDofList& dof);

	// helper function for unpacking element dofs
	void UnpackLM(FEElement& el, const FEDofList& dof, vector<int>& lm);
};
