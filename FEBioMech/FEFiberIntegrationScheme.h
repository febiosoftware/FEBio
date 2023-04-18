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
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "febiomech_api.h"

//----------------------------------------------------------------------------------
// This is an iterator class that can be used to loop over all integration points of
// a fiber integration scheme. It must implement the Next functions.
// During the Next call, it should also update the fiber vector and weight.
// Next should return false if there is no more integration point
class FEFiberIntegrationSchemeIterator
{
public:
	FEFiberIntegrationSchemeIterator() {}
	virtual ~FEFiberIntegrationSchemeIterator() {}

	// Move to the next integration point
	// This also updates the m_fiber and m_weight members
	virtual bool Next() = 0;

	// check if the iterator is valid
	virtual bool IsValid() = 0;

public:
	vec3d	m_fiber;		// current fiber vector at integration point
	double	m_weight;		// current integration weight
};

//----------------------------------------------------------------------------------
// Base clase for integration schemes for continuous fiber distributions.
// The purpose of this class is mainly to provide an interface to the integration schemes
// for the FEBio input file. The code will use the GetIterator function to create an
// iterator that can be used to loop over all the integration points of the scheme and to
// evaluate the fiber vector and weights at each point.
class FEBIOMECH_API FEFiberIntegrationScheme : public FEMaterialProperty
{
public:
    FEFiberIntegrationScheme(FEModel* pfem);
    
	// Creates an iterator for the scheme. 
	// In general, the integration scheme may depend on the material point.
	// The passed material point pointer will be zero when evaluating the integrated fiber density
	virtual FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp = 0) = 0;

	FECORE_BASE_CLASS(FEFiberIntegrationScheme)
};
