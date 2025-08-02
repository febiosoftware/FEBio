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
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Geodesic dome integration scheme for continuous fiber distributions
//
class FEFiberIntegrationTriangle : public FEFiberIntegrationScheme
{
	class Iterator;

public:
    FEFiberIntegrationTriangle(FEModel* pfem);
    ~FEFiberIntegrationTriangle();
	
	//! Initialization
	bool Init() override;
    
	// serialization
	void Serialize(DumpStream& ar) override;

	// create iterator
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp) override;

protected:
	void InitIntegrationRule();
    
public: // parameters
	int             m_nres;	// resolution
    int             m_nre;  // resolution entry

protected:
    int             m_nint; // number of integration points
	double          m_cth[2000];
	double          m_sth[2000];
	double          m_cph[2000];
	double          m_sph[2000];
	double          m_w[2000];
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
