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
#include "FEModelLoad.h"
#include "FESurface.h"
#include "FESolver.h"
#include "FETimeInfo.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FECORE_API FESurfaceLoad : public FEModelLoad
{
public:
	FESurfaceLoad(FEModel* pfem);
	virtual ~FESurfaceLoad(void);

	//! Set the surface to apply the load to
	virtual void SetSurface(FESurface* ps);

	bool Init() override;

	//! Get the surface
	FESurface& GetSurface() { return *m_psurf; }

	void Serialize(DumpStream& ar) override;

public:
    virtual double ScalarLoad(FESurfaceMaterialPoint& mp) { return 0; }
    virtual vec3d VectorLoad(FESurfaceMaterialPoint& mp) { return vec3d(0,0,0); }

	// TODO: Can I get rid of this?
	// This is needed to update the mesh after some surface loads, which aren't really
	// surface loads modify boundary conditions
	void ForceMeshUpdate();

protected:
	FESurface*	m_psurf;

	FECORE_BASE_CLASS(FESurfaceLoad)
};
