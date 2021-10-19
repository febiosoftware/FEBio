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
#include "FEBioImport.h"
#include "FECore/FESurfacePairConstraint.h"
#include <map>

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEFileSection
{
public:
	FEBioBoundarySection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseBCFix         (XMLTag& tag);
	void ParseBCPrescribe   (XMLTag& tag);
	void ParseContactSection(XMLTag& tag);
	void ParseConstraints   (XMLTag& tag);
	void ParseSpringSection (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* psi);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);

protected:
	void BuildNodeSetMap();

	void AddFixedBC(FENodeSet* set, int bc);

protected:
	std::map<std::string, FENodeSet*>	m_NodeSet;	// map for faster lookup of node sets
};

//-----------------------------------------------------------------------------
// older formats (1.2)
class FEBioBoundarySection1x : public FEBioBoundarySection
{
public:
	FEBioBoundarySection1x(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// version 2.0
class FEBioBoundarySection2 : public FEBioBoundarySection
{
public:
	FEBioBoundarySection2(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix      (XMLTag& tag);
	void ParseBCPrescribe(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// version 2.5
class FEBioBoundarySection25 : public FEBioBoundarySection
{
public:
	FEBioBoundarySection25(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix      (XMLTag& tag);
	void ParseBCPrescribe(XMLTag& tag);
	void ParseBCRigid    (XMLTag& tag);
	void ParseRigidBody  (XMLTag& tag);
	void ParseBC         (XMLTag& tag);

	void ParsePeriodicLinearConstraint  (XMLTag& tag); // version 2.5 (temporary construction)
	void ParsePeriodicLinearConstraint2O(XMLTag& tag); // version 2.5 (temporary construction)
	void ParseMergeConstraint           (XMLTag& tag); // version 2.5
};
