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
#include "FEElasticMaterial2O.h"
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>
#include "FE2OMicroConstraint.h"
#include "FEMicroMaterial.h"
#include "FERVEModel2O.h"

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint2O : public FEMaterialPointData
{
public:
	//! constructor
	FEMicroMaterialPoint2O(FEMaterialPointData* mp);

	//! create a shallow copy
	FEMaterialPointData* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	FEMicroModel2O m_rve;				//!< local copy of the rve		
	int		m_elem_id;		//!< element ID
	int		m_gpt_id;		//!< Gauss point index (0-based)
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial2O : public FEElasticMaterial2O
{
public:
	FEMicroMaterial2O(FEModel* pfem);
	~FEMicroMaterial2O(void);

public:
	std::string		m_szrve;		//!< filename for RVE file
	std::string		m_szbc;			//!< name of nodeset defining boundary
	int				m_rveType;		//!< RVE type
	double			m_scale;		//!< geometry scale factor
	FERVEModel2O	m_mrve;			//!< the parent RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	void Stress(FEMaterialPoint& mp, mat3d& P, tens3drs& Q) override;

	//! calculate tangent stiffness at material point
	void Tangent(FEMaterialPoint &mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J) override;
	
	//! data initialization
	bool Init() override;

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() override;

public:
	int Probes() { return (int) m_probe.size(); }
	FERVEProbe& Probe(int i) { return *m_probe[i]; }

protected:
	std::vector<FERVEProbe*>	m_probe;

public:
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
