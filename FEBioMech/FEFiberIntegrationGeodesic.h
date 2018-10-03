//
//  FEFiberIntegrationGeodesic.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/19/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#pragma once
#include "FEFiberIntegrationScheme.h"
#include "geodesic.h"

//----------------------------------------------------------------------------------
// Geodesic dome integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGeodesic : public FEFiberIntegrationScheme
{
	class Iterator;

public:
    FEFiberIntegrationGeodesic(FEModel* pfem);
    ~FEFiberIntegrationGeodesic();
	
	//! Initialization
	bool Init() override;
    
	// serialization
	void Serialize(DumpStream& ar) override;

	// get iterator
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp) override;

protected:
	void InitIntegrationRule();  

private: // parameters
	int             m_nres;	// resolution

protected:
    int             m_nint; // number of integration points
	double          m_cth[NSTH];
	double          m_sth[NSTH];
	double          m_cph[NSTH];
	double          m_sph[NSTH];
	double          m_w[NSTH];

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
