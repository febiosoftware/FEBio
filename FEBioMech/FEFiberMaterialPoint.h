//
//  FEFiberMaterialPoint.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberMaterialPoint__
#define __FEBioMech__FEFiberMaterialPoint__

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for elastic fiber materials
class FEFiberMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEFiberMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
	//! copy material point data
	FEMaterialPoint* Copy();
    
	//! Initialize material point data
	void Init();
    
	//! Serialize data to archive
	void Serialize(DumpStream& ar);
    
public:
    vec3d   m_n0;    //!< fiber direction, reference configuration, global CS
};


#endif /* defined(__FEBioMech__FEFiberMaterialPoint__) */
