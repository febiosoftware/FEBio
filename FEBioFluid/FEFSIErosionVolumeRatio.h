//
//  FEFSIErosionVolumeRatio.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/15/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFSIErosionVolumeRatio_hpp
#define FEFSIErosionVolumeRatio_hpp

#include <FECore/FEMeshAdaptor.h>

class FEFSIErosionVolumeRatio : public FEMeshAdaptor
{
public:
    FEFSIErosionVolumeRatio(FEModel* fem);
    
    bool Apply(int iteration) override;
    
private:
    double     m_minJ;
    int        m_maxElems;
    int        m_maxIters;
    int        m_metric;
    
    DECLARE_FECORE_CLASS()
};

#endif /* FEFSIErosionVolumeRatio_hpp */
