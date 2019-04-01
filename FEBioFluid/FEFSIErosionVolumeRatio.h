#pragma once
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
