#pragma once
#include <FEBioMix/FEMultiphasicStandard.h>

class FEBIOMIX_API FEFCDMaterialPoint : public FESolutesMaterialPoint
{
public:
    FEFCDMaterialPoint(FEMaterialPoint* ppt);
    
    void Init(bool bflag);
    
public:
    double    m_cFr;
};

class FEMultiphasicFCD : public FEMultiphasicStandard
{
public:
    FEMultiphasicFCD(FEModel* pfem) : FEMultiphasicStandard(pfem){}
    
    FEMaterialPoint* CreateMaterialPointData() override;
    
    //! fixed charge density
    double FixedChargeDensity(FEMaterialPoint& pt) override;
    
};
