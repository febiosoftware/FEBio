#include "stdafx.h"
#include "FECallBack.h"
#include "FEModel.h"

REGISTER_SUPER_CLASS(FECallBack, FECALLBACK_ID);

bool do_FECallBack_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FECallBack* pCB = (FECallBack*) (pd);
	return pCB->Execute(*pfem, nwhen);
}

FECallBack::FECallBack(FEModel* fem, int when) : FECoreBase(fem, FECALLBACK_ID)
{
	fem->AddCallback(do_FECallBack_cb, when, this);
}

FECallBack::~FECallBack()
{
}
