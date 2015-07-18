#include "stdafx.h"
#include "FECallBack.h"
#include "FEModel.h"

void do_FECallBack_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FECallBack* pCB = (FECallBack*) (pd);
	pCB->Execute(*pfem, nwhen);
}

FECallBack::FECallBack(FEModel* pfem, int when)
{
	pfem->AddCallback(do_FECallBack_cb, when, this);
}

FECallBack::~FECallBack()
{
}
