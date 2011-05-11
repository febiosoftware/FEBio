#pragma once

class FEModel;

template <typename T> class FEBioFactory_T
{
public:
	virtual T* Create(FEModel* pfem) = 0;
};
