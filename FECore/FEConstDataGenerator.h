#pragma once
#include "FEDataGenerator.h"
#include "fecore_type.h"

template<class T> class FEConstDataGenerator : public FEDataGenerator
{
public:
	FEConstDataGenerator(FEModel* fem) : FEDataGenerator(fem), m_val(0.0) {}

	void value(const vec3d& r, T& d) override { d = m_val; }

	void BuildParamList() override
	{
		AddParameter(m_val, "value");
	}

private:
	T	m_val;
};
