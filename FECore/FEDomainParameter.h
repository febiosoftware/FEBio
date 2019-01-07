#pragma once
#include "FEParam.h"
#include "FEMaterialPoint.h"

//! The domain parameter is a mechanism for accessing material point data
//! indirectly through the domain. 
//! Domains keep lists of domain parameters that can be queried.
//! Notice that these classes return FEParamValue so that they can be used in optimization
class FECORE_API FEDomainParameter
{
public:
	FEDomainParameter(const std::string& name);
	virtual ~FEDomainParameter();

	void setName(const std::string& name);
	const std::string& name() const;

	//! derived classes must override this function
	virtual FEParamValue value(FEMaterialPoint& mp) = 0;

private:
	std::string	m_name;
};
