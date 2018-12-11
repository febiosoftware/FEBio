#include "stdafx.h"
#include "FEDomainParameter.h"

FEDomainParameter::FEDomainParameter(const std::string& name) : m_name(name)
{

}

FEDomainParameter::~FEDomainParameter()
{

}

void FEDomainParameter::setName(const std::string& name)
{
	m_name = name;
}

const std::string& FEDomainParameter::name() const
{
	return m_name;
}
