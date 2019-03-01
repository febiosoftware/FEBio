#pragma once
#include <string>
#include "fecore_api.h"

class DumpStream;

class FECORE_API FEItemList
{
public:
	FEItemList() {}
	virtual ~FEItemList() {}

	void SetName(const std::string& name) { m_name = name;  }
	const std::string& GetName() const { return m_name; }

	void Serialize(DumpStream& ar);

private:
	std::string	m_name;
};
