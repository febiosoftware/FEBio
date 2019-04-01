#pragma once
#include <vector>
#include "fecore_api.h"

//-----------------------------------------------------------------------------
enum IOResult { IO_ERROR, IO_OK, IO_END };

//-----------------------------------------------------------------------------
// Output archive
class FECORE_API Archive
{
public:
	Archive();
	virtual ~Archive();

	// virtual function for writing data to an archive
	virtual void WriteData(int nid, std::vector<float>& data);
};
