// Archive.h: interface for the Archive class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <vector>

//-----------------------------------------------------------------------------
enum IOResult { IO_ERROR, IO_OK, IO_END };

//-----------------------------------------------------------------------------
// Output archive
class Archive  
{
public:
	Archive();
	virtual ~Archive();

	// virtual function for writing data to an archive
	virtual void WriteData(int nid, std::vector<float>& data);
};
