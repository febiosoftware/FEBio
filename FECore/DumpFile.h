#pragma once

#include <stdio.h>
#include "DumpStream.h"

//-----------------------------------------------------------------------------
//! Class for serializing data to a binary archive.

//! This class is used to read data from or write
//! data to a binary file. The class defines several operators to 
//! simplify in- and output.
//! \sa FEM::Serialize()

class FECORE_API DumpFile : public DumpStream
{
public:
	// overloaded from DumpStream
	size_t write(const void* pd, size_t size, size_t count);
	size_t read(void* pd, size_t size, size_t count);
	void clear(){}

public:
	DumpFile(FEModel& fem);
	virtual ~DumpFile();

	//! Open archive for reading
	bool Open(const char* szfile);

	//! Open archive for writing
	bool Create(const char* szfile);

	//! Open archive for appending
	bool Append(const char* szfile);

	//! Close archive
	void Close();

	//! See if the archive is valid
	bool IsValid() { return (m_fp != 0); }

	//! Flush the archive
	void Flush() { fflush(m_fp); }

	//! get the current index
	int GetDataIndex() const { return m_nindex; }

protected:
	FILE*		m_fp;		//!< The actual file pointer
	int			m_nindex;	//!< file index (gives amount of bytes written or read in so far)
};
