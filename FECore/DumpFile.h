/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



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
	size_t write(const void* pd, size_t size, size_t count) override;
	size_t read(void* pd, size_t size, size_t count) override;
	void clear() override {}
	bool EndOfStream() const override;

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

	size_t Size() { return m_size; }

protected:
	FILE*		m_fp;		//!< The actual file pointer
	size_t		m_size;
};
