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
#include "DumpStream.h"

//-----------------------------------------------------------------------------
//! The dump stream allows a class to record its internal state to a memory object
//! so that it can be restored later.
//! This can be used for storing the FEModel state during running restarts
class FECORE_API DumpMemStream : public DumpStream
{
public:
	DumpMemStream(FEModel& fem);
	~DumpMemStream();

public: // overloaded from base class
	size_t write(const void* pd, size_t size, size_t count);
	size_t read(void* pd, size_t size, size_t count);
	void clear();
	void Open(bool bsave, bool bshallow);

	size_t size() const { return m_nsize; }
	size_t reserved() const { return m_nreserved; }
	bool EndOfStream() const;

protected:
	void grow_buffer(size_t l);
	void set_position(size_t l);

private:
	char*	m_pb;			//!< pointer to buffer
	char*	m_pd;			//!< position to insert a new value
	size_t	m_nsize;		//!< size of stream
	size_t	m_nreserved;	//!< size of reserved buffer
};
