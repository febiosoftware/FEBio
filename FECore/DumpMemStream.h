#pragma once
#include "DumpStream.h"

//-----------------------------------------------------------------------------
//! The dump stream allows a class to record its internal state to a memory object
//! so that it can be restored later.
//! This can be used for storing the FEModel state during running restarts
class DumpMemStream : public DumpStream
{
public:
	DumpMemStream(FEModel& fem);
	~DumpMemStream();

public: // overloaded from base class
	size_t write(const void* pd, size_t size, size_t count);
	size_t read(void* pd, size_t size, size_t count);
	void clear();
	void Open(bool bsave, bool bshallow);
	void check();

protected:
	void grow_buffer(int l);
	void set_position(int l);

private:
	char*	m_pb;			//!< pointer to buffer
	char*	m_pd;			//!< position to insert a new value
	int		m_nsize;		//!< size of stream
	int		m_nreserved;	//!< size of reserved buffer
};
