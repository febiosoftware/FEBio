// version.h : include file that contains the version numbers
//

#if !defined(AFX_VERSION_H__5901DABB_91FB_C34E_9011_12397479QBEE__INCLUDED_)
#define AFX_VERSION_H__5901DABB_91FB_C34E_9011_12397479QBEE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////
// FEBio version numbers
// VERSION is the main version number. This number is only incremented when
// major modifications or additions where added to the code
// SUBVERSION is only incremented when minor modifications or 
// additions where added to the code.
// SUBSUBVERSION is incremented when bugs are fixed.
//
// IMPORTANT NOTE: License files can only be used for FEBio versions 1.3.0 and up
//

#define VERSION			2
#define SUBVERSION		8
#define SUBSUBVERSION	0
#ifdef SVN
#include "svnrev.h"
#else
#define SVNREVISION 0
#endif
///////////////////////////////////////////////////////////////////////////////
// Restart file version
// This is the version number of the restart dump file format.
// It is incremented when the structure of this file is modified.
//

#define RSTRTVERSION		0x06

#endif // !defined(AFX_VERSION_H__5901DABB_91FB_C34E_9011_12397479QBEE__INCLUDED_)
