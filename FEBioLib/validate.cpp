#include "stdafx.h"
#include "febio.h"
#include "validate.h"

#ifndef FEBIOLM
int GetLicenseKeyStatus()
{
	return 0;
}
#else
#include <febiolm.h>

#endif
