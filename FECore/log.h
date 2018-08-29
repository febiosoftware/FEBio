#pragma once

#include "Logfile.h"
#include "FECoreKernel.h"

//! get the one-and-only log file (defined in LogFile.cpp)
//extern Logfile& felog;
#define felog (FECoreKernel::GetLogfile())
