/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
// units are defined as strings that use the following
// characters to represent the physical quantities that
// define the SI base units.
// (Note that the symbols for time and temperature differ
// from the conventional SI dimension symbols)
//
// L = length
// M = mass
// t = time
// T = temperature
// l = electric current
// n = substance
//
// In addition, the following symbols for derived units are also defined:
//
// F = force
// P = pressure
// E = energy
// W = power
// d = angles in degrees
// r = angles in radians
//
// or use one of these predefined constants
#define UNIT_NONE	""
#define UNIT_LENGTH "L"
#define UNIT_MASS   "M"
#define UNIT_TIME   "t"
#define UNIT_TEMPERATURE "T"
#define UNIT_CURRENT "I"
#define UNIT_SUBSTANCE	"n"
#define UNIT_FORCE "F"
#define UNIT_MOMENT "F.L"
#define UNIT_PRESSURE "P"
#define UNIT_ENERGY	"E"
#define UNIT_POWER "W"
#define UNIT_VOLTAGE "V"
#define UNIT_CONCENTRATION "c"
#define UNIT_RELATIVE_TEMPERATURE "R"

#define UNIT_DEGREE "d"
#define UNIT_RADIAN	"r"

#define UNIT_ACCELERATION "L/t^2"
#define UNIT_ANGULAR_MOMENTUM "E.t"
#define UNIT_ANGULAR_VELOCITY "r/t"
#define UNIT_AREA   "L^2"
#define UNIT_COUPLE_VISCOSITY "F.t/r"
#define UNIT_CURRENT_CONDUCTIVITY "I/V.L^2"
#define UNIT_CURRENT_DENSITY "I/L^2"
#define UNIT_DENSITY "M/L^3"
#define UNIT_DENSITY_RATE "M/L^3.t"
#define UNIT_DIFFUSIVITY  "L^2/t"
#define UNIT_ENERGY_AREAL_DENSITY "E/L^2"
#define UNIT_ENERGY_DENSITY "E/L^3"
#define UNIT_ENERGY_FLUX "W/L^2"
#define UNIT_FARADAY_CONSTANT "I.t/n"
#define UNIT_FILTRATION "L^2/F.t"
#define UNIT_FLOW_CAPACITANCE "L^5/F"
#define UNIT_FLOW_RATE "L^3/t"
#define UNIT_FLOW_RESISTANCE "F.t/L^3"
#define UNIT_GAS_CONSTANT "E/n.T"
#define UNIT_GRADIENT "1/L"
#define UNIT_LINEAR_MOMENTUM "F.t"
#define UNIT_MASS_FLOW_RATE "M/t"
#define UNIT_MOLAR_AREAL_CONCENTRATION "n/L^2"
#define UNIT_MOLAR_FLUX "n/L^2.t"
#define UNIT_MOLAR_MASS "M/n"
#define UNIT_MOLAR_VOLUME "L^3/n"
#define UNIT_PERMEABILITY    "L^4/F.t"
#define UNIT_POWER_DENSITY "W/L^3"
#define UNIT_RECIPROCAL_LENGTH "1/L"
#define UNIT_RECIPROCAL_TIME "1/t"
#define UNIT_ROTATIONAL_VISCOSITY "P.t/r"
#define UNIT_SPECIFIC_ENERGY "E/M"
#define UNIT_SPECIFIC_ENTROPY "E/M.T"
#define UNIT_SPECIFIC_FORCE "F/M"
#define UNIT_SPECIFIC_MOMENT "F.L/M"
#define UNIT_STIFFNESS "F/L"
#define UNIT_THERMAL_CONDUCTIVITY "W/L.T"
#define UNIT_VELOCITY "L/t"
#define UNIT_VISCOSITY "P.t"
#define UNIT_VOLUME "L^3"
