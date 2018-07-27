//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StupidReadUnsignedInt.h"
#include <sstream>

registerMooseObject("MooseTestApp", StupidReadUnsignedInt);

template <>
InputParameters
validParams<StupidReadUnsignedInt>()
{
  InputParameters params = validParams<UserObject>();
  params.addRequiredRangeCheckedParam<unsigned int>("dump_unsigned", "dump_unsigned>0", "Number of elements in the angular direction");
  return params;
}

StupidReadUnsignedInt::StupidReadUnsignedInt(const InputParameters & params)
  : GeneralUserObject(params),
    _dumb_unsigned(getParam<unsigned int>("dump_unsigned"))
{

}

void
StupidReadUnsignedInt::execute()
{
  _console << _dumb_unsigned << std::endl;
}
