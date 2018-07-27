//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef StupidReadUnsignedInt_H
#define StupidReadUnsignedInt_H

#include "GeneralUserObject.h"

class StupidReadUnsignedInt;

template <>
InputParameters validParams<StupidReadUnsignedInt>();

class StupidReadUnsignedInt : public GeneralUserObject
{
public:
  StupidReadUnsignedInt(const InputParameters & params);

  virtual void initialize(){};
  virtual void execute();
  virtual void finalize(){};

protected:
  unsigned int _dumb_unsigned;
};

#endif /* StupidReadUnsignedInt_H */
