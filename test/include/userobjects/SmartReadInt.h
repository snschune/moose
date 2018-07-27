//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SmartReadInt_H
#define SmartReadInt_H

#include "GeneralUserObject.h"

class SmartReadInt;

template <>
InputParameters validParams<SmartReadInt>();

class SmartReadInt : public GeneralUserObject
{
public:
  SmartReadInt(const InputParameters & params);

  virtual void initialize(){};
  virtual void execute();
  virtual void finalize(){};

protected:
  int _smart_signed;
};

#endif /* SmartReadInt_H */
