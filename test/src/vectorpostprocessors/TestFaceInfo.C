//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TestFaceInfo.h"

// MOOSE includes
#include "MooseMesh.h"

registerMooseObject("MooseApp", TestFaceInfo);

defineLegacyParams(TestFaceInfo);

InputParameters
TestFaceInfo::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  return params;
}

TestFaceInfo::TestFaceInfo(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _face_info(_fe_problem.mesh().faceInfo()),
    _face_id(declareVector("id")),
    _face_area(declareVector("area")),
    _left_element_id(declareVector("left_elem")),
    _right_element_id(declareVector("right_elem")),
    _left_element_side(declareVector("left_side")),
    _right_element_side(declareVector("right_side")),
    _nx(declareVector("nx")),
    _ny(declareVector("ny")),
    _nz(declareVector("nz")),
    _face_cx(declareVector("face_cx")),
    _face_cy(declareVector("face_cy")),
    _face_cz(declareVector("face_cz")),
    _left_cx(declareVector("left_cx")),
    _left_cy(declareVector("left_cy")),
    _left_cz(declareVector("left_cz")),
    _right_cx(declareVector("right_cx")),
    _right_cy(declareVector("right_cy")),
    _right_cz(declareVector("right_cz"))
{
}

void
TestFaceInfo::execute()
{
  for (unsigned int j = 0; j < _face_info.nFaces(); ++j)
  {
    _face_id.push_back(j);
    _face_area.push_back(_face_info.area(j));
    _left_element_id.push_back(_face_info.leftElem(j)->id());
    _right_element_id.push_back(_face_info.rightElem(j)->id());
    _left_element_side.push_back(_face_info.leftSideID(j));
    _right_element_side.push_back(_face_info.rightSideID(j));

    Point normal = _face_info.normal(j);
    _nx.push_back(normal(0));
    _ny.push_back(normal(1));
    _nz.push_back(normal(2));
    Point fc = _face_info.faceCentroid(j);
    _face_cx.push_back(fc(0));
    _face_cy.push_back(fc(1));
    _face_cz.push_back(fc(2));
    Point lc = _face_info.leftCentroid(j);
    _left_cx.push_back(lc(0));
    _left_cy.push_back(lc(1));
    _left_cz.push_back(lc(2));
    Point rc = _face_info.rightCentroid(j);
    _right_cx.push_back(rc(0));
    _right_cy.push_back(rc(1));
    _right_cz.push_back(rc(2));

    if (_face_info.leftElem(j) != _face_info.elements(j).first)
      mooseError("Mismatch of elements and leftElem functions");
    if (_face_info.rightElem(j) != _face_info.elements(j).second)
      mooseError("Mismatch of elements and leftElem functions");
    if (_face_info.leftSideID(j) != _face_info.sideIDs(j).first)
      mooseError("Mismatch of sideIDs and leftSideID functions");
    if (_face_info.rightSideID(j) != _face_info.sideIDs(j).second)
      mooseError("Mismatch of sideIDs and rightSideID functions");
    if (_face_info.centroids(j).first != lc)
      mooseError("Mismatch of centroids and leftCentroid functions");
    if (_face_info.centroids(j).second != rc)
      mooseError("Mismatch of centroids and rightCentroid functions");
  }
}
