//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// #define DEBUG_RAY_IF true
// #define DEBUG_RAY_MESH_IF true
// #define DEBUG_RAY_INTERSECTIONS

#include "MooseError.h"

#ifndef DEBUG_RAY_IF
#define debugRay(...) ((void)0)
#else
#define debugRay(...)                                                                              \
  if (DEBUG_RAY_IF)                                                                                \
  {                                                                                                \
    std::ostringstream ss;                                                                         \
    moose::internal::mooseStreamAll(ss, __VA_ARGS__);                                              \
    libMesh::err << "[" << _pid << "/" << _tid << "/" << (*_current_ray)->id() << "/" << __LINE__  \
                 << "]: " << ss.str() << std::endl                                                 \
                 << std::flush;                                                                    \
  }
#endif

#ifndef DEBUG_RAY_INTERSECTIONS
#define debugRaySimple(...) ((void)0)
#else
#define debugRaySimple(...)                                                                        \
  if (debug)                                                                                       \
  {                                                                                                \
    std::ostringstream ss;                                                                         \
    moose::internal::mooseStreamAll(ss, __VA_ARGS__);                                              \
    libMesh::err << ss.str() << std::endl << std::flush;                                           \
  }
#endif

#ifndef DEBUG_RAY_MESH_IF
#define possiblySaveDebugRayMesh() ((void)0)
#else
#define possiblySaveDebugRayMesh()                                                                 \
  if (DEBUG_RAY_MESH_IF && _debug_mesh)                                                            \
  {                                                                                                \
    _debug_mesh->prepare_for_use();                                                                \
    _debug_mesh->write("debug_ray_" + std::to_string((*_current_ray)->id()) + "_pc" +              \
                       std::to_string((*_current_ray)->processorCrossings()) + "_pid" +            \
                       std::to_string(_pid) + ".e");                                               \
    delete _debug_mesh;                                                                            \
  }
#endif

#ifndef DEBUG_RAY_MESH_IF
#define possiblyAddDebugRayMeshPoint(start, end) ((void)0)
#else
#define possiblyAddDebugRayMeshPoint(start, end)                                                   \
  if (_debug_mesh && DEBUG_RAY_MESH_IF)                                                            \
  {                                                                                                \
    _debug_mesh->add_point(start);                                                                 \
    _debug_node_count++;                                                                           \
    _debug_mesh->add_point(end);                                                                   \
    _debug_node_count++;                                                                           \
                                                                                                   \
    Elem * elem = Elem::build(EDGE2).release();                                                    \
    elem->subdomain_id() = 0;                                                                      \
                                                                                                   \
    elem = _debug_mesh->add_elem(elem);                                                            \
                                                                                                   \
    elem->set_node(0) = _debug_mesh->node_ptr((_debug_node_count - 2));                            \
    elem->set_node(1) = _debug_mesh->node_ptr((_debug_node_count - 2) + 1);                        \
  }
#endif
