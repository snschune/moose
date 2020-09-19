//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PatchSidesetGenerator.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "CastUniquePointer.h"
#include "MooseUtils.h"
#include "MooseMeshUtils.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/linear_partitioner.h"
#include "libmesh/centroid_partitioner.h"
#include "libmesh/parmetis_partitioner.h"
#include "libmesh/hilbert_sfc_partitioner.h"
#include "libmesh/morton_sfc_partitioner.h"
#include "libmesh/enum_elem_type.h"

// libmesh elem types
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad9.h"

#include <set>
#include <limits>

registerMooseObject("HeatConductionApp", PatchSidesetGenerator);

InputParameters
PatchSidesetGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<boundary_id_type>("sideset",
                                            "The sideset that will be divided into patches");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "n_patches", "n_patches>0", "Number of patches");

  MooseEnum partitioning(
      "default=-3 metis=-2 parmetis=-1 linear=0 grid=2 centroid hilbert_sfc morton_sfc", "default");
  params.addParam<MooseEnum>(
      "partitioner",
      partitioning,
      "Specifies a mesh partitioner to use when splitting the mesh for a parallel computation.");
  MooseEnum direction("x y z radial");
  params.addParam<MooseEnum>("centroid_partitioner_direction",
                             direction,
                             "Specifies the sort direction if using the centroid partitioner. "
                             "Available options: x, y, z, radial");

  params.addParamNamesToGroup("partitioner centroid_partitioner_direction", "Partitioning");

  params.addClassDescription(
      "Divides the given sideset into smaller patches of roughly equal size.");

  return params;
}

PatchSidesetGenerator::PatchSidesetGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _n_patches(getParam<unsigned int>("n_patches")),
    _sideset(getParam<boundary_id_type>("sideset")),
    _partitioner_name(getParam<MooseEnum>("partitioner"))
{
}

std::unique_ptr<MeshBase>
PatchSidesetGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  _mesh->errorIfDistributedMesh("PatchSidesetGenerator");

  // Get a reference to our BoundaryInfo object for later use
  BoundaryInfo & boundary_info = mesh->get_boundary_info();

  // get dimensionality
  _dim = mesh->mesh_dimension() - 1;

  // get a list of all sides; vector of tuples (elem, loc_side, side_set)
  auto side_list = boundary_info.build_active_side_list();

  // create a dim - 1 dimensional mesh
  auto boundary_mesh =
      libmesh_make_unique<libMesh::ReplicatedMesh>(comm(), mesh->mesh_dimension() - 1);
  boundary_mesh->set_mesh_dimension(mesh->mesh_dimension() - 1);
  boundary_mesh->set_spatial_dimension(mesh->mesh_dimension());

  // nodes in the new mesh by boundary_node_id (index)
  std::vector<Node *> boundary_nodes;
  // a map from the node numbering on the volumetric mesh to the numbering
  // on the boundary_mesh
  std::map<dof_id_type, dof_id_type> mesh_node_id_to_boundary_node_id;
  // a local counter keeping track of how many entries have been added to boundary_nodes
  dof_id_type boundary_node_id = 0;
  // a map from new element id in the boundary mesh to the element id/side/sideset
  // tuple it came from
  std::map<dof_id_type, std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
      boundary_elem_to_mesh_elem;
  for (auto & side : side_list)
  {
    if (std::get<2>(side) == _sideset)
    {
      // the original volumetric mesh element
      const Elem * elem = mesh->elem_ptr(std::get<0>(side));

      // the boundary element
      std::unique_ptr<const Elem> boundary_elem = elem->side_ptr(std::get<1>(side));

      // an array that saves the boundary node ids of this elem in the right order
      std::vector<dof_id_type> bnd_elem_node_ids(boundary_elem->n_nodes());

      // loop through the nodes in boundary_elem
      for (unsigned int j = 0; j < boundary_elem->n_nodes(); ++j)
      {
        const Node * node = boundary_elem->node_ptr(j);

        // Is this node a new node?
        if (mesh_node_id_to_boundary_node_id.find(node->id()) ==
            mesh_node_id_to_boundary_node_id.end())
        {
          // yes, it is new, need to add it to the mesh_node_id_to_boundary_node_id map
          mesh_node_id_to_boundary_node_id.insert(
              std::pair<dof_id_type, dof_id_type>(node->id(), boundary_node_id));

          // this adds this node to the boundary mesh and puts it at the right position
          // in the boundary_nodes array
          Point pt(*node);
          boundary_nodes.push_back(boundary_mesh->add_point(pt, boundary_node_id));

          // keep track of the boundary node for setting up the element
          bnd_elem_node_ids[j] = boundary_node_id;

          // increment the boundary_node_id counter
          ++boundary_node_id;
        }
        else
          bnd_elem_node_ids[j] = mesh_node_id_to_boundary_node_id.find(node->id())->second;
      }

      // all nodes for this element have been added, so we can add the element to the
      // boundary mesh
      Elem * new_bnd_elem = boundaryElementHelper(*boundary_mesh, boundary_elem->type());

      // keep track of these new boundary elements in boundary_elem_to_mesh_elem
      boundary_elem_to_mesh_elem.insert(
          std::pair<dof_id_type, std::tuple<dof_id_type, unsigned short int, boundary_id_type>>(
              new_bnd_elem->id(), side));

      // set the nodes & subdomain_id of the new element by looping over the
      // boundary_elem and then inserting its nodes into new_bnd_elem in the
      // same order
      for (unsigned int j = 0; j < boundary_elem->n_nodes(); ++j)
      {
        dof_id_type old_node_id = boundary_elem->node_ptr(j)->id();
        if (mesh_node_id_to_boundary_node_id.find(old_node_id) ==
            mesh_node_id_to_boundary_node_id.end())
          mooseError("Node id", old_node_id, " not linked to new node id.");
        dof_id_type new_node_id = mesh_node_id_to_boundary_node_id.find(old_node_id)->second;
        new_bnd_elem->set_node(j) = boundary_nodes[new_node_id];
      }
    }
  }

  // partition the boundary mesh
  boundary_mesh->prepare_for_use();
  _n_boundary_mesh_elems = boundary_mesh->n_elem();
  if (_partitioner_name == "grid")
    partition(*boundary_mesh);
  else
  {
    MooseMesh::setPartitioner(*boundary_mesh, _partitioner_name, false, _pars, *this);
    boundary_mesh->partition(_n_patches);
  }

  // prepare sideset names and boundary_ids added to mesh
  std::vector<BoundaryName> sideset_names =
      sidesetNameHelper(boundary_info.get_sideset_name(_sideset));

  std::vector<boundary_id_type> boundary_ids =
      MooseMeshUtils::getBoundaryIDs(*mesh, sideset_names, true);

  mooseAssert(sideset_names.size() == _n_patches,
              "sideset_names must have as many entries as user-requested number of patches.");
  mooseAssert(boundary_ids.size() == _n_patches,
              "boundary_ids must have as many entries as user-requested number of patches.");

  // loop through all elements in the boundary mesh and assign the side of
  // the _original_ element to the new sideset
  for (const auto & elem : boundary_mesh->active_element_ptr_range())
  {
    if (boundary_elem_to_mesh_elem.find(elem->id()) == boundary_elem_to_mesh_elem.end())
      mooseError("Element in the boundary mesh with id ",
                 elem->id(),
                 " not found in boundary_elem_to_mesh_elem.");

    auto side = boundary_elem_to_mesh_elem.find(elem->id())->second;

    mooseAssert(elem->processor_id() < boundary_ids.size(),
                "Processor id larger than number of patches.");
    boundary_info.add_side(
        std::get<0>(side), std::get<1>(side), boundary_ids[elem->processor_id()]);
  }

  // make sure new boundary names are set
  for (unsigned int j = 0; j < boundary_ids.size(); ++j)
  {
    boundary_info.sideset_name(boundary_ids[j]) = sideset_names[j];
    boundary_info.nodeset_name(boundary_ids[j]) = sideset_names[j];
  }

  return mesh;
}

void
PatchSidesetGenerator::partition(MeshBase & mesh) const
{
  if (_partitioner_name != "grid")
    return;

  // subdivide by length along the curve
  if (_dim == 1)
    mooseError("1D damnit");
  else
  {
    // compute the target area for each partition
    Real target_area = 0;
    for (auto & elem_ptr : mesh.active_element_ptr_range())
      target_area += elem_ptr->volume();
    target_area /= _n_patches;

    // initial partition
    auto begin = mesh.active_local_elements_begin();
    auto end = mesh.active_local_elements_end();

    // find element with lowest degree
    unsigned int min_degree = 1000;
    const Elem * min_degree_elem = nullptr;
    for (auto it = begin; it != end; ++it)
    {
      unsigned int degree = 0;
      const Elem * elem = *it;
      for (unsigned int j = 0; j < elem->n_sides(); ++j)
        if (elem->neighbor_ptr(j))
          ++degree;

      if (degree < min_degree)
      {
        min_degree = degree;
        min_degree_elem = elem;
      }
    }

    // initial point distribution, add min_degree elem as first center
    std::vector<const Elem *> bubble_centers = {mostDistantElement({min_degree_elem})};
    for (unsigned int j = 1; j < _n_patches; ++j)
    {
      const Elem * new_center = mostDistantElement(bubble_centers);
      bubble_centers.push_back(new_center);
    }

    // main bubble growth and recentering loop
    std::vector<Bubble> bubbles;
    bool keep_going = true;
    while (keep_going)
    {
      // grow bubbles from existing centers
      growBubbles(bubble_centers, bubbles);

      // find new centroid elements for all bubbles
      std::vector<const Elem *> new_bubble_centers(_n_patches);
      for (unsigned int j = 0; j < _n_patches; ++j)
        new_bubble_centers[j] = bubbles[j].bubbleCentroidElem();

      // has anything changed?; if so keep going
      keep_going = false;
      for (unsigned int j = 0; j < _n_patches; ++j)
        if (bubble_centers[j] != new_bubble_centers[j])
          keep_going = true;

      // copy over bubble centers
      bubble_centers = new_bubble_centers;
    }

//for (unsigned int j = 0; j < _n_patches; ++j)
//{
//  std::cout << "Bubble " << j << ": ";
//  for (auto & e : bubbles[j].members())
//    std::cout << e->id() << " ";
//  std::cout << std::endl;
//}

    // compare centroid elems new vs old
    //for (unsigned int j = 0; j < _n_patches; ++j)
    //  std::cout << bubble_centers[j]->id() << " " << new_bubble_centers[j]->id() << std::endl;

    // assign partition
    for (unsigned int j = 0; j < _n_patches; ++j)
      for (auto & e : bubbles[j].members())
        mesh.elem_ptr(e->id())->processor_id() = j;
  }
}

void
PatchSidesetGenerator::growBubbles(const std::vector<const Elem *> & seeds,
                                   std::vector<Bubble> & bubbles) const
{
  // prepare bubbles & insert seeds
  bubbles.clear();
  bubbles.resize(seeds.size());
  std::set<const Elem *> assigned_elems;
  for (unsigned int j = 0; j < seeds.size(); ++j)
  {
    bubbles[j] = PatchSidesetGenerator::Bubble(seeds[j]);
    assigned_elems.insert(seeds[j]);
  }

  // initialize frontier
  for (auto & b : bubbles)
    b.computeFrontier(assigned_elems);

  // main loop for growing bubbles
  while (assigned_elems.size() < _n_boundary_mesh_elems)
  {
    // find bubble with smallest workload
    unsigned int min_b = 0;
    Real min_workload = std::numeric_limits<double>::max();
    for (unsigned int j = 0; j < bubbles.size(); ++j)
      if (bubbles[j].workload() < min_workload && bubbles[j].frontierSize() > 0)
      {
        min_workload = bubbles[j].workload();
        min_b = j;
      }

    // find the element with shortest distance from seed
    const Elem * candidate = bubbles[min_b].closestFrontierElem();

    // add to
    assigned_elems.insert(candidate);

    // modify all bubbles
    for (unsigned int j = 0; j < bubbles.size(); ++j)
      if (j == min_b)
        bubbles[j].addToMember(candidate, assigned_elems);
      else
        bubbles[j].removeFromFrontier(candidate);
  }
}

const Elem *
PatchSidesetGenerator::mostDistantElement(const std::vector<const Elem *> & seeds) const
{
  Real max_distance = 0;
  const Elem * most_distant_elem = nullptr;

  // fill visited and queue with the seeds &
  std::set<const Elem *> visited;
  std::queue<const Elem *> next_visits;
  for (auto & e : seeds)
  {
    visited.insert(e);
    next_visits.push(e);
  }

  // main loop of BFS
  while (!next_visits.empty())
  {
    // select element from queue and remove it
    const Elem * elem = next_visits.front();
    next_visits.pop();

    // Is this element the most distant from the seeds?
    Real distance = std::numeric_limits<double>::max();
    for (auto & e : seeds)
    {
      Real d = (elem->centroid() - e->centroid()).norm();
      if (d < distance)
        distance = d;
    }
    if (distance > max_distance)
    {
      max_distance = distance;
      most_distant_elem = elem;
    }

    // look at the neighbors of elem
    for (unsigned int j = 0; j < elem->n_sides(); ++j)
    {
      const Elem * neigh = elem->neighbor_ptr(j);
      if (neigh && visited.find(neigh) == visited.end())
      {
        // remember that neigh was added and will be visited
        visited.insert(neigh);
        // inert into the queue so we actually visit neigh
        next_visits.push(neigh);
      }
    }
  }
  return most_distant_elem;
}

std::vector<BoundaryName>
PatchSidesetGenerator::sidesetNameHelper(const std::string & base_name) const
{
  std::vector<BoundaryName> rv;
  for (unsigned int j = 0; j < _n_patches; ++j)
  {
    std::stringstream ss;
    ss << base_name << "_" << j;
    rv.push_back(ss.str());
  }
  return rv;
}

Elem *
PatchSidesetGenerator::boundaryElementHelper(MeshBase & mesh, libMesh::ElemType type) const
{
  switch (type)
  {
    case 0:
      return mesh.add_elem(new libMesh::Edge2);
    case 1:
      return mesh.add_elem(new libMesh::Edge3);
    case 2:
      return mesh.add_elem(new libMesh::Edge4);
    case 3:
      return mesh.add_elem(new libMesh::Tri3);
    case 4:
      return mesh.add_elem(new libMesh::Tri6);
    case 5:
      return mesh.add_elem(new libMesh::Quad4);
    case 6:
      return mesh.add_elem(new libMesh::Quad8);
    case 7:
      return mesh.add_elem(new libMesh::Quad9);
    default:
      mooseError("Unsupported element type (libMesh elem_type enum): ", type);
  }
}

PatchSidesetGenerator::Bubble::Bubble() : _workload(0), _seed(nullptr) { _member_elems.clear(); }

PatchSidesetGenerator::Bubble::Bubble(const Elem * seed) : _workload(seed->volume()), _seed(seed)
{
  _member_elems.clear();
  _member_elems.insert(seed);
}

void
PatchSidesetGenerator::Bubble::computeFrontier(const std::set<const Elem *> & assigned_elems)
{
  _frontier.clear();
  for (auto & elem : _member_elems)
    for (unsigned int j = 0; j < elem->n_sides(); ++j)
    {
      const Elem * neigh = elem->neighbor_ptr(j);
      if (neigh && assigned_elems.find(neigh) == assigned_elems.end() &&
          _member_elems.find(neigh) == _member_elems.end())
        _frontier.insert(neigh);
    }
}

void
PatchSidesetGenerator::Bubble::removeFromFrontier(const Elem * elem)
{
  const auto & it = std::find(_frontier.begin(), _frontier.end(), elem);
  if (it != _frontier.end())
    _frontier.erase(it);
}

void
PatchSidesetGenerator::Bubble::addToMember(const Elem * elem,
                                           const std::set<const Elem *> & assigned_elems)
{
  // remove from frontier
  const auto & it = std::find(_frontier.begin(), _frontier.end(), elem);
  if (it == _frontier.end())
    ::mooseError("Elem ", elem->id(), " is not in the frontier of this bubble.");

  _frontier.erase(it);
  // add to member elems
  if (std::find(_member_elems.begin(), _member_elems.end(), elem) != _member_elems.end())
    ::mooseError("Elem ", elem->id(), " is about to be added to bubble but already exists.");
  _member_elems.insert(elem);
  _workload += elem->volume();

  // update the frontier
  for (unsigned int j = 0; j < elem->n_sides(); ++j)
  {
    const Elem * neigh = elem->neighbor_ptr(j);
    if (neigh && assigned_elems.find(neigh) == assigned_elems.end() &&
        _member_elems.find(neigh) == _member_elems.end())
      _frontier.insert(neigh);
  }
}

const Elem *
PatchSidesetGenerator::Bubble::closestFrontierElem() const
{
  Real distance = std::numeric_limits<double>::max();
  const Elem * closest_frontier_elem = nullptr;
  for (auto & elem : _frontier)
  {
    Real d = (_seed->centroid() - elem->centroid()).norm();
    if (d < distance)
    {
      distance = d;
      closest_frontier_elem = elem;
    }
  }
  return closest_frontier_elem;
}

const Elem *
PatchSidesetGenerator::Bubble::bubbleCentroidElem() const
{
  // compute centroid
  Point centroid(0, 0, 0);
  Real volume = 0;
  for (auto & e : _member_elems)
  {
    centroid += e->volume() * e->centroid();
    volume += e->volume();
  }
  centroid /= volume;

  // find closest elem in _member_elems
  const Elem * centroid_elem = _seed;
  Real distance = std::numeric_limits<double>::max();
  for (auto & e : _member_elems)
  {
    Real d = (centroid - e->centroid()).norm();
    if (d < distance)
    {
      distance = d;
      centroid_elem = e;
    }
  }

  return centroid_elem;
}
