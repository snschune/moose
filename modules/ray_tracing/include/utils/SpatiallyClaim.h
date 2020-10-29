// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html
//
// #pragma once
//
// // MOOSE includes
// #include "MooseMesh.h"
//
// // libMesh includes
// #include "libmesh/elem.h"
// #include "libmesh/mesh_tools.h"
// #include "libmesh/parallel_algebra.h"
// #include "libmesh/parallel_sync.h"
// #include "libmesh/remote_elem.h"
//
// // System includes
// #include <functional>
// #include <unordered_map>
//
// /**
//  * Spatially claims a vector of objects.
//  */
// template <typename ObjectType, typename ReceiveContext>
// class SpatiallyClaim
// {
// public:
//   /**
//    * Constructor.
//    * @param mesh The MooseMesh
//    * @param all_objects The vector of objects that need to be claimed
//    * @param local_objects Insertion point for objects that have been claimed
//    * @param do_exchange Whether or not an exhange is needed, i.e., if all_objects still needs to
//    be
//    * filled by objects on other processors
//    * @param context A context to be used in receiving packed data (if necessary)
//    */
//   SpatiallyClaim(MooseMesh & mesh,
//                  const std::vector<ObjectType> & all_objects,
//                  std::vector<ObjectType> & local_objects,
//                  const bool do_exchange = false,
//                  ReceiveContext * context = nullptr)
//     : _mesh(mesh),
//       _comm(_mesh.getMesh().comm()),
//       _pid(_comm.rank()),
//       _all_objects(all_objects),
//       _local_objects(local_objects),
//       _do_exchange(do_exchange),
//       _context(context)
//   {
//   }
//
//   /**
//    * Claims objects in _all_objects into _local_objects
//    */
//   void claim()
//   {
//     preClaim();
//
//     _local_objects.clear();
//
//     // Grab the point locator
//     _point_locator = _mesh.getMesh().sub_point_locator();
//     _point_locator->enable_out_of_mesh_mode();
//
//     // Exchange: filter objects into processors that _may_ claim them
//     std::unordered_map<processor_id_type, std::vector<ObjectType>> objs_to_send;
//     if (_do_exchange)
//       for (auto & obj : _all_objects)
//       {
//         const auto point = getPoint(obj);
//         for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
//           if (_inflated_bboxes[pid].contains_point(point))
//             objs_to_send[pid].push_back(obj);
//       }
//
//     // Functor for possibly claiming a vector of objects
//     std::function<void(processor_id_type, const std::vector<ObjectType> &)> claim_functor =
//         [&](processor_id_type /* pid */, const std::vector<ObjectType> & objs) {
//           for (auto & obj : objs)
//             possiblyClaim(obj);
//         };
//
//     if (_do_exchange) // Send the relevant objects to everyone and then claim
//       pushData(objs_to_send, claim_functor);
//     else // Already have the relevant objects, just claim
//       claim_functor(_pid, _all_objects);
//
//     postClaim();
//   }
//
//   /**
//    * Initializes the bounding boxes.
//    *
//    * Same as reinit().
//    */
//   void init()
//   {
//     buildBoundingBoxes();
//     buildPointNeighbors();
//   }
//
// protected:
//   /**
//    * Entry point for acting on an object before it is attempted to be claimed.
//    */
//   void prePossiblyClaim(const ObjectType & /* obj */){};
//   /**
//    * Entry point for acting on an object after it is claimed.
//    */
//   void postClaim(ObjectType & /* obj */, const Elem * /* elem */){};
//   /**
//    * Entry point at the beginning of claiming all objects.
//    */
//   void preClaim(){};
//   /**
//    * Entry point at the end of claiming all objects.
//    */
//   void postClaim(){};
//
//   /**
//    * Builds the bounding boxes (_bbox, _inflated_bboxes, _inflated_neighbor_bboxes).
//    */
//   void buildBoundingBoxes()
//   {
//     // Local bounding box
//     _bbox = MeshTools::create_local_bounding_box(_mesh.getMesh());
//     _global_bbox = _bbox;
//
//     // Gather the bounding boxes of all processors
//     std::vector<std::pair<Point, Point>> bb_points = {static_cast<std::pair<Point,
//     Point>>(_bbox)}; _comm.allgather(bb_points, true); _inflated_bboxes.resize(_comm.size()); for
//     (processor_id_type pid = 0; pid < _comm.size(); ++pid)
//     {
//       const BoundingBox pid_bbox = static_cast<BoundingBox>(bb_points[pid]);
//       _inflated_bboxes[pid] = inflateBoundingBox(pid_bbox);
//       _global_bbox.union_with(pid_bbox);
//     }
//
//     // Find intersecting (neighbor) bounding boxes
//     _inflated_neighbor_bboxes.clear();
//     for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
//     {
//       // Skip this processor
//       if (pid == _pid)
//         continue;
//       // Insert if the searched processor's bbox intersects my bbox
//       const auto & pid_bbox = _inflated_bboxes[pid];
//       if (_bbox.intersects(pid_bbox))
//         _inflated_neighbor_bboxes.emplace_back(pid, pid_bbox);
//     }
//   }
//
//   /**
//    * Build the map of elements to all of their point neighbors
//    *
//    * TODO: Move this eventually into MooseMesh, MeshBase, or FEProblemBase
//    */
//   void buildPointNeighbors()
//   {
//     _elem_point_neighbors.clear();
//     const auto & node_to_elem_map = _mesh.nodeToElemMap();
//
//     for (const auto & elem : _mesh.getMesh().active_element_ptr_range())
//     {
//       auto & fill = _elem_point_neighbors[elem->id()];
//
//       const auto first_order_type = Elem::first_order_equivalent_type(elem->type());
//       if (first_order_type == QUAD4)
//         fill.reserve(8);
//       else if (first_order_type == HEX8)
//         fill.reserve(26);
//       else if (first_order_type == TRI3)
//         fill.reserve(10);
//       else if (first_order_type == TET4)
//         fill.reserve(102); // Beware: this is really really bad!
//
//       for (unsigned int v = 0; v < elem->n_vertices(); ++v)
//       {
//         const auto & node = elem->node_ptr(v);
//         for (const auto & neighbor_id : node_to_elem_map.at(node->id()))
//         {
//           if (neighbor_id == elem->id())
//             continue;
//
//           const auto & neighbor = _mesh.elemPtr(neighbor_id);
//           if (std::count(fill.begin(), fill.end(), neighbor) == 0)
//             fill.emplace_back(neighbor);
//         }
//       }
//     }
//   }
//
//   /**
//    * Possibly claim an object.
//    */
//   void possiblyClaim(const ObjectType & obj)
//   {
//     const auto & point = getPoint(obj);
//     const auto id = getID(obj);
//
//     prePossiblyClaim(obj);
//
//     const auto elem = claimPoint(point, id, (*_point_locator)(point));
//     if (elem)
//     {
//       _local_objects.push_back(obj);
//       postClaim(_local_objects.back(), elem);
//     }
//   }
//
//   /**
//    * Try to claim a spatial point.
//    *
//    * @param point The point to claim
//    * @param id An ID associated with the point
//    * @param elem The local element to first consider for this processor's ownership
//    * @return The element that contains the point if we claim the point, nullptr if we don't claim
//    it
//    */
//   const Elem * claimPoint(const Point & point, const dof_id_type id, const Elem * elem)
//   {
//     if (elem)
//     {
//       // Looking for smallest (even ID object) or largest (odd ID object) elem id
//       const bool smallest = id % 2 == 0;
//
//       // Start with the element we found, as it is a valid candidate
//       const Elem * extremum_elem = elem;
//
//       // All point neighbors for this element
//       mooseAssert(_elem_point_neighbors.count(elem->id()), "Not in point neighbor map");
//       const auto & neighbors = _elem_point_neighbors.at(elem->id());
//
//       // Find element that matches the extremum criteria
//       for (const auto & neighbor : neighbors)
//       {
//         mooseAssert(neighbor->active(), "Inactive neighbor");
//
//         if ((smallest && neighbor->id() < extremum_elem->id()) || // satisfies
//             (!smallest && neighbor->id() > extremum_elem->id()))  // ...one of the id checks
//           if (neighbor->contains_point(point))                    // and also contains the point
//             extremum_elem = neighbor;
//       }
//
//       // Claim the object if we own the extremum elem
//       if (extremum_elem->processor_id() == _pid)
//       {
//         mooseAssert(extremum_elem->active(), "Inactive element");
//         return extremum_elem;
//       }
//     }
//
//     return nullptr;
//   }
//
//   /**
//    * Pushes the data in data (whose key is the processor that each vector of ObjectType that
//    needs
//    * to be communicated) and allows for acting on the data using functor after the data has been
//    * received.
//    *
//    * This _MUST_ be specialized if you are using SpatiallyClaim with the exhange feature enabled
//    * (_do_exhange = true), as this method depends on ObjectType.
//    *
//    * In cases where ObjectType is castable to a StandardType, using a TIMPI
//    * push_parallel_vector_data call will do:
//    * TIMPI::push_parallel_vector_data(_comm, data, functor);
//    *
//    * In cases where ObjectType is more complex, packing routines must be created and
//    * TIMPI::push_parallel_packed_range can be used instead.
//    */
//   void
//   pushData(const std::unordered_map<processor_id_type, std::vector<ObjectType>> & /* data */,
//            std::function<void(processor_id_type, const std::vector<ObjectType> &)> & /* functor
//            */)
//   {
//     mooseError("Unimplemented SpatiallyClaim::pushData");
//   }
//
//   /**
//    * Gets the point associated with an object. This needs to be specialized for object type.
//    */
//   static Point getPoint(const ObjectType & /* obj */)
//   {
//     mooseError("Unimplemented SpatiallyClaim::getPoint");
//   }
//
//   /**
//    * Gets an ID associated with an object. This is not required, but should be specialized for
//    * object type if possible to aid in parallel efficiency.
//    */
//   static dof_id_type getID(const ObjectType & /* obj */) { return 0; }
//
//   /**
//    * Creates an inflated bounding box.
//    */
//   static BoundingBox inflateBoundingBox(BoundingBox bbox, Real multiplier = 0.01)
//   {
//     Real amount = multiplier * (bbox.max() - bbox.min()).norm();
//     Point inflation(amount, amount, amount);
//     bbox.first -= inflation;
//     bbox.second += inflation;
//     return bbox;
//   }
//
//   /// The mesh
//   MooseMesh & _mesh;
//   /// The communicator
//   const libMesh::Parallel::Communicator & _comm;
//   /// This processor ID
//   const processor_id_type _pid;
//
//   /// The objects that need to be searched to possibly claimed
//   const std::vector<ObjectType> & _all_objects;
//   /// The local objects that are claimed
//   std::vector<ObjectType> & _local_objects;
//   /// Whether or not the data needs to be initially exchanged
//   const bool _do_exchange;
//   /// Context used for receiving a packed range (not used for StandardType castable objects)
//   ReceiveContext * _context;
//
//   /// The point locator
//   std::unique_ptr<PointLocatorBase> _point_locator = nullptr;
//
//   /// The bounding box for this processor
//   BoundingBox _bbox;
//   /// The global bounding box
//   BoundingBox _global_bbox;
//   /// The inflaed bounding boxes for all processors
//   std::vector<BoundingBox> _inflated_bboxes;
//   /// Inflated bounding boxes that are neighboring to this processor (pid : bbox for each entry)
//   std::vector<std::pair<processor_id_type, BoundingBox>> _inflated_neighbor_bboxes;
//
//   /// Map of point neighbors for each element
//   std::unordered_map<dof_id_type, std::vector<const Elem *>> _elem_point_neighbors;
// };
