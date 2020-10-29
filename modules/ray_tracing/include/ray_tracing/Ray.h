//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Local includes
#include "RayTracingCommon.h"

// MOOSE includes
#include "MooseError.h"

// libMesh Includes
#include "libmesh/parallel.h"

// Forward declarations
namespace libMesh
{
class Elem;
}
class RayTracingStudy;
// Friend access to Ray
class TraceRay;
// Friend access to ChangeDirectionKey for accessing changeDirection()
class RayBoundaryConditionBase;
// Friend access to ChangeStartDirectionKey for accessing changeStartDirection()
class RayKernelBase;
// Friend access to NonResetCountersKey for accessing constructor/reset without counter reset
namespace MooseUtils
{
template <class T>
class SharedPool;
}

/// Type for a Ray's ID
typedef unsigned long int RayID;
/// Type for a Ray's data
#ifdef SINGLE_PRECISION_RAY
typedef float RayData;
#else
typedef libMesh::Real RayData;
#endif
/// Type for the index into the data and aux data on a Ray
typedef unsigned int RayDataIndex;

/**
 * Basic datastructure for a ray that will traverse the mesh.
 */
class Ray
{
public:
  /**
   * Class that is used as a parameter to changeDirection() that allows only
   * RayBC methods to call changeDirection()
   */
  class ChangeDirectionKey
  {
    friend class RayBoundaryConditionBase;
    ChangeDirectionKey() {}
    ChangeDirectionKey(const ChangeDirectionKey &) {}
  };

  /**
   * Class that is used as a parameter to changeStartDirection() that allows only
   * RayKernelBase methods to call changeStartDirection()
   */
  class ChangeStartDirectionKey
  {
    friend class RayKernelBase;
    ChangeStartDirectionKey() {}
    ChangeStartDirectionKey(const ChangeStartDirectionKey &) {}
  };

  /**
   * Class that is used as a parameter to the constructor/reset method that does
   * not reset a Ray's counters.
   */
  class NonResetCountersKey
  {
    friend class Parallel::Packing<std::shared_ptr<Ray>>;
    friend class MooseUtils::SharedPool<Ray>;
    NonResetCountersKey() {}
    NonResetCountersKey(const NonResetCountersKey &) {}
  };

  /**
   * Class that is used as a parameter to the public resetCounters() method that
   * only allows a RayTracingStudy to call resetCounters().
   */
  class ResetCountersKey
  {
    friend class RayTracingStudy;
    ResetCountersKey() {}
    ResetCountersKey(const ResetCountersKey &) {}
  };

  /**
   * Basic constructor for a Ray
   */
  Ray();

  /**
   * Completely reset a Ray.
   *
   * Also used by SharedPool to completely reset a Ray.
   */
  void reset();

  /**
   * Ray constructor with a size of data to be zeroed
   * @param data_size Size of data to initialize with zeros
   * @param aux_data_size Size of aux data to initialize with zeros
   * @param starting_elem Element the Ray starts in
   */
  Ray(const std::size_t data_size, const std::size_t aux_data_size);

  /**
   * Ray reset with a size of data to be zeroed
   *
   * Used by the SharedPool to reset a Ray with these forwarded arguments
   *
   * @param data_size Size of data to initialize with zeros
   * @param aux_data_size Size of aux data to initialize with zeros
   */
  void reset(const std::size_t data_size, const std::size_t aux_data_size);

  /**
   * Ray constructor that only sizes the data and does not reset the counters
   *
   * THIS IS ONLY USED INTERNALLY. You cannot and should not use it:
   * While this method is public, it is protected by the NonResetCountersKey,
   * which makes it only possible to be called by the Parallel unpacking routines.
   * This is to force users to only acquire/construct new Rays with a clean
   * history (intersections, direction, trajectory information, etc) in order
   * to not produce undefined behavior in tracing.
   */
  Ray(const std::size_t data_size, const std::size_t aux_data_size, const NonResetCountersKey &);

  /**
   * Ray reset that only resizes the data and does not reset the resetCounters
   *
   * Used by the SharedPool to reset a Ray with these forwarded arguments
   *
   * THIS IS ONLY USED INTERNALLY. You cannot and should not use it:
   * While this method is public, it is protected by the NonResetCountersKey,
   * which makes it only possible to be called by the Parallel unpacking routines.
   * This is to force users to only acquire/construct new Rays with a clean
   * history (intersections, direction, trajectory information, etc) in order
   * to not produce undefined behavior in tracing.
   */
  void
  reset(const std::size_t data_size, const std::size_t aux_data_size, const NonResetCountersKey &);

  /**
   * Copy constructor that copies another Ray.
   *
   * Does reset the counters.
   */
  Ray(const Ray & ray);

  /**
   * Reset a Ray by copying from another.
   *
   * Used by the SharedPool to reset a Ray with the trajectory
   * information from another Ray.
   *
   * Does reset the counters.
   */
  void reset(const Ray & ray);

  /**
   * Copy assignment operator.
   *
   * This WILL reset the counters on the Ray.
   */
  Ray & operator=(const Ray & other);
  /**
   * Equality operator.
   *
   * This will perform "fuzzy" equals checks on the points and data.
   */
  bool operator==(const Ray & other) const;
  /**
   * Non-equal operator.
   *
   * This will perform "fuzzy" equals checks on the points and data.
   */
  bool operator!=(const Ray & other) const;

  /// Invalid index into a Ray's data
  static const RayDataIndex INVALID_RAY_DATA_INDEX = static_cast<RayDataIndex>(-1);
  /// Invalid Ray ID
  static const RayID INVALID_RAY_ID = static_cast<RayID>(-1);

  /**
   * Sets the Ray's ID
   */
  void setID(const RayID id);
  /**
   * Gets the Ray's ID
   */
  RayID id() const { return _id; }
  /**
   * Whether or not the Ray's ID is invalid
   */
  bool invalidID() const { return _id == INVALID_RAY_ID; }
  /**
   * Invalidate the Ray's ID
   */
  void invalidateID();

  /**
   * Gets the Ray's start point
   *
   * Can only be called before a Ray has started to be traced!
   */
  const Point & startPoint() const;
  /**
   * Gets whether or not the Ray's start point is invalid
   *
   * Can only be called before a Ray has started to be traced!
   */
  bool invalidStartPoint() const;

  /**
   * Gets the point that the Ray is currently at.
   *
   * Before a Ray is traced, this is the starting point of the Ray.
   * While a Ray is being traced, this is the furthest point that the Ray
   * has travelled. During RayKernel execution, this is the end of the segment.
   * After a Ray has been traced, this is the point where the Ray was killed.
   */
  const Point & currentPoint() const { return _current_point; }
  /**
   * Whether or not the point that the Ray is currently at is valid.
   */
  bool invalidCurrentPoint() const { return _current_point == RayTracingCommon::invalid_point; }

  /**
   * Sets the Ray to start from \p start and to be killed when it hits \p end.
   *
   * Internally, this sets the max distance of the Ray to || start - end || and
   * sets endSet() true to keep track of the fact that user set an end point.
   *
   * Can only be called before a Ray has started to be traced!
   */
  void setStartEnd(const Point & start, const Point & end);

  /**
   * Sets the Ray's start (as the current point) and direction of travel.
   *
   * Can only be called before a Ray has started to be traced!
   */
  void setStartDirection(const Point & start, const Point & direction);

  /**
   * This method is for internal use only. It is intended to be called only by
   * RayBoundaryConditionBase::changeRayDirection().
   *
   * If you wish to change a Ray's direction mid-trace in a RayBC, see
   * RayBoundaryConditionBase::changeRayDirection() instead.
   *
   * ChangeDirectionKey is constructable only by RayBC objects on purpose to limit usage of this
   * method.
   */
  void changeDirection(const Point & direction, const ChangeDirectionKey);

  /**
   * This method is for internal use only. It is intended to be called only by
   * RayKernelBase::changeRay().
   *
   * If you wish to change a Ray's direction mid-trace in a RayKernel, see
   * RayKernelBase::changeRay() instead.
   *
   * ChangeStartDirectionKey is constructable only by RayKernelBase objects on purpose to limit
   * usage of this method.
   */
  void
  changeStartDirection(const Point & start, const Point & direction, const ChangeStartDirectionKey);

  /**
   * Gets the Ray's direction
   */
  const Point & direction() const { return _direction; }
  /**
   * Whether or not the Ray's direction is set to invalid.
   */
  bool invalidDirection() const { return _direction == RayTracingCommon::invalid_point; }

  /**
   * Gets a writeable reference to the Ray's data
   */
  std::vector<RayData> & data() { return _data; }
  /**
   * Gets a read only reference to the Ray's data
   */
  const std::vector<RayData> & data() const { return _data; }
  /**
   * Gets a writeable reference to the Ray's data at an index
   */
  RayData & data(const std::size_t i)
  {
    mooseAssert(_data.size() > i, "Accessing Ray data out of range");
    return _data[i];
  }
  /**
   * Gets a read only reference to the Ray's data at an index
   */
  const RayData & data(const std::size_t i) const
  {
    mooseAssert(_data.size() > i, "Accessing Ray data out of range");
    return _data[i];
  }

  /**
   * Gets a writeable reference to the Ray's auxilary data
   */
  std::vector<RayData> & auxData() { return _aux_data; }
  /**
   * Gets a read only reference to the Ray's auxilary data
   */
  const std::vector<RayData> & auxData() const { return _aux_data; }
  /**
   * Gets a writeable reference to a component of the Ray's auxilary data
   */
  RayData & auxData(const std::size_t i)
  {
    mooseAssert(_aux_data.size() > i, "Accessing Ray aux data out of range");
    return _aux_data[i];
  }
  /**
   * Gets a read only reference to a component of the Ray's auxilary data
   */
  const RayData & auxData(const std::size_t i) const
  {
    mooseAssert(_aux_data.size() > i, "Accessing Ray aux data out of range");
    return _aux_data[i];
  }
  /**
   * Sets a component of the Ray's auxilary data
   */
  void setAuxData(const std::size_t i, const RayData value)
  {
    mooseAssert(_aux_data.size() > i, "Accessing Ray aux data out of range");
    _aux_data[i] = value;
  }
  /**
   * Sets multiple components of the Ray's auxilary data starting at the start_i position
   */
  void setAuxData(const std::size_t start_i, const std::vector<RayData> & values);

  /**
   * Sets the element that the Ray should start in
   *
   * Can only be called before a Ray has started to be traced!
   */
  void setStartingElem(const Elem * current_elem);
  /**
   * Gets the element that the Ray should start in
   *
   * Can only be called before a Ray has started to be traced!
   */
  const Elem * startingElem() const;
  /**
   * Invalidate the Ray's starting elem. Also invalides the starting incoming side.
   *
   * Can only be called before a Ray has started to be traced!
   */
  void invalidateStartingElem();

  /**
   * Gets the current element that the Ray is in.
   *
   * Before tracing, this is the starting element for the Ray.
   *
   * During tracing:
   * - It is valid within RayKernels, because a Ray can only operate on a
   *   single element per segment.
   * - When used on boundaries (RayBCs), it is the element that the
   *   Ray actually traced through. When on a boundary, a RayBC may be
   *   applied to multiple elements when hitting a vertex or edge.
   *   Therefore, this will be only one of said elements.
   */
  const Elem * currentElem() const { return _current_elem; }

  /**
   * Sets the Ray's starting incoming side on startingElem()
   *
   * Can only be called before a Ray has started to be traced!
   */
  void setStartingIncomingSide(const unsigned short side);
  /**
   * Gets the Ray's starting incoming side on startingElem()
   *
   * Can only be called before a Ray has started to be traced!
   */
  unsigned short startingIncomingSide() const;
  /**
   * Whether or not the Ray's current incoming side is invalid
   *
   * Can only be called before a Ray has started to be traced!
   */
  bool invalidStartingIncomingSide() const;
  /**
   * Invalidate the Ray's current incoming side
   *
   * Can only be called before a Ray has started to be traced!
   */
  void invalidateStartingIncomingSide();
  /**
   * Get a Ray's current incoming side
   *
   * Before tracing, this is the starting incoming side (if any).
   *
   * During and after tracing, this is ONLY guaranteed to be valid
   * while executing RayKernels!
   */
  unsigned short currentIncomingSide() const { return _current_incoming_side; }
  /**
   * Whether or not the Ray's current incoming side is invalid
   *
   * Before tracing, this is the starting incoming side (if any).
   *
   * During and after tracing, this is ONLY guaranteed to be valid
   * while executing RayKernels!
   */
  bool invalidCurrentIncomingSide() const
  {
    return _current_incoming_side == RayTracingCommon::invalid_side;
  }

  /**
   * Sets the maximum distance this Ray should travel.
   *
   * If setting a Ray's trajectory with setStartEnd(), the max distance
   * is set internally to be || end - start ||.
   *
   * Can only be called before a Ray has started to be traced!
   */
  void setMaxDistance(const Real max_distance);

  /**
   * Whether or not the user has set and end point for this Ray.
   *
   * This is done by limiting the distance of the Ray in its set direction.
   */
  bool endSet() const { return _end_set; }
  /**
   * Whether or not the Ray is at the user-defined end point.
   *
   * This is only valid when the user set the trajectory of the Ray
   * with setStartEnd(), which internally set its maximum distance
   * to the straight-line distance between start and end and set
   * endSet() == true.
   */
  bool atEnd() const;

  /**
   * Gets the number of times this Ray has crossed a processor
   */
  unsigned int processorCrossings() const { return _processor_crossings; }

  /**
   * Gets the number of intersections this Ray has done
   */
  unsigned int intersections() const { return _intersections; }

  /**
   * Gets the distance this Ray has traveled
   */
  Real distance() const { return _distance; }
  /**
   * Gets the max distance this Ray is allowed to travel
   *
   * This may be set internally to || end - start || in the case that
   * the user initialized the Ray with setStartEnd().
   */
  Real maxDistance() const { return _max_distance; }

  /**
   * Whether or not this Ray should continue
   */
  bool shouldContinue() const { return _should_continue; }
  /**
   * Sets whether or not this Ray should continue
   */
  void setShouldContinue(const bool should_continue) { _should_continue = should_continue; }

  /**
   * Whether or not this Ray has had its trajectory changed
   */
  bool trajectoryChanged() const { return _trajectory_changed; }
  /**
   * Gets the number of trajectory changes this Ray has had
   */
  unsigned int trajectoryChanges() const { return _trajectory_changes; }

  /**
   * Helper function for getting information about the Ray
   *
   * The RayTracingStudy may be optionally provided to add additional information
   * (the registered names associated with each data and aux data entry)
   */
  std::string getInfo(const RayTracingStudy * study = nullptr) const;

  /**
   * Reset all of the internal counters
   *
   * THIS IS ONLY USED INTERNALLY. You cannot and should not use it:
   * It is protected by the ResetCountersKey which is only constructable by the
   * RayTracingStudy. This is so that only the RayTracingStudy can call this
   * method at times when it is appropriate to reset the counters in order
   * to avoid undefined behavior in the ray-tracing routines.
   */
  void resetCounters(const ResetCountersKey);

  /**
   * Whether or not a Ray has begun tracing
   */
  bool hasTraced() const { return distance() || processorCrossings(); }

private:
  /**
   * Changes the Ray's ID
   */
  void changeID(const RayID id) { _id = id; }

  /**
   * Sets the Ray's current point
   */
  void setCurrentPoint(const Point & current_point) { _current_point = current_point; }

  /**
   * Sets the Ray's direction
   */
  void setDirection(const Point & direction) { _direction = direction; }

  /**
   * Change a Ray's current elem
   */
  void setCurrentElem(const Elem * current_elem) { _current_elem = current_elem; }

  /**
   * Change a Ray's incoming side
   */
  void setCurrentIncomingSide(const unsigned short current_incoming_side)
  {
    _current_incoming_side = current_incoming_side;
  }

  /**
   * Set whether or not this Ray has had its trajectory changed
   */
  void setTrajectoryChanged(const bool trajectory_changed)
  {
    _trajectory_changed = trajectory_changed;
  }

  /**
   * Increment the Ray's processor crossing counter
   */
  void addProcessorCrossing() { ++_processor_crossings; }

  /**
   * Increment the Ray's intersection counter
   */
  void addIntersection() { ++_intersections; }

  /**
   * Increment the Ray's trajectory change counter
   */
  void addTrajectoryChange() { ++_trajectory_changes; }

  /**
   * Adds to the distance this Ray has traveled
   */
  void addDistance(const Real add_distance) { _distance += add_distance; }
  /**
   * Changes the Ray's max distance to be traveled
   */
  void changeMaxDistance(const Real max_distance) { _max_distance = max_distance; }

  /**
   * Produces a useful error if a Ray has started tracing
   */
  void errorIfTracing(const std::string & reason) const;

  /**
   * Reset all of the internal counters
   */
  void resetCounters();

  /**
   * Helper for the equality operators.
   */
  bool equalityHelper(const Ray & other, const bool equal) const;

  /// A unique ID for this Ray
  RayID _id;

  /**
   * Current point of the Ray.
   *
   * Before tracing, this is the starting point for the Ray.
   * During tracing, this is the furthest ahead that a Ray has traced. For example,
   * when on a segment in a RayKernel, this will be end of said segment.
   * After tracing, this is where the Ray ended.
   */
  Point _current_point;

  /// Direction of the Ray
  Point _direction;

  /**
   * Current element that the Ray is in.
   *
   * Before tracing, this is the starting element for the Ray.
   *
   * During tracing:
   * - It is valid within RayKernels, because a Ray can only operate on a
   *   single element per segment.
   * - When used on boundaries (RayBCs), it is the element that the
   *   Ray actually traced through. When on a boundary, a RayBC may be
   *   applied to multiple elements when hitting a vertex or edge.
   *   Therefore, this will be only one of said elements.
   */
  const Elem * _current_elem;

  /**
   * The side of _current_elem that the Ray is incoming on (if any).
   *
   * Before tracing, this is the starting incoming side (if any).
   *
   * During tracing, this is ONLY guaranteed to be valid during
   * the execution of RayKernels!
   */
  unsigned short _current_incoming_side;

  /// Whether or not the user has set an end point for this Ray (via limiting its distance)
  bool _end_set;

  /// Number of times this Ray has been communicated
  unsigned int _processor_crossings;

  /// Number of intersections done for this Ray
  unsigned int _intersections;

  /// Number of times this Ray has had its trajectory changed
  unsigned int _trajectory_changes;
  /// Whether or not this Ray had its trajectory changed (not sent in parallel)
  bool _trajectory_changed;

  /// Total distance this Ray has traveled
  Real _distance;
  /// Maximum distance the Ray is allowed to travel
  Real _max_distance;

  /// Wether or not the Ray should continue to be traced (not sent in parallel!)
  bool _should_continue;

  /// The data that is carried with the Ray
  std::vector<RayData> _data;

  /// Auxilary data that is carried with the ray
  std::vector<RayData> _aux_data;

  /// Extra padding to avoid false sharing
  long padding[8];

  // TraceRay is the only object that should be executing Rays and therefore needs access
  friend class TraceRay;
  // Packing needs access to changing the internal counters during the trace
  friend class Parallel::Packing<std::shared_ptr<Ray>>;
};

/**
 * The following methods are specializations for using the Parallel::packed_range_* routines
 * for a vector of Rays
 */
namespace libMesh
{
namespace Parallel
{
template <>
class Packing<std::shared_ptr<Ray>>
{
public:
  typedef Real buffer_type;

  static unsigned int mixed_size;
  static unsigned int packed_size(typename std::vector<Real>::const_iterator in);
  static unsigned int packable_size(const std::shared_ptr<Ray> & ray, const void *);
  static unsigned int size(const std::size_t data_size, const std::size_t aux_data_size);

  template <typename Iter, typename Context>
  static void pack(const std::shared_ptr<Ray> & object, Iter data_out, const Context * context);

  template <typename BufferIter, typename Context>
  static std::shared_ptr<Ray> unpack(BufferIter in, Context * context);
};

} // namespace Parallel

} // namespace libMesh
