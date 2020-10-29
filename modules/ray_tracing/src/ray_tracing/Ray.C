//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Ray.h"

// Local includes
#include "RayTracingStudy.h"
#include "RayTracingPackingUtils.h"

Ray::Ray()
  : _id(INVALID_RAY_ID),
    _current_point(RayTracingCommon::invalid_point),
    _direction(RayTracingCommon::invalid_point),
    _current_elem(nullptr),
    _current_incoming_side(RayTracingCommon::invalid_side),
    _end_set(false),
    _processor_crossings(0),
    _intersections(0),
    _trajectory_changes(0),
    _trajectory_changed(false),
    _distance(0),
    _max_distance(std::numeric_limits<Real>::max()),
    _should_continue(true),
    _data(),
    _aux_data()
{
}

void
Ray::reset()
{
  resetCounters();
  _id = INVALID_RAY_ID;
  _current_point = RayTracingCommon::invalid_point;
  _direction = RayTracingCommon::invalid_point;
  _current_elem = nullptr;
  _current_incoming_side = RayTracingCommon::invalid_side;
  _end_set = false;
  _max_distance = std::numeric_limits<Real>::max();
  _data.clear();
  _aux_data.clear();
}

Ray::Ray(const std::size_t data_size, const std::size_t aux_data_size)
  : _id(INVALID_RAY_ID),
    _current_point(RayTracingCommon::invalid_point),
    _direction(RayTracingCommon::invalid_point),
    _current_elem(nullptr),
    _current_incoming_side(RayTracingCommon::invalid_side),
    _end_set(false),
    _max_distance(std::numeric_limits<Real>::max()),
    _should_continue(true),
    _data(data_size, 0),
    _aux_data(aux_data_size, 0)
{
  resetCounters();
}

void
Ray::reset(const std::size_t data_size, const std::size_t aux_data_size)
{
  resetCounters();
  _id = INVALID_RAY_ID;
  _current_point = RayTracingCommon::invalid_point;
  _direction = RayTracingCommon::invalid_point;
  _current_elem = nullptr;
  _current_incoming_side = RayTracingCommon::invalid_side;
  _end_set = false;
  _max_distance = std::numeric_limits<Real>::max();
  _data.resize(data_size);
  std::fill(_data.begin(), _data.end(), 0);
  _aux_data.resize(aux_data_size);
  std::fill(_aux_data.begin(), _aux_data.end(), 0);
}

Ray::Ray(const std::size_t data_size, const std::size_t aux_data_size, const NonResetCountersKey &)
  : _data(data_size), _aux_data(aux_data_size)
{
}

Ray::Ray(const Ray & other)
  : _id(other._id),
    _current_point(other._current_point),
    _direction(other._direction),
    _current_elem(other._current_elem),
    _current_incoming_side(other._current_incoming_side),
    _end_set(other._end_set),
    _processor_crossings(0),
    _intersections(0),
    _trajectory_changes(0),
    _trajectory_changed(false),
    _distance(0),
    _max_distance(other._max_distance),
    _should_continue(true),
    _data(other._data),
    _aux_data(other._aux_data)
{
}

void
Ray::reset(const Ray & other)
{
  *this = other;
}

Ray &
Ray::operator=(const Ray & other)
{
  resetCounters();

  _id = other._id;
  _current_point = other._current_point;
  _direction = other._direction;
  _current_elem = other._current_elem;
  _current_incoming_side = other._current_incoming_side;
  _end_set = other._end_set;
  _max_distance = other._max_distance;
  _data = other._data;
  _aux_data = other._aux_data;

  return *this;
}

bool
Ray::operator==(const Ray & other) const
{
  return equalityHelper(other, true);
}

bool
Ray::operator!=(const Ray & other) const
{
  return equalityHelper(other, false);
}

bool
Ray::equalityHelper(const Ray & other, const bool equal) const
{
  if (this == &other)
    return equal;

  if (_id != other._id)
    return !equal;
  if (invalidCurrentPoint() != other.invalidCurrentPoint() ||
      !_current_point.absolute_fuzzy_equals(other._current_point))
    return !equal;
  if (invalidDirection() != other.invalidDirection() ||
      !_direction.absolute_fuzzy_equals(other._direction))
    return !equal;
  if (_current_elem != other._current_elem)
    return !equal;
  if (_current_incoming_side != other._current_incoming_side)
    return !equal;
  if (_end_set != other._end_set)
    return !equal;
  if (_processor_crossings != other._processor_crossings)
    return !equal;
  if (_intersections != other._intersections)
    return !equal;
  if (_trajectory_changes != other._trajectory_changes)
    return !equal;
  if (_trajectory_changed != other._trajectory_changed)
    return !equal;
  if (!MooseUtils::absoluteFuzzyEqual(_distance, other._distance))
    return !equal;
  if (!MooseUtils::absoluteFuzzyEqual(_max_distance, other._max_distance))
    return !equal;
  if (_should_continue != other._should_continue)
    return !equal;
  if (_data.size() != other._data.size())
    return !equal;
  for (std::size_t i = 0; i < _data.size(); ++i)
    if (!MooseUtils::absoluteFuzzyEqual(_data[i], other._data[i]))
      return !equal;
  if (_aux_data.size() != other._aux_data.size())
    return !equal;
  for (std::size_t i = 0; i < _aux_data.size(); ++i)
    if (!MooseUtils::absoluteFuzzyEqual(_aux_data[i], other._aux_data[i]))
      return !equal;

  return equal;
}

void
Ray::reset(const std::size_t data_size,
           const std::size_t aux_data_size,
           const NonResetCountersKey &)
{
  _data.resize(data_size);
  _aux_data.resize(aux_data_size);
}

void
Ray::setAuxData(const std::size_t start_i, const std::vector<RayData> & values)
{
  mooseAssert(start_i + values.size() <= _aux_data.size(), "Ray aux data size too small");
  for (std::size_t i = 0; i < values.size(); ++i)
    _aux_data[start_i + i] = values[i];
}

void
Ray::setStartEnd(const Point & start, const Point & end)
{
  errorIfTracing("Cannot set a Ray's start/end");

  if (start.absolute_fuzzy_equals(end))
    mooseError("Cannot set start/end for a Ray when start == end\n\n", getInfo());

  const auto distance = (end - start).norm();
  _current_point = start;
  _max_distance = distance;
  _direction = (end - start) / distance;
  _end_set = true;
}

void
Ray::setStartDirection(const Point & start, const Point & direction)
{
  errorIfTracing("Cannot set a Ray's start/direction");

  if (direction.absolute_fuzzy_equals(Point(0, 0, 0)))
    mooseError("Cannot set zero vector direction for a Ray\n\n", getInfo());

  _current_point = start;
  _direction = direction.unit();
  _end_set = false;
}

bool
Ray::atEnd() const
{
  if (!_end_set)
    mooseError("Cannot use Ray::atEnd() when the Ray trajectory was not set with setStartEnd().");

  return MooseUtils::absoluteFuzzyEqual(_distance, _max_distance);
}

void
Ray::changeDirection(const Point & direction, const ChangeDirectionKey)
{
  if (direction.absolute_fuzzy_equals(Point(0, 0, 0)))
    mooseError("Cannot set zero vector direction for a Ray\n\n", getInfo());

  _direction = direction.unit();
  _trajectory_changed = true;
}

void
Ray::changeStartDirection(const Point & start,
                          const Point & direction,
                          const ChangeStartDirectionKey)
{
  if (direction.absolute_fuzzy_equals(Point(0, 0, 0)))
    mooseError("Cannot set zero vector direction for a Ray\n\n", getInfo());

  _current_point = start;
  _direction = direction.unit();
  _trajectory_changed = true;
}

const Point &
Ray::startPoint() const
{
  errorIfTracing("Cannot get a Ray's start point");
  return _current_point;
}

bool
Ray::invalidStartPoint() const
{
  errorIfTracing("Cannot see if a Ray's start point is invalid");
  return _current_point == RayTracingCommon::invalid_point;
}

void
Ray::setStartingElem(const Elem * elem)
{
  errorIfTracing("Cannot set a Ray's starting element");
  _current_elem = elem;
}

const Elem *
Ray::startingElem() const
{
  errorIfTracing("Cannot get a Ray's starting elem");
  return _current_elem;
}

void
Ray::invalidateStartingElem()
{
  errorIfTracing("Cannot invalidate a Ray's starting element");
  _current_elem = nullptr;
  _current_incoming_side = RayTracingCommon::invalid_side;
}

void
Ray::setStartingIncomingSide(const unsigned short side)
{
  errorIfTracing("Cannot set a Ray's starting incoming side");
  _current_incoming_side = side;
}

unsigned short
Ray::startingIncomingSide() const
{
  errorIfTracing("Cannot get a Ray's starting incoming side");
  return _current_incoming_side;
}

bool
Ray::invalidStartingIncomingSide() const
{
  errorIfTracing("Cannot get a Ray's starting incoming side");
  return _current_incoming_side == RayTracingCommon::invalid_side;
}

void
Ray::invalidateStartingIncomingSide()
{
  errorIfTracing("Cannot invalidate a Ray's starting incoming side");
  _current_incoming_side = RayTracingCommon::invalid_side;
}

void
Ray::setMaxDistance(const Real max_distance)
{
  errorIfTracing("Cannot set a Ray's max distance");

  if (_end_set)
    mooseError("Cannot set the max distance for a Ray after calling setStartEnd().",
               "\nThe max distance is set internally to account for the user-set end point.\n\n",
               getInfo());

  _max_distance = max_distance;
}

void
Ray::errorIfTracing(const std::string & reason) const
{
  if (hasTraced())
    mooseError(reason, " after it has started tracing\n\n", getInfo());
}

void
Ray::resetCounters()
{
  _processor_crossings = 0;
  _intersections = 0;
  _trajectory_changes = 0;
  _distance = 0;
  _trajectory_changed = false;
  _should_continue = true;
}

void
Ray::resetCounters(const ResetCountersKey)
{
  resetCounters();
}

void
Ray::setID(const RayID id)
{
  errorIfTracing("Cannot set a Ray's ID");
  _id = id;
}

void
Ray::invalidateID()
{
  errorIfTracing("Cannot invalidate a Ray's ID");
  _id = INVALID_RAY_ID;
}

std::string
Ray::getInfo(const RayTracingStudy * study /* = nullptr */) const
{
  std::ostringstream oss;

  oss << "Ray information";
  if (study)
  {
    oss << " with " << study->type() << " '" << study->name() << "'";
    oss << " on pid " << study->comm().rank();
  }
  oss << "\n";
  oss << "  this = " << this << "\n";
  oss << "  id() = ";
  if (invalidID())
    oss << "INVALID_RAY_ID\n";
  else
    oss << id() << "\n";
  if (study && study->useRayRegistration() && !invalidID())
    oss << "  study->registeredRayName(id()) = " << study->registeredRayName(id()) << "\n";
  oss << "  currentPoint() = ";
  if (invalidCurrentPoint())
    oss << "invalid point\n";
  else
    oss << currentPoint() << "\n";
  oss << "  direction() = " << direction() << "\n";
  oss << "  currentIncomingSide() = ";
  if (invalidCurrentIncomingSide())
    oss << "invalid side\n";
  else
    oss << currentIncomingSide() << "\n";
  oss << "  endSet() = " << (endSet() ? "true" : "false") << "\n";
  if (endSet())
    oss << "  atEnd() = " << (atEnd() ? "true" : "false") << "\n";
  oss << "  distance() = " << distance() << "\n";
  oss << "  maxDistance() = " << maxDistance() << "\n";
  if (currentElem())
  {
    oss << "  currentElem()->id() = " << currentElem()->id() << "\n";
    oss << "  currentElem()->processor_id() = " << currentElem()->processor_id() << "\n";
  }
  else
    oss << "  currentElem()->id() = invalid\n";
  oss << "  processorCrossings() = " << processorCrossings() << "\n";
  oss << "  intersections() = " << intersections() << "\n";
  oss << "  trajectoryChanges() = " << trajectoryChanges() << "\n";
  oss << "  shouldContinue() = " << (shouldContinue() ? "true" : "false") << "\n";
  oss << "  trajectoryChanged() = " << (trajectoryChanged() ? "true" : "false") << "\n";
  oss << "  data() = ";
  if (data().empty())
    oss << "empty";
  if (!study || study->rayDataSize() != data().size())
  {
    for (const auto & entry : data())
      oss << entry << " ";
  }
  else
  {
    const auto & names = study->rayDataNames();
    for (std::size_t i = 0; i < data().size(); ++i)
      oss << "\n    '" << names[i] << "' = " << data(i);
  }
  oss << "\n";
  oss << "  auxData() = ";
  if (auxData().empty())
    oss << "empty";
  if (!study || study->rayAuxDataSize() != auxData().size())
  {
    for (const auto & entry : auxData())
      oss << entry << " ";
  }
  else
  {
    const auto & names = study->rayAuxDataNames();
    for (std::size_t i = 0; i < auxData().size(); ++i)
      oss << "\n    '" << names[i] << "' = " << auxData(i);
  }
  oss << "\n";

  return oss.str();
}

namespace libMesh
{
namespace Parallel
{

unsigned int Packing<std::shared_ptr<Ray>>::mixed_size =
    RayTracingPackingUtils::mixedPackSize<buffer_type>(
        (unsigned short)0, (bool)true, (unsigned int)0, (unsigned int)0, (unsigned int)0);

unsigned int
Packing<std::shared_ptr<Ray>>::size(const std::size_t data_size, const std::size_t aux_data_size)
{
  // Data lengths, current point and direction, current element id, distance, max_distance, ray ID
  unsigned int size = 12;
  // Current incoming side, end_set, processor crossings, intersections, trajectory changes
  // (packed into as few buffer_type as possible)
  size += mixed_size;

#ifdef SINGLE_PRECISION_RAY
  if (data_size)
    size += RayTracingPackingUtils::reinterpretCopySize<RayData, buffer_type>(data_size);
  if (aux_data_size)
    size += RayTracingPackingUtils::reinterpretCopySize<RayData, buffer_type>(aux_data_size);
#else
  size += data_size + aux_data_size;
#endif

  return size;
}

unsigned int
Packing<std::shared_ptr<Ray>>::packed_size(typename std::vector<buffer_type>::const_iterator in)
{
  const std::size_t data_size = static_cast<std::size_t>(*in++);
  const std::size_t aux_data_size = static_cast<std::size_t>(*in++);

  return size(data_size, aux_data_size);
}

unsigned int
Packing<std::shared_ptr<Ray>>::packable_size(const std::shared_ptr<Ray> & ray, const void *)
{
  return size(ray->data().size(), ray->auxData().size());
}

template <>
std::shared_ptr<Ray>
Packing<std::shared_ptr<Ray>>::unpack(std::vector<buffer_type>::const_iterator in,
                                      RayTracingStudy * study)
{
  // Grab the data size
  const std::size_t data_size = static_cast<std::size_t>(*in++);
  const std::size_t aux_data_size = static_cast<std::size_t>(*in++);

  // Needed to access the Ray constructor/reset that doesn't reset the counters
  const Ray::NonResetCountersKey key;
  std::shared_ptr<Ray> ray = study->acquireRay(/* tid = */ 0, data_size, aux_data_size, key);

  // Current Point
  ray->_current_point(0) = *in++;
  ray->_current_point(1) = *in++;
  ray->_current_point(2) = *in++;

  // Direction
  ray->_direction(0) = *in++;
  ray->_direction(1) = *in++;
  ray->_direction(2) = *in++;

  // Current Element
  RayTracingPackingUtils::unpack(ray->_current_elem, *in++, &study->meshBase());

  // Current incoming size, end set, processor crossings, intersections, trajectory changes
  // (unpaced from as few buffer_type as possible - 5 values from 2 Reals)
  RayTracingPackingUtils::mixedUnpack<buffer_type>(in,
                                                   ray->_current_incoming_side,
                                                   ray->_end_set,
                                                   ray->_processor_crossings,
                                                   ray->_intersections,
                                                   ray->_trajectory_changes);

  // Distance
  ray->_distance = *in++;

  // Max distance
  ray->_max_distance = *in++;

  // ID
  RayTracingPackingUtils::unpack(*in++, ray->_id);

#ifdef SINGLE_PRECISION_RAY
  RayTracingPackingUtils::reinterpretUnpackCopy<buffer_type>(ray->_data, in);
  RayTracingPackingUtils::reinterpretUnpackCopy<buffer_type>(ray->_aux_data, in);
#else
  // Copy out data
  RayTracingPackingUtils::unpackCopy(ray->_data, in);
  RayTracingPackingUtils::unpackCopy(ray->_aux_data, in);
#endif

  ray->_should_continue = true;
  ray->_trajectory_changed = false;

  return ray;
}

template <>
void
Packing<std::shared_ptr<Ray>>::pack(const std::shared_ptr<Ray> & ray,
                                    std::back_insert_iterator<std::vector<buffer_type>> data_out,
                                    const RayTracingStudy * study)
{
  // Storing the data size first makes it easy to verify and reserve space
  data_out = static_cast<buffer_type>(ray->_data.size());
  data_out = static_cast<buffer_type>(ray->_aux_data.size());

  // Current Point
  data_out = ray->_current_point(0);
  data_out = ray->_current_point(1);
  data_out = ray->_current_point(2);

  // Direction
  data_out = ray->_direction(0);
  data_out = ray->_direction(1);
  data_out = ray->_direction(2);

  // Current element
  data_out = RayTracingPackingUtils::pack<buffer_type>(ray->_current_elem, &study->meshBase());

  // Current incoming size, end set, processor crossings, intersections, trajectory changes
  // (packed into as few buffer_type as possible - 5 values into 2 Reals
  RayTracingPackingUtils::mixedPack<buffer_type>(data_out,
                                                 ray->_current_incoming_side,
                                                 ray->_end_set,
                                                 ray->_processor_crossings,
                                                 ray->_intersections,
                                                 ray->_trajectory_changes);

  // Distance
  data_out = ray->_distance;
  // Max distance
  data_out = ray->_max_distance;

  // ID
  data_out = RayTracingPackingUtils::pack<buffer_type>(ray->id());

  // Copy out data
#ifdef SINGLE_PRECISION_RAY
  RayTracingPackingUtils::reinterpretPackCopy<buffer_type>(ray->_data, data_out);
  RayTracingPackingUtils::reinterpretPackCopy<buffer_type>(ray->_aux_data, data_out);
#else
  std::copy(ray->_data.begin(), ray->_data.end(), data_out);
  std::copy(ray->_aux_data.begin(), ray->_aux_data.end(), data_out);
#endif
}

} // namespace Parallel

} // namespace libMesh
