//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ParallelStudyMethod.h"

// MOOSE includes
#include "InputParameters.h"
#include "SharedPool.h"
#include "MooseEnum.h"
#include "CircularBuffer.h"
#include "LIFOBuffer.h"

// Local includes
#include "SendBuffer.h"
#include "ReceiveBuffer.h"

// libMesh Includes
#include "libmesh/parallel_object.h"

template <typename Object, typename Context>
class ParallelStudy : public ParallelObject
{
public:
  static InputParameters validParams();

  ParallelStudy(const libMesh::Parallel::Communicator & comm,
                Context & context,
                const InputParameters & params);

  /**
   * Pre-execute method that MUST be called before execute() and before adding objects.
   */
  void preExecute();
  /**
   * Execute method.
   */
  void execute();

  /**
   * Adds an object to the buffer to be executed. This will move the object into the buffer
   * (with std::move), therefore the object will be nullptr after this call.
   *
   * During pre-execution (between preExecute() and execute()), this method can ONLY
   * be called on thread 0.
   *
   * During execute(), this method is thread safe and can be used to add objects during execution.
   */
  void moveObjectToBuffer(std::shared_ptr<Object> & object, const THREAD_ID tid);
  /**
   * Adds objects to the buffer to be executed. This will move the objects into the buffer
   * (with std::move), therefore all valid objects in \p objects will be nullptr after this call.
   *
   * During pre-execution (between preExecute() and execute()), this method can ONLY
   * be called on thread 0.
   *
   * During execute(), this method is thread safe and can be used to add objects during execution.
   */
  void moveObjectsToBuffer(std::vector<std::shared_ptr<Object>> & objects, const THREAD_ID tid);

  /**
   * Acquire an object from the pool.
   */
  template <typename... Args>
  typename MooseUtils::SharedPool<Object>::PtrType acquireObject(const THREAD_ID tid,
                                                                 Args &&... args)
  {
    return _threaded_object_pools[tid].acquire(std::forward<Args>(args)...);
  }

  /**
   * Gets the receive buffer.
   */
  const ReceiveBuffer<Object, Context> & receiveBuffer() const { return *_receive_buffer; }
  /**
   * Gets the threaded during execution objects
   */
  const std::vector<std::vector<std::shared_ptr<Object>>> & threadedDuringExecutionObjects() const
  {
    return _threaded_during_execution_objects;
  }
  /**
   * Gets the working buffer.
   */
  const MooseUtils::Buffer<std::shared_ptr<Object>> & workingObjectBuffer() const
  {
    return *_working_object_buffer;
  }

  /**
   * Gets tht total number of send buffer pools created.
   */
  unsigned long long int sendBufferPoolCreated() const;
  /**
   * Gets the total number of objects sent from this processor.
   */
  unsigned long long int objectsSent() const;
  /**
   * Gets the total number of buffers sent from this processor.
   */
  unsigned long long int buffersSent() const;
  /**
   * Gets the total number of objects created in all of the threaded pools
   */
  unsigned long long int poolObjectsCreated() const;

  /**
   * Gets the total number of objects started from this processor.
   */
  unsigned long long int localObjectsStarted() const { return _local_objects_started; }
  /**
   * Gets the total number of objects finished
   */
  unsigned long long int totalObjectsFinished() const { return _total_objects_finished; }
  /**
   * Gets the total number of objects executed on this processor.
   */
  unsigned long long int localObjectsExecuted() const { return _local_objects_executed; }
  /**
   * Gets the total number of chunks executed on this processor.
   */
  unsigned long long int localChunksExecuted() const { return _local_chunks_executed; }

  /**
   * Whether or not this object is currently in execute().
   */
  bool currentlyExecuting() const { return _currently_executing; }
  /**
   * Whether or not this object is between preExecute() and execute().
   */
  bool currentlyPreExecuting() const { return _currently_pre_executing; }

  /**
   * Gets the max buffer size
   */
  unsigned int maxBufferSize() const { return _max_buffer_size; }
  /**
   * Gets the chunk size
   */
  unsigned int chunkSize() const { return _chunk_size; }

  /**
   * Gets the number of iterations to wait before communicating
   */
  unsigned int clicksPerCommunication() const { return _clicks_per_communication; }
  /**
   * Gets the number of iterations to wait before communicating with root
   */
  unsigned int clicksPerRootCommunication() const { return _clicks_per_root_communication; }
  /**
   * Gets the number of iterations to wait before checking for new objects
   */
  unsigned int clicksPerReceive() const { return _clicks_per_receive; }

  /**
   * Gets the method
   */
  ParallelStudyMethod method() const { return _method; }

  /**
   * Reserve \p size entries in the buffer.
   *
   * This can only be used during the pre-execution phase (between preExecute() and execute()).
   *
   * This is particularly useful when one wants to move many objects into the buffer using
   * moveObjectToBuffer() and wants to allocate the space ahead of time.
   */
  void reserveBuffer(const std::size_t size);

protected:
  /**
   * Enum for providing useful errors during object addition in addObjectError().
   */
  enum AddObjectError
  {
    ADDITION_DURING_EXECUTION_DISABLED,
    PRE_EXECUTION_AND_EXECUTION_ONLY,
    PRE_EXECUTION_ONLY,
    PRE_EXECUTION_THREAD_0_ONLY
  };

  /**
   * Pure virtual to be overridden that executes a single object on a given thread
   */
  virtual void executeObject(const std::shared_ptr<Object> & object, const THREAD_ID tid) = 0;
  /**
   * Pure virtual to be overridden that is called after an object is completed (when it is executed
   * and its objectProcessorId() does not change)
   */
  virtual void onCompleteObject(const std::shared_ptr<Object> & object) = 0;
  /**
   * Pure virtual to be overridden that specifies whether or not an object should be executed.
   *
   * If this is false, there still exists an opportunity to communicate an object (without doing
   * work on it on this processor) with shouldCommunicateObject().
   * If this is true, executeObject() will be called on the object.
   */
  virtual bool shouldExecuteObject(const std::shared_ptr<Object> & object) const = 0;
  /**
   * Pure virtual to be overriden that specifices whether or not an object should be communicated.
   * \p to_pid is to be filled with the processor that the object should be communicated to.
   *
   * If this returns false, it is implied that the object is done with execution. onCompleteObject()
   * will be called on the object and ownership of the object by the working buffer is released.
   */
  virtual bool shouldCommunicateObject(const std::shared_ptr<Object> & object,
                                       processor_id_type & to_pid) const = 0;

  /**
   * Virtual that allows for the customiation of error text for adding objects.
   */
  virtual std::string addObjectError(const std::shared_ptr<Object> & object,
                                     const AddObjectError error) const;

  /// The context
  Context & _context;
  /// This rank
  const processor_id_type _pid;
  /// Name for this object for use in error handling
  const std::string _name;
  /// The study method
  const ParallelStudyMethod _method;

private:
  /**
   * Flushes all objects out of the send buffers
   */
  void flushSendBuffers();

  /**
   * Execute objects using SMART
   */
  void smartExecute();
  /**
   * Execute objects using HARM
   */
  void harmExecute();
  /**
   * Execute objects using BS
   */
  void bsExecute();

  /**
   * Reeive packets of objects from other processors and execute them
   */
  bool receiveAndExecute();

  /**
   * Execute a chunk of objects and buffer
   */
  void executeAndBuffer(const std::size_t chunk_size);

  /// Minimum size of a SendBuffer
  const unsigned int _min_buffer_size;
  /// Number of objects to buffer before communication
  const unsigned int _max_buffer_size;
  /// Multiplier for the buffer size for growing the buffer
  const Real _buffer_growth_multiplier;
  /// Multiplier for the buffer size for shrinking the buffer
  const Real _buffer_shrink_multiplier;
  /// Number of objects to execute at once during communication
  const unsigned int _chunk_size;
  /// Whether or not to allow the addition of objects during the execution phase
  const bool _allow_addition_during_execution;

  /// Iterations to wait before communicating
  const unsigned int _clicks_per_communication;
  /// Iterations to wait before communicating with root
  const unsigned int _clicks_per_root_communication;
  /// Iterations to wait before checking for new objects
  const unsigned int _clicks_per_receive;

  /// MessageTag for sending objects
  Parallel::MessageTag _object_buffer_tag;
  /// Pools for re-using destructed objects
  std::vector<MooseUtils::SharedPool<Object>> _threaded_object_pools;
  /// Threaded temprorary storage for objects added during execution
  std::vector<std::vector<std::shared_ptr<Object>>> _threaded_during_execution_objects;
  /// Buffer for executing objects
  std::unique_ptr<MooseUtils::Buffer<std::shared_ptr<Object>>> _working_object_buffer;
  /// The receive buffer
  std::unique_ptr<ReceiveBuffer<Object, Context>> _receive_buffer;
  /// Send buffers for each processor
  std::unordered_map<processor_id_type, std::unique_ptr<SendBuffer<Object, Context>>> _send_buffers;

  /// Number of chunks executed on this processor
  unsigned long long int _local_chunks_executed;
  /// Number of objects executed on this processor
  unsigned long long int _local_objects_executed;
  /// Number of objects started on this processor
  unsigned long long int _local_objects_started;
  /// Number of objects started everywhere
  unsigned long long int _total_objects_started;
  /// Number of objects finished on this processor
  unsigned long long int _local_objects_finished;
  /// Number of object finished everywhere
  unsigned long long int _total_objects_finished;

  /// Whether or not this object is currently within execute()
  bool _currently_executing;
  /// Whether or not this object is between preExecute() and execute()
  bool _currently_pre_executing;
};

template <typename Object, typename Context>
ParallelStudy<Object, Context>::ParallelStudy(const libMesh::Parallel::Communicator & comm,
                                              Context & context,
                                              const InputParameters & params)
  : ParallelObject(comm),
    _context(context),
    _pid(comm.rank()),
    _name("temp_name"),

    _method((ParallelStudyMethod)(int)(params.get<MooseEnum>("method"))),
    _min_buffer_size(params.isParamSetByUser("min_buffer_size")
                         ? params.get<unsigned int>("min_buffer_size")
                         : params.get<unsigned int>("send_buffer_size")),
    _max_buffer_size(params.get<unsigned int>("send_buffer_size")),
    _buffer_growth_multiplier(params.get<Real>("buffer_growth_multiplier")),
    _buffer_shrink_multiplier(params.get<Real>("buffer_shrink_multiplier")),
    _chunk_size(params.get<unsigned int>("chunk_size")),
    _allow_addition_during_execution(params.get<bool>("allow_addition_during_execution")),

    _clicks_per_communication(params.get<unsigned int>("clicks_per_communication")),
    _clicks_per_root_communication(params.get<unsigned int>("clicks_per_root_communication")),
    _clicks_per_receive(params.get<unsigned int>("clicks_per_receive")),

    _object_buffer_tag(Parallel::MessageTag(100)),
    _threaded_object_pools(libMesh::n_threads()),
    _threaded_during_execution_objects(_allow_addition_during_execution ? libMesh::n_threads() : 0),

    _currently_executing(false),
    _currently_pre_executing(false)
{
  // Create the work buffer
  const auto buffer_type = params.get<MooseEnum>("work_buffer_type");
  if (buffer_type == "lifo")
    _working_object_buffer = libmesh_make_unique<MooseUtils::LIFOBuffer<std::shared_ptr<Object>>>();
  else if (buffer_type == "circular")
    _working_object_buffer =
        libmesh_make_unique<MooseUtils::CircularBuffer<std::shared_ptr<Object>>>();
  else
    mooseError("Unknown work buffer type");

  // Create the receive buffer
  _receive_buffer = libmesh_make_unique<ReceiveBuffer<Object, Context>>(
      comm, context, *_working_object_buffer, _method, _clicks_per_receive, _object_buffer_tag);

  if (_method != ParallelStudyMethod::SMART && _allow_addition_during_execution)
    mooseError(_name,
               ": When allowing object addition during "
               "execution\n('allow_addition_during_execution = true'), the method must be SMART");

#ifndef LIBMESH_HAVE_OPENMP
  if (libMesh::n_threads() != 1)
    mooseWarning(_name, ": Threading will not be used without OpenMP");
#endif
}

template <typename Object, typename Context>
InputParameters
ParallelStudy<Object, Context>::validParams()
{
  auto params = emptyInputParameters();

  params.addRangeCheckedParam<unsigned int>(
      "send_buffer_size", 100, "send_buffer_size > 0", "The size of the send buffer");
  params.addRangeCheckedParam<unsigned int>(
      "chunk_size",
      100,
      "chunk_size > 0",
      "The number of objects to process at one time during execution");
  params.addRangeCheckedParam<unsigned int>("clicks_per_communication",
                                            10,
                                            "clicks_per_communication > 0",
                                            "Iterations to wait before communicating");
  params.addRangeCheckedParam<unsigned int>("clicks_per_root_communication",
                                            10,
                                            "clicks_per_root_communication > 0",
                                            "Iterations to wait before communicating with root");
  params.addRangeCheckedParam<unsigned int>("clicks_per_receive",
                                            1,
                                            "clicks_per_receive > 0",
                                            "Iterations to wait before checking for new objects");

  params.addParam<unsigned int>("min_buffer_size",
                                "The initial size of the SendBuffer and the floor for shrinking "
                                "it.  This defaults to send_buffer_size if not set (i.e. the "
                                "buffer won't change size)");
  params.addParam<Real>("buffer_growth_multiplier",
                        2.,
                        "How much to grow a SendBuffer by if the buffer completely fills and "
                        "dumps.  Will max at send_buffer_size");
  params.addRangeCheckedParam<Real>("buffer_shrink_multiplier",
                                    0.5,
                                    "0 < buffer_shrink_multiplier <= 1.0",
                                    "Multiplier (between 0 and 1) to apply to the current buffer "
                                    "size if it is force dumped.  Will stop at "
                                    "min_buffer_size.");

  params.addParam<bool>(
      "allow_addition_during_execution",
      false,
      "Whether or not to allow the addition of objects to the buffer during execution");

  MooseEnum methods("smart harm bs", "smart");
  params.addParam<MooseEnum>("method", methods, "The algorithm to use");

  MooseEnum work_buffers("lifo circular", "circular");
  params.addParam<MooseEnum>("work_buffer_type", work_buffers, "The work buffer type to use");

  return params;
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::executeAndBuffer(const std::size_t chunk_size)
{
  // If chunk_size > the number of objects left, this will properly grab all of them
  const auto begin = _working_object_buffer->beginChunk(chunk_size);
  const auto end = _working_object_buffer->endChunk(chunk_size);

  _local_chunks_executed++;

#ifdef LIBMESH_HAVE_OPENMP
#pragma omp parallel
#endif
  {
    const THREAD_ID tid =
#ifdef LIBMESH_HAVE_OPENMP
        omp_get_thread_num();
#else
        0;
#endif

#ifdef LIBMESH_HAVE_OPENMP
#pragma omp for schedule(dynamic, 20) nowait
#endif
    for (auto it = begin; it < end; ++it)
    {
      const std::shared_ptr<Object> & object = *it;

      if (!object || !shouldExecuteObject(object))
        continue;

      executeObject(object, tid);
      ++_local_objects_executed;
    }
  }

  processor_id_type to_pid;
  for (auto it = begin; it < end; ++it)
  {
    std::shared_ptr<Object> & object = *it;

    if (!object)
      continue;

    // Communicate object to another processor
    if (shouldCommunicateObject(object, to_pid))
    {
      if (comm().size() < to_pid)
        mooseError(_name, ": Object to be communicated has invalid pid");
      if (to_pid == _pid)
        mooseError(_name, ": Object to be communicated has same pid");

      // Get the send buffer for the proc this object is going to
      auto find_pair = _send_buffers.find(to_pid);
      // Need to create a send buffer for said processor
      if (find_pair == _send_buffers.end())
        _send_buffers
            .emplace(to_pid,
                     libmesh_make_unique<SendBuffer<Object, Context>>(comm(),
                                                                      _context,
                                                                      to_pid,
                                                                      _method,
                                                                      _min_buffer_size,
                                                                      _max_buffer_size,
                                                                      _buffer_growth_multiplier,
                                                                      _buffer_shrink_multiplier,
                                                                      _object_buffer_tag))
            .first->second->addObject(object);
      // Send buffer exists for this processor
      else
        find_pair->second->addObject(object);
    }
    // We've done work on this object and it's still on this processor: it's done
    else
    {
      ++_local_objects_finished;
      onCompleteObject(object);
    }

    // Remove the working buffer's ownership (use count) of this object
    object.reset();
  }

  // Remove the objects we just worked on from the buffer
  _working_object_buffer->eraseChunk(chunk_size);

  // Add any objects that were added during execution (added to _threaded_during_execution_objects)
  // to the working buffer to be executed
  if (_allow_addition_during_execution)
  {
    // Allocate space ahead of time
    std::size_t num_objects = 0;
    for (const auto & objects : _threaded_during_execution_objects)
      for (const std::shared_ptr<Object> & object : objects)
        if (object)
          ++num_objects;

    if (num_objects)
    {
      // We don't ever want to decrease the capacity, so only set it if we need more entries
      if (_working_object_buffer->capacity() < _working_object_buffer->size() + num_objects)
        _working_object_buffer->setCapacity(_working_object_buffer->size() + num_objects);

      // Move the objects
      for (auto & objects : _threaded_during_execution_objects)
      {
        for (std::shared_ptr<Object> & object : objects)
          if (object)
          {
            _working_object_buffer->move(object);
            mooseAssert(!object, "Object was not moved");
          }

        objects.clear();
      }

      // Variable that must be set when adding objects so that the algorithm can keep count
      // of how many objects still need to be executed
      _local_objects_started += num_objects;
    }
  }

  if (_method == ParallelStudyMethod::HARM)
    flushSendBuffers();
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::flushSendBuffers()
{
  for (auto & send_buffer_iter : _send_buffers)
    send_buffer_iter.second->forceSend();
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::reserveBuffer(const std::size_t size)
{
  if (!_currently_pre_executing)
    mooseError(_name, ": Can only reserve in object buffer during pre-execution");

  // We don't ever want to decrease the capacity, so only set if we need more entries
  if (_working_object_buffer->capacity() < size)
    _working_object_buffer->setCapacity(size);
}

template <typename Object, typename Context>
bool
ParallelStudy<Object, Context>::receiveAndExecute()
{
  bool executed_some = false;

  if (_receive_buffer->currentlyReceiving() && _method == ParallelStudyMethod::SMART)
    _receive_buffer->cleanupRequests();
  else
    _receive_buffer->receive();

  while (!_working_object_buffer->empty())
  {
    executed_some = true;

    // Switch between tracing a chunk and buffering with SMART
    if (_method == ParallelStudyMethod::SMART)
    {
      // Look for extra work first so that these transfers can be finishing while we're executing
      // Start receives only if our work buffer is decently sized
      _receive_buffer->receive(_working_object_buffer->size() > (2 * _chunk_size));

      // Execute some objects
      executeAndBuffer(_chunk_size);
    }
    // Execute all of them and then buffer with the other methods
    else
      executeAndBuffer(_working_object_buffer->size());
  }

  return executed_some;
}

template <typename T>
void
nonblockingSum(const libMesh::Parallel::Communicator & comm,
               T & r,
               T & o,
               libMesh::Parallel::Request & req)
{
  if (comm.size() > 1)
    libmesh_call_mpi(MPI_Iallreduce(
        &r, &o, 1, libMesh::Parallel::StandardType<T>(&r), MPI_SUM, comm.get(), req.get()));
  else
    o = r;
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::smartExecute()
{
  mooseAssert(_method == ParallelStudyMethod::SMART, "Should be called with SMART only");

  // Request for the sum of the objects started
  Parallel::Request started_request;
  // Request for the sum of the objects finished
  Parallel::Request finished_request;

  // Temp for use in sending the current value in a nonblocking sum instead of an updated value
  unsigned long long int temp;

  // Whether or not to make the started request first, or after every finished request.
  // When allow adding objects during the execution phase, the starting object counts could
  // change after right now, so we must update them after each finished request is complete.
  // When not allowing generation during propagation, we know the counts up front.
  const bool started_request_first = !_allow_addition_during_execution;

  // Get the number of objects that were started in the whole domain, if applicable
  if (started_request_first)
    nonblockingSum(comm(), _local_objects_started, _total_objects_started, started_request);

  // Whether or not the started request has been made
  bool made_started_request = started_request_first;
  // Whether or not the finished request has been made
  bool made_finished_request = false;

  // Good time to get rid of whatever's currently in our SendBuffers
  flushSendBuffers();

  // Use these to try to delay some forced communication
  unsigned int non_executing_clicks = 0;
  unsigned int non_executing_root_clicks = 0;
  bool executed_some = true;

  // Keep executing the objects until they've all completed
  while (true)
  {
    executed_some = receiveAndExecute();

    if (executed_some)
    {
      non_executing_clicks = 0;
      non_executing_root_clicks = 0;
    }
    else
    {
      non_executing_clicks++;
      non_executing_root_clicks++;
    }

    if (non_executing_clicks >= _clicks_per_communication)
    {
      non_executing_clicks = 0;

      flushSendBuffers();
    }

    if (non_executing_root_clicks >= _clicks_per_root_communication)
    {
      non_executing_root_clicks = 0;

      // We need the starting object sum first but said request isn't complete yet
      if (started_request_first && !started_request.test())
        continue;

      // At this point, we need to make a request for the finished object sum
      if (!made_finished_request)
      {
        made_finished_request = true;
        temp = _local_objects_finished;
        nonblockingSum(comm(), temp, _total_objects_finished, finished_request);
        continue;
      }

      // We have the finished object sum
      if (finished_request.test())
      {
        // The starting object sum must be requested /after/ we have finishing counts and we
        // need to make the request for said sum
        if (!made_started_request)
        {
          made_started_request = true;
          temp = _local_objects_started;
          nonblockingSum(comm(), temp, _total_objects_started, started_request);
          continue;
        }

        // The starting object sum must be requested /after/ we have finishing sum and we
        // don't have the starting sum yet
        if (!started_request_first && !started_request.test())
          continue;

        // Started count is the same as the finished count - we're done!
        if (_total_objects_started == _total_objects_finished)
          return;

        // Next time around we should make a finished sum request
        made_finished_request = false;
        // If we need the starting sum after the finishing sum, we need those now as well
        if (!started_request_first)
          made_started_request = false;
      }
    }
  }
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::harmExecute()
{
  mooseAssert(_method == ParallelStudyMethod::HARM, "Should be called with HARM only");
  if (_allow_addition_during_execution)
    mooseError(_name, ": The addition of objects during execution is not supported by HARM");

  // Request for the total number of objects started
  Parallel::Request objects_started_request;
  // Requests for sending the number of finished objects to every other processor
  std::vector<Parallel::Request> objects_finished_requests(comm().size());
  // Whether or not the finished requests have been sent to each processor
  std::vector<bool> objects_finished_requests_sent(comm().size(), false);
  // Values of objects killed on this processor that are being sent to other processors
  std::vector<unsigned long long int> objects_finished_requests_temps(comm().size(), 0);
  // Objects finished by each processor
  std::vector<unsigned long long int> objects_finished_per_proc(comm().size(), 0);
  // Tag for sending objects finished
  Parallel::MessageTag objects_finished_requests_tag = Parallel::MessageTag(21000);

  // Get the number of objects that were started in the whole domain
  nonblockingSum(comm(), _local_objects_started, _total_objects_started, objects_started_request);

  // All objects have been executed, so time to communicate
  flushSendBuffers();

  // HARM only does some communication based on times through the loop.
  // This counter will be used for that
  unsigned int communication_clicks = 0;

  Parallel::Status objects_finished_probe_status;
  int objects_finished_probe_flag;

  // Keep bouncing the objects around until they've all completed
  while (true)
  {
    receiveAndExecute();

    flushSendBuffers();

    if (communication_clicks > comm().size())
    {
      // Receive messages about objects being finished
      do
      {
        MPI_Iprobe(MPI_ANY_SOURCE,
                   objects_finished_requests_tag.value(),
                   comm().get(),
                   &objects_finished_probe_flag,
                   objects_finished_probe_status.get());

        if (objects_finished_probe_flag)
        {
          auto proc = objects_finished_probe_status.source();
          comm().receive(proc, objects_finished_per_proc[proc], objects_finished_requests_tag);
        }
      } while (objects_finished_probe_flag);

      _total_objects_finished = std::accumulate(objects_finished_per_proc.begin(),
                                                objects_finished_per_proc.end(),
                                                _local_objects_finished);

      // Reset
      communication_clicks = 0;
    }

    // Send messages about objects being finished
    for (processor_id_type pid = 0; pid < comm().size(); ++pid)
      if (pid != _pid &&
          (!objects_finished_requests_sent[pid] || objects_finished_requests[pid].test()) &&
          _local_objects_finished > objects_finished_requests_temps[pid])
      {
        objects_finished_requests_temps[pid] = _local_objects_finished;
        comm().send(pid,
                    objects_finished_requests_temps[pid],
                    objects_finished_requests[pid],
                    objects_finished_requests_tag);
        objects_finished_requests_sent[pid] = true;
      }

    // All procs agree on the number of objects started and we've finished all the objects started
    if (objects_started_request.test() && _total_objects_started == _total_objects_finished)
    {
      // Need to call the post wait work for all of the requests
      for (processor_id_type pid = 0; pid < comm().size(); ++pid)
        if (pid != _pid)
          objects_finished_requests[pid].wait();

      return;
    }

    communication_clicks++;
  }
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::bsExecute()
{
  mooseAssert(_method == ParallelStudyMethod::BS, "Should be called with BS only");

  if (_allow_addition_during_execution)
    mooseError(_name, ": The addition of objects during execution is not supported by BS");

  Parallel::Request objects_started_request;
  Parallel::Request objects_finished_request;

  // Temp for use in sending the current value in a nonblocking sum instead of an updated value
  unsigned long long int temp;

  // Get the number of objects that were started in the whole domain
  nonblockingSum(comm(), _local_objects_started, _total_objects_started, objects_started_request);

  // Keep bouncing the objects around until they've all completed
  while (true)
  {
    bool receiving = false;
    bool sending = false;

    Parallel::Request some_left_request;
    unsigned int some_left = 0;
    unsigned int all_some_left = 1;

    do
    {
      _receive_buffer->receive();
      flushSendBuffers();

      receiving = _receive_buffer->currentlyReceiving();

      sending = false;
      for (auto & send_buffer : _send_buffers)
        sending = sending || send_buffer.second->currentlySending() ||
                  send_buffer.second->currentlyBuffered();

      if (!receiving && !sending && some_left_request.test() && all_some_left)
      {
        some_left = receiving || sending;

        nonblockingSum(comm(), some_left, all_some_left, some_left_request);
      }
    } while (receiving || sending || !some_left_request.test() || all_some_left);

    executeAndBuffer(_working_object_buffer->size());

    comm().barrier();

    if (objects_started_request.test() && objects_finished_request.test())
    {
      if (_total_objects_started == _total_objects_finished)
        return;

      temp = _local_objects_finished;
      nonblockingSum(comm(), temp, _total_objects_finished, objects_finished_request);
    }
  }
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::preExecute()
{
  // Buffers should be empty
  if (!_working_object_buffer->empty())
    mooseError(_name, ": Working buffer is not empty in preExecute()");
  mooseAssert(std::all_of(_working_object_buffer->data().begin(),
                          _working_object_buffer->data().end(),
                          [](const auto & object) { return !object; }),
              "Raw working buffer has objects in preExecute()");
  for (const auto & threaded_buffer : _threaded_during_execution_objects)
    if (!threaded_buffer.empty())
      mooseError(_name, ": Threaded buffer is not empty in preExecute()");
  if (_receive_buffer->currentlyReceiving())
    mooseError(_name, ": Receive buffer is not empty in preExecute()");
  for (const auto & map_pair : _send_buffers)
    if (map_pair.second->currentlySending() || map_pair.second->currentlyBuffered())
      mooseError(_name, ": Send buffer is not empty in preExecute()");

  // Clear communication buffers
  for (auto & send_buffer_pair : _send_buffers)
    send_buffer_pair.second->clear();
  _send_buffers.clear();
  _receive_buffer->clear();

  // Clear counters
  _local_chunks_executed = 0;
  _local_objects_executed = 0;
  _local_objects_started = 0;
  _total_objects_started = 0;
  _local_objects_finished = 0;
  _total_objects_finished = 0;

  _currently_pre_executing = true;
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::execute()
{
  if (!_currently_pre_executing)
    mooseError(_name, ": preExecute() was not called before execute()");

  _currently_pre_executing = false;
  _currently_executing = true;

  switch (_method)
  {
    case ParallelStudyMethod::SMART:
      smartExecute();
      break;
    case ParallelStudyMethod::HARM:
      harmExecute();
      break;
    case ParallelStudyMethod::BS:
      bsExecute();
      break;
    default:
      mooseError("Unknown ParallelStudyMethod");
  }

  _currently_executing = false;

  // Sanity checks on if we're really done
  comm().barrier();

  if (!_working_object_buffer->empty())
    mooseError(_name, ": Working object buffer is not empty after execution");
  mooseAssert(std::all_of(_working_object_buffer->data().begin(),
                          _working_object_buffer->data().end(),
                          [](const auto & object) { return !object; }),
              "Raw working buffer has objects after execution");
  for (const auto & threaded_buffer : _threaded_during_execution_objects)
    if (!threaded_buffer.empty())
      mooseError(_name, ": Threaded object buffer is not empty after execution");
  if (_receive_buffer->currentlyReceiving())
    mooseError(_name, ": Receive buffer is not empty after execution");
  for (const auto & map_pair : _send_buffers)
    if (map_pair.second->currentlySending() || map_pair.second->currentlyBuffered())
      mooseError(_name, ": Send buffer is not empty after execution");
}

template <typename Object, typename Context>
std::string
ParallelStudy<Object, Context>::addObjectError(const std::shared_ptr<Object> & /* object */,
                                               const AddObjectError error) const
{
  if (error == AddObjectError::ADDITION_DURING_EXECUTION_DISABLED)
    return _name + ": The addition of objects during execution is disabled.\n\nThis feature can be "
                   "enabled with the '_allow_addition_during_execution' private parameter";
  if (error == AddObjectError::PRE_EXECUTION_AND_EXECUTION_ONLY)
    return _name + ": Can only add objects in the pre-execution and execution phase (between "
                   "preExecute() and the end of execute())";
  if (error == AddObjectError::PRE_EXECUTION_ONLY)
    return _name + ": Can only add objects in the pre-execution phase (between preExecute() and "
                   "execute())";
  if (error == AddObjectError::PRE_EXECUTION_THREAD_0_ONLY)
    return _name + ": Can only add objects in the pre-execution phase (between preExecute() and "
                   "execute()) on tid 0 (not thread safe)";

  mooseError("Unknown AddObjectError");
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::moveObjectToBuffer(std::shared_ptr<Object> & object,
                                                   const THREAD_ID tid)
{
  if (_currently_executing)
  {
    if (!_allow_addition_during_execution)
      mooseError(addObjectError(nullptr, AddObjectError::ADDITION_DURING_EXECUTION_DISABLED));
  }
  else if (!_currently_pre_executing)
  {
    if (_allow_addition_during_execution)
      mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_AND_EXECUTION_ONLY));
    else
      mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_ONLY));
  }
  else if (tid != 0)
    mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_THREAD_0_ONLY));

  if (object)
  {
    // Objects added during pre-execution go directly into the work buffer (not thread safe)
    if (_currently_pre_executing)
    {
      ++_local_objects_started; // must ALWAYS increment when adding work to the working buffer
      _working_object_buffer->move(object);
    }
    // Objects added during execution go into a temprorary threaded vector (is thread safe) to be
    // moved into the working buffer when possible
    else
      _threaded_during_execution_objects[tid].emplace_back(std::move(object));

    mooseAssert(!object, "Object was not moved");
  }
}

template <typename Object, typename Context>
void
ParallelStudy<Object, Context>::moveObjectsToBuffer(std::vector<std::shared_ptr<Object>> & objects,
                                                    const THREAD_ID tid)
{
  if (_currently_executing)
  {
    if (!_allow_addition_during_execution)
      mooseError(addObjectError(nullptr, AddObjectError::ADDITION_DURING_EXECUTION_DISABLED));
  }
  else if (!_currently_pre_executing)
  {
    if (_allow_addition_during_execution)
      mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_AND_EXECUTION_ONLY));
    else
      mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_ONLY));
  }
  else if (tid != 0)
    mooseError(addObjectError(nullptr, AddObjectError::PRE_EXECUTION_THREAD_0_ONLY));

  // Get object size before hand so we can resize
  std::size_t num_objects = 0;
  for (const std::shared_ptr<Object> & object : objects)
    if (object)
      ++num_objects;

  // See comments in moveObjectToBuffer() as to why this distinction exists
  if (_currently_pre_executing)
  {
    if (_working_object_buffer->capacity() < _working_object_buffer->size() + num_objects)
      _working_object_buffer->setCapacity(_working_object_buffer->size() + num_objects);
    _local_objects_started += num_objects;
  }
  else
    _threaded_during_execution_objects[tid].reserve(_threaded_during_execution_objects[tid].size() +
                                                    num_objects);

  // Move the objects
  for (std::shared_ptr<Object> & object : objects)
    if (object)
    {
      if (_currently_pre_executing)
        _working_object_buffer->move(object);
      else
        _threaded_during_execution_objects[tid].emplace_back(std::move(object));

      mooseAssert(!object, "Object was not moved");
    }
}

template <typename Object, typename Context>
unsigned long long int
ParallelStudy<Object, Context>::sendBufferPoolCreated() const
{
  unsigned long long int total = 0;

  for (const auto & buffer : _send_buffers)
    total += buffer.second->bufferPoolCreated();

  return total;
}

template <typename Object, typename Context>
unsigned long long int
ParallelStudy<Object, Context>::objectsSent() const
{
  unsigned long long int total_sent = 0;

  for (const auto & buffer : _send_buffers)
    total_sent += buffer.second->objectsSent();

  return total_sent;
}

template <typename Object, typename Context>
unsigned long long int
ParallelStudy<Object, Context>::buffersSent() const
{
  unsigned long long int total_sent = 0;

  for (const auto & buffer : _send_buffers)
    total_sent += buffer.second->buffersSent();

  return total_sent;
}

template <typename Object, typename Context>
unsigned long long int
ParallelStudy<Object, Context>::poolObjectsCreated() const
{
  unsigned long long int num_created = 0;

  for (const auto & pool : _threaded_object_pools)
    num_created += pool.num_created();

  return num_created;
}
