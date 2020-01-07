
Finite Volumes Notes
=====================

Assembly
---------

This is info about how DG kernels get their values computed and stored into assembly:

* Assembly holds vectors that collect residual values.

* TaggingInterface has local member variables that collect residual values for
  specific tags.  Subclasses store info in _local_re and then call its
  accumulateTaggedLocalResidual function to copy accumulated local vec
  residual values into the assembly residual vector that is accessed by
  pointer/reference.  This usually is called in e.g. Kernel::computeResidual,
  DGKernel::computeElemNeighResidual, etc.

* The compute residual loop has an onInternalSide function that is called on
  each elements' sides.  This is for DGKernels.  The DGKernel's regular
  compute residual functions only get the computed residual stored in assembly
  for the current element of the current internal side.  So the onInternalSide
  loop manually has to get this residual contribution added to the neighbor
  element too by calling FEProblem::addresidualNeighbor which takes the
  current accumulated residual values in Assembly and adds them to the
  neighbor residual vector.

This separation of residual collection is silly.  Element and neighbor
residual accumulation into assembly should not be performed in two separate
ways in very different places.  

Face residual calculations for FV will always be terms that underwent the
divergence theorem conversion to a surface integral and so will always be
dotted with the unit normal vector.  So we can skip recomputing the residual
contribution to each neighboring element of a face by just computing it for
one element and then inverting the residual value's sign for the other
element.


