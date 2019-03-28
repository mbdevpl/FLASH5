#if 0
  Use these macros to set the value associated with each point to write to UNK.
  For a given point, set the point's value to REFINE_TO_L3 if it intended that
  AMReX refine the block containing the point up to level 3.

  NOTE: These mass values were chosen in a non-physical way for programming
  convenience.  In particular, not that at the finest level of refinement,
  mass = density, which is not physically correct for the domain specified in
  flash.par.
#endif

#define REFINE_TO_L1 1
#define REFINE_TO_L2 2
#define REFINE_TO_L3 3
#define REFINE_TO_L4 4
#define REFINE_TO_L5 5

