An alpha version of the rearchitected FLASH multiphysics, multiscale simulation code - __FLASH5__ - is now available in this Github repository. Many of the physics capabilities and some of the infrastructural capabilities of the current full official release of FLASH4 are not included in this version.  Major new changes in FLASH5 are:

* The primary AMR package is [AMReX](https://amrex-codes.github.io) which is configured to provide an octree-like mesh.  Support for Paramesh still exists. Uniform Grid is also supported, though some of its features are deprecated.
* There are several fundamental changes to the code architecture.
   * Looping over blocks is replaced by smart iterators.
   * Tiling (within blocks) is supported.
   * It is possible to iterate over block/tiles one level at a time or all at once.
   * Many interfaces for querying the Grid unit are deprecated in favor of making the corresponding information available through the iterators and block/tile objects.
   * Knowledge about blocks and mesh is further abstracted away from physics units. There is no entity called "blockID" unlike in FLASH4. Any paramesh functionality depending upon the knowledge of blockID is now internal to the paramesh implementation of the Grid unit.

Some current known limitations are:
* Restarting from checkpoints doesn't work.
* Non-trival tiling only works with AMReX. 
* Threading with OpenMPI does not work with Paramesh.
* Threading with OpenMPI is not properly tested.
