<img width="1138" alt="neuron_img" src="https://user-images.githubusercontent.com/46076747/209508271-8afe5961-a7e7-4a80-afa3-ecbe90d60003.png">
# trees-prep

## rwsim

- Modified implementation of lookup table for realistic voxelized geometries of neuron skeletons
- Conversion of swc to voxelized Geometries

## In Progress

- Simulation Code

Using a voxelized geometry has its limitations when modeling the complex structures of neurons. As the quality of a neuron skeleton increases the voxel size necessary to accurately represent the skeleton decreases. This leads to overwhelming volume sizes that make simulations a computational nightmare. This repo uses a modified lookup table to precomputed cellular components in any given voxel. This eliminates runtime computations during Monte Carlo simulations.

A 3d sparse matrix indicates whether any elements are in a given voxel. For voxels containing elelemts, a linear index of a 1d cell array is stored and will return a list of elements within a that voxel. Nesting the sparse matrix and linear indexed lookup table provides a very effiecient method to determine what elements are in a given voxel.

