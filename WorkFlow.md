# WorkFlow

## Section 1

The purpose of section 1 is to read an swc file, translate the data to a voxelized geometry and then generate a lookup table. The swc is preprocessed using the trees toolbox. Preprocessing as of now involves resampling the tree for consistent connection distances and scaling the swc file to what I think is [um].

### PrepSim

This is the main function of section 1.

Args:

- voxelscale:
voxelscale will scale the voxel size determined by the minimum distance between nodes in our graph.

- filename:
filename is the name and path (global path or root path) to the target swc file.

- **scale:
scale is one of the arguments I need help with. It's intented use is to scale the swc coordinates from whatever they are in neuroglancer to micrometers [um].**

The first step we take is reading the swc file and scaling it. This is a redundancy given our preprocessing in trees prep can scale the data. We pass 1 to account for this.

Next we calculate the distance of each connection except the root node using dists = setdists(tree).

#### initbounds

[b, swc, boundSize, pairs, vsize, ranges] = initbounds(tree, dists, voxelscale);

Args:

- tree

- dists

- voxelscale

initbounds creates a copy of tree called swc to modify its data structure and values.

The first change is setting the parent of the first node to itself. This allows us to reference its parent without and index error and saves a lot of time writing clauses in our code for if (id ~= 1).

We then use hist counts to threshold the Radii.
Threshold is set at the average of the first bin.

we isolate the [x,y,z] values and determine the range as min and max [xyz+||-r]. *r is the thresholded radii.*

Next we set the vsize variable as min(r)*voxelscale.
We translate the xyz values using minimum values in ranges and scale the coords using vsize.

***Note***: a padding of 10 voxels is added to the coords after this change.

Radius is then scaled with vsize and ranges variable is reset to max,min of voxelized [xyz+||-r].

the bounds of our geometry is defined as the max values of ranges

***Note***: a padding of 10 voxels is added to the bounds.

XYZR are rewritten into swc.

The bounds of each pair is then determined from calcBounds(swc) using the same process of [xyz(i)+||-r(i)]

Output:

- b: bounds of each pair.

- swc: scaled&translated swc.

- boundSize: bounds of geometry.

- pairs: pairs of swc.

- vsize: voxel scalar.

- ranges: range of voxelized geometry.

The next step in prepSim is to generate the Lookup Table.

Args:

- boundSize: bounds of geometry.

- b: bounds of each pair.

The first step in generating our lookup table is to create a cell array to store node ids and initialize our LUT as zeros(boundSize).

Next we need to enumerate all voxels within the range of each pair.

We use cellfun to call extract_Range for each pair.

#### extract_Range

Args:

- sub: range of pair i

I believe I've reinvented mesh grid or a similar function but this essentially takes the range of values of one dimension and creates copies for the other dimensions, thus enumerating all combinations of [x,y,z];

Output:

- ci: all subscript indicies of voxels within range

Accessing values using subscripts is very very slow in matlab when the matrix is large.

We convert each subscript to linear indices using cf_sub2ind.

Args:

- bsize: bounds of geometry
- ci: subscripted indicies for pair i.

Output:

- inds: linear indicies of all voxels pertaining to pair i

Now that we have the linear indicies of all voxels each pair spans, we just have to iterate over each pair and add its value to the lookup table.

We don't want to overwrite values so instead of storing the pair index in our lookup table we store a reference index to our cell array.

If the value is zero in our lookup table then no pairs have overlapped with this voxel. We store the next index of our cell array and store the pair id at that index in our cell array

If the value is not zero then a pair has already overlapped with this voxel. We append our pair id to the list stored in our cell array at the index from our lut.

Output:

- A: cell array containing pairs for given voxel.

- LUT: Lookup Table containing indices of cell array.

The last step in prepSim is to remove empty cells from A. Growing the cell array by 1 each time we find a voxel that hasn't been visited is very inefficient. Instead I have it growing by a large amount so we only need to resize a few times over the course of our lookup table generation.

Output:

- LUT: Lookup table containing indicies of cell array.

- B: cell array giving node ids not pair ids

- pairs: pairs of swc

- boundSize: bounds of geometry

- swc: translated&scaled swc

- memoized_distance: distance of each pair

- A: cell array of pair ids for given voxel

- vsize: voxel size

- range: range of our voxelized geometry
