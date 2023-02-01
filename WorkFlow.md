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


## Section 2

We will take our simulation parameters and create a simulation object.

We do this by initializing a random_walker_sim class

### random_walker_sim  

Args:

- LUT: Our lookup table.
- B: node ids {childList, parentList} for voxels.
- pairs: pairs of elements in swc.
- boundSize: bounds of geometry.
- swcmat: swc in matrix format.
- step_size: step size of sim.
- perm_prob: permeation probability of sim.
- step_num: number of steps.
- init_in: Boolean repr if particles must start inside the cell
- particle_num: number of particles to initialize
- memoized_distance: distance between connections.

We first set our arguments to our object props.
The only additional step in the constructor is to initialize our particles.

#### init_particles(iter,chunksize)

Args:

- obj: simulation obj.
- iter: number of iterations in our simulation.
- chunksize: depreciated argument.

We create a cell array of size particle_num x 1.

We the loop over each element of that cell array and create a randomwalker object with a random initial position.

Output:

- particles: cell array of initialized particles.

After initializing particles we return the simulation object

Output:

- sim
  - PARAMS:
    - lookup_table: Lookup table of sim.

    - index_array: index array of sim.

    - pairs: pairs of swc.

    - swc: swcfile in matrix format.

    - step_size: step size of sim.

    - perm_prob: permeation probability.

    - boundSize: bounds of geometry.

    - rwpath: empty matrix that will store positions of particles.

    - logical_alloc: depreciated variable.

    - currstate: array of booleans repr states of particles (inside||outside).

    - particles: cell array of particle objects;

    - particle_num: number of particles.

    - memoized_distance: distances between each pair

    - chunkSize: depreciated.

    - path: path to simulation output.

## Section 3

After we have created a simulation object we are ready to run the simulation.

### eventloop

Args: 

- obj: simulation object.

- iter: number of steps the simulation will perform for each particle (redundant)

eventloop is the main function of section 3.

It starts by initializing the output files and paths.

### obj.setsimulationpath

Args:

- iter: number of steps in sim.

- obj.particle_num: number of particles in sim.

Create a new directory for this simulations results.

With larger step_nums and particle_nums saving the results to a variable can exceed the memory of local machine.  Instead we use a mat file to write to specific elements of file without reading full file into workspace.

We initialize a matfile to store all particle positions per step. SIZE:
(iter, particle_num, 3)

We also initialize a matfile that stores our initial particle positions. SIZE: (particle_num,3)

Output:

- path: path to simulation results 

- fileName: fileName of matObj.

- matObj: handle for matfile of all particle positions per step.

- initPosObj: handle for matfile of particle inital positions.

As mentioned before the simulation results can't be stored entirely in memory. We define chunk_iter as the number of particles we iterate over before writing to matfile.

The next step of event loop is to initialize our parallel pool. We have several constants that can be static on our workers and therefore only require initialization once.

We use variable PARVARS to store these constants.

finally we can iterate over the particles in a batch.

We extract the data of our batch from the matobj and set it to zero.

We then extract the relevant particles in our batch.

We use parfor to iterate over each particle in batch.

Create variables for the particles path and extract constants from worker variable.

We initialize the random unit vectors for each step.

iterate over each step accessing the random direction and perform the computations for the step.

### cellgap2

Args:

- particle: [x y z state flag] of particle.

- swc: swc matrix.

- index_array: index array of sim.

- boundSize: bounds of geometry.

- lookup_table: lookup table of sim.

- step: step its on.

- perm_prob: permeation probability.

- crand: current random vector.


We start by defining position,state,and flag from our particle argument.

The function enters a control sequence.

If the flag is true we just collided with a boundary and didn't permeate. We should stand still for 1 step.

Otherwise we find the next position of the particle

#### setnext

Args:

- position: current position of particle

- step: step size.

- crand: random unit vector.

We take the product of our vector and our step and add it to the current position of the particle.

Output:

- next: next position of particle.

After determining the next position of the particle we must decide whether the particle can make this step.

#### checkpos

Args:

- next: next position of particle.

- swc: swc matrix

- index_array: index array of particles.

- state: state of current position

- boundSize: bounds of geometry

- lookup_table: lookup table of simulation.

the first step in check pos is to extract the [child, parent] ids of any segments near the location of the particle. 

#### preprocesses

Args:

- next: next position.

- boundSize: bounds of geometry.

- lookup_table: lookup table of sim.

We convert the coordinate location to voxel of our next position.

Generate a logical array of elements within range 0 < n < boundSize.

Extract the voxels that exist within the range of our lookup table.

These voxels are in subscript indicies of our lookup table. Subscript indices are slow in matlab so we convert them to linear indices.

We return the nonzero elements of our lookup table at these linear indicies.

Output:

- Indicies: Indicies obtained from our lookup table that reference data in our index array.

There are two outcomes from preprocesses:

Either there are elements of the cell near the particle or there are not.

if there are we check each connection to determine whether a boundary is crossed by executing this step.

if the are no elements near the particle the no boundaries were crossed but the particle is not inside the cell.

CASE 1: There were elements near particle.

#### check_connections

Args:

- indicies: indicies of our index array.

- A: index array.

- swc: swc matrix.

- next: next position.

- currstate: current state of particle.

Extract child and parent ids from index array.

Extract child parent coords from index array.

Init variables x,y,z of next.

Iterate over node pairs checking connections using pointbasedswc2v.

set nextinside to result of iteration or its previous value.

if the particle is currently inside a connection and we find a connection that the particle will be inside if the step is taken we can break our loop and return.

Output:

- currinside: current state

- nextinside: next state

We use insideLogic to enumerate the boolean logic of the particle state.

#### insideLogic

Args:

- currinside: current state

- nextinside: next state

Simple Control Sequence enumerating the possible states.

if both are inside: no boundaries crossed.

if the current position is inside and the next isn't: boundary was crossed: OUT.

if the current position isn't inside and the next is: boundary was crossed: IN.

if the both positions are outside then no boundary was crossed.

Output:

- inside: boolean for determining whether to freely step or apply and control sequence on step.

***This can be made much better we can create an enum or class for the type of step it is***

The determination of whether to control the step or freely step is returned to cellgap2.

Output:

- Collided: boolean for control step or free step.

if there wasn't any collision, we freely step and set our current position to next.

If there was a collision we determine whether the particle permeates using rand and our permeation probability.

if it permeates we update the position and update the state.

if it doesn't we do not update position and state. we set our flag to true.

Lastly we rewrite the position, state, and flag to our particle and return it.

Output:

- Particle: [x,y,z,state,flag];

The last component of a step is to store the coordinate positions of the particle for that step.

After iteration over each step we write the particle positions for that particles simulation to our temporary data variable.

After the parfor has iterated over each particle in that chunk we write the position data to our mat file.

Once each chunk has been completed we are finished.
