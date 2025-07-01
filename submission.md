## FEM: A Finite Element Manipulation-Based Method for Physics Simulations

## Design Choices
- In this project, I implemented an FEM-based simulation to animate deformable solid objects using tetrahedral meshes, computing internal elastic and damping forces, apply gravitational forces, and resolving collisions that occur as a result of these forces when moving the simulation forward in time.

Below is a descrption of how each of the follow features were implemented
- Extracting the Surface Mesh: Explicitly, I extracted the surface mesh by iterating through every face on every tetrahedral mesh, and identifying faces that belonged only to one tetrahedron, as this indicates that this faces lives on the exterior of an object, on its surface. Any internal tetrahedral face that lives inside a mesh will belong to two tetrahedrons. I then used the ordering of the vertices used in the hard coded single tet that we were provided to pass the vertices and faces in the correct order to opengl. I also used that ordering to compute face normals (which were later used when applying the stress tensor).

- Computing and applying internal forces: Internal elastic force was computed using the deformation gradient F, which is obtained from matrix "P" and "B" (in the O'Brien and Hodgins paper) which are the deformed and undeformed tetrahedron vertex positions. Green's strain tensor, F*F-1 - I is then used to determine elastic stress using the sigma_E equation from the Lecture slides), and viscuous stress is applied by computing the velocity gradient within each tetrahedral element and using it to derive the damping stress tensor. The elastic and viscuous stress are then added, and then converted into per-node forces via area-weighted surface normals.

- Collision Resultion: Collisions between the ground and a deformable object are resolved by checking the y-component of each vertex's position to see if it is at the same y-component of the ground height. If so, then a collision is detected and the deformable object's velocity is decomposed to a normal and tangential to ground component. The normal component is scaled by the restitution coefficient and the tangential component is scaled by friction, and the scaled velocity components are summed back together to create the updated velocity for a given vertex. 

- Integration method: The midpoint integration method was implemented to move the simulation forward in time. Forces are computed at an intermaediate half-step, and then each particle's acceleration is derived from its mass and applied forces and then used to update the velociy and position. Forces are applied again, and this second half-step we resolve collisions now that we have a better estimate of acceleration, velocity, and position having looked ahead a time-step.

- For extra credit I implemented a method called export_obj(). What this does is write an obj file for the mesh at every frame. I then brought the files into blender as an obj sequence, where I rendered it. I also applied subdivision surface to the sphere to see what happened and it ended up looking great. You can see it in my video.

    
## Links to the starting lines of the implemenetation of the following features
- [Extract the surface mesh from your tetrahedral mesh](src/simulation.cpp#L26)
    - [Compute and apply force due to gravity](src/simulation.cpp#L201)
    - [Compute and apply internal elastic forces](src/simulation.cpp#L225)
    - [Compute and apply internal viscous damping forces](src/simulation.cpp#L225)
    - [Resolve collisions](src/simulation.cpp#L212)
    - [Explicit midpoint method](src/simulation.cpp#L264)
    - [Any extra features](src/simulation.cpp#L1)

## Results

![Alt Text](blender_videos/single-tet0001-0240.mp4)
![Alt Text](blender_videos/sphere0001-0700.mp4)

