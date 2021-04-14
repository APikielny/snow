To run:
I didn't modify the command-line arguments.  
You can edit the lambda, mu, phi and psi parameters as well as gravity, the acceleration constant and the mesh file to load in simulation.h

Implementation:
Surface Mesh Extraction:
To extract the surface mesh from the tetrahedral mesh I store the faces in a vector and then create a second vector of faces that I populate based on whether or not a face "pair" exists. I decided that a face "pair" is a face that shares all three of its vertices with another face which means that both faces are internal. I check to see if the face that I'm adding to my face list shares vertices with a face that's already in the list - if it does, I remove the existing face from the list and keep looping. If it doesn't, I add the face to the list. If you try to load the ellipsoid or cube meshes, however, it errors when trying to remove a face pair from the list - I'm not sure why this happens and was not able to debug it. I was able to render those shapes by just commenting out line 123 in simulation.cpp. 

Stress/Strain Computation:
To help with the storage of values, I created tet structs with fields to store face areas and normals, particle indices and beta matrices. I created particle structs with fields to store position, velocity, midpointPosition, midpointVelocity and force. I implemented the midpoint integrator and computed stress and strain using the formulas from class. 

Results: 
I recorded three videos in the results folder. The "ellipsoid with sphere+ground" and "not exploding sphere" have the parameters as set in the README except psi and phi are increased by a factor of 2. These demonstrate collision, gravity, and stress/strain forces. The "exploding sphere" video shows what happens when I reduce the phi and psi values to 10 - it bounces but quickly becomes unstable. 

