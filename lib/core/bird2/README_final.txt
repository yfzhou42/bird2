The final project is based on the bird2 code base where the collision simulation is inherited from our bird2 
implementation (which can be buggy).

The compile and run command are the same as in bird2. 

In UI:
To simulate fracture using Force Absorption approach, keep "Fracture Enabled" on and turn off "Lagrangian 
Break ENabled". 
The brittleness of the rigid body is simultaed using "Fracture Max Strain", which needs to be specified 
before the object is loaded.

To simulate fracture using Lagrangian Multiplier approach, keep "Fracture Enabled" on and turn on "Lagrangian 
Break Enabled".
The brittleness of the rigid body is simultaed using "Lagrangian Max".

The "Max Fracture Iteration" simulates how many generations a body can possibly shatter, so that we won't have
infinitely small subdivisions.
The "Max Strain Multiplier" expoenntially increases the hardness of breaking a rigid body after each generation
of shattering. 