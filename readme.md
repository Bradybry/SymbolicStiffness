This was my simple solution to checking my homework and engaging with the material for my structural analysis class. I started with a lazy implementation that leverages SymPy to solve the system rather than partitioning and solving with just NumPy functions. In the process, I discovered that my lazy implementation enables structures to be solved symbolically, which is a huge benefit because several of my HWs have required that and I have been left with difficulty checking assignments.

The code is hacky and undocumented. Slightly OO so yay. My philosophy was to hack something together that would allow me to ensure that I am doing my HW right every step of the way. For that reason, the program calculates things that you wouldn't necessarily need or want if you were just worried about performance. But that means that you can access all of the transformation matrices, local and global stiffnesses for each element, etc. 

I didn't not implement any good output functions, so you will have to manually grab the final displacements, reactions, and internal forces. I think this is a slight advantage because most other implementations leave you a singularly formatted answer that must be parsed, and only then can you leverage the results accordingly in a different environment. Feel free to make you own output methods and visualizations, I would certainly appreciate it. 


This has not been extensively tested and probably is really broken. It definitely does not scale well but should be able to handle simple symbolic structures and relatively large 2-d defined structures. 

The program can only handle nodal forces atm, I plan to add in distributed forces along the members once that has been covered in class. All of the class methods are weakly implemented. For example, apply loads require that you overwrite all of the forces at the node even if some of them remain unknown. I would suggest bypassing the methods in such cases and manually editing the attributes of the nodes as needed. 

I am publishing this because I want to give this tool to anyone in a similar position and I would also love input from others on what I am doing really wrong etc. I will definitely refactor and reimplement several things throughout the semester as we learn more about what is possible in the class. For example, I plan to reimplement Structure as an abstract class and create two child classes, one that solves things quickly via partitioning and is not dependent on sympy and another that employs the current method. Moving towards abstract classes will probably be an advantage so that when we move onto 3d structures it will be a less complicated transition and I can maintain 2d frame support. 

I will not take any criticism too harshly, thanks.

Check out the examples Ipython notebook to get up and going

Dependencies are numpy and sympy and Python 3.6.
