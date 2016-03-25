#MIMMO TASKS

##1. MimmoNamespace e BaseManipulation

 - 1.1. 6 todos in basemanipulation  [&#10004;]
 - 1.2. Comments 
 - 1.3. line 124 and 149 static_cast unique pointer? If no virtuals are present in the base class, dynamic cast on reference fails. The only way to do this job, in template context, is using static cast (checked at compiler runtime) and hoping theuser knows how is doing. But, static cast works. [&#10004;]
 - 1.4. add mimmonamespace to each mimmo class [&#10004;]
 - 1.5  members ndeg/displ no more required in base class. [&#10004;]
 - 1.6  Removal of internal pins only of a BaseManipulation class, without external removalPin method-> is meaningful? When a class is destroyed, pins with other class survive eventually. how to manage this?  remove all input/output pins. Party ended. [&#10004;]
 - 1.7  Copy Operator and Constructor Done managing m_input, m_result (IOdata modified in copyOp and CopyConstr) [&#10004;]

##2. Apply

 - 2.1. Update execut with new input/output pins  [&#10004;]
 - 2.2. Comments

##3. GenericInput

 - 3.1. Template read implementation [&#10004;]
 - 3.1. Template set input implementation [&#10004;]
 - 3.3. Comments

##4. Info [&#10004;]

 - 4.1. recoding manipulation objects to avoid the use of info (?) DOES NOT EXIST ANYMORE [&#10004;]


##5. FFD Lattice

- 5.1. modify execute in order to use m_result of base class  [&#10004;]
- 5.2. use input of base class and avoid the use of displ (m_displ for the moment set as member of lattice, input of base class is temporary)  [&#10004;] 
- 5.3 	Atomize set/get methods for paramaters & output of classes, according to new pin inout style  [&#10004;]
- 5.4  Adjusting set/get of primitive class BasicMesh and included Class BasicShape  [&#10004;]
- 5.6. Comments 


##6. Bend

- 6.1. modify execute to use input/output of base class (no displ and ndeg) [&#10004;]
- 6.2. useInfo maybe useless [&#10004;]
- 6.3. Comments

##7. Mask

- 7.1. modify execute in order to use input/output of base class (deprecated use of displ)  [&#10004;]
- 7.2. useInfo maybe useless [&#10004;]
- 7.3. Comments

##8. OutputDoF

- 8.1. modify execute in order to use input/output of base class [&#10004;]
- 8.2. modify execute to avoid parent/child relationship (deprecated) [&#10004;]
- 8.3. useInfo maybe useless [&#10004;]
- 8.4  must become a GenericOutput class [&#10004;]
- 8.5  template data to ofstream maybe doesn't work, try to use it [&#10004;]
- 8.5. Comments

##9. Rotation and Translation Box 

- 9.1. modify execute in order to use input/output of base class or new parameter values for rotation and translation [&#10004;]
- 9.2. modify execute to save the result in result and then use input/output pins [&#10004;]
- 9.3. useInfo maybe useless [&#10004;]
- 9.4. Comments

##10. InOut Pin block  [&#10004;]

- 10.1 Comments  [&#10004;]
- 10.1 Template .tpp file  [&#10004;]

##11. Mimmo Tests

- 10.1. make them coherent with new basemanipulation objects design  [&#10004;]
- 10.2. Comments 

##12. Chain

- 11.1	managing total DOF of your manipulation  
- 11.2	Useful to bypass some blocks in the chain: Switch mechanism or break/unbreak pins? to be decided and scheduled
- 11.3  Loop detection and blocking  [&#10004;]
- 11.4  Tracking of DOF along the chain (11.1?)
- 11.5 	Comments

##12. Lattice - new object

- 12.1 Make the FFDLattice derived from here ?
- 12.2 Comments

##13. SHAPE

- 13.2 Check setRefSystem with one axis as input (not orthogonal result?)

##14. RBF MORPHERS

- 14.1 Add more geometric interfaces
- 14.2 Comments

##15. MISC

- 15.1 Object that transform coordinates/DOF needed (?).

