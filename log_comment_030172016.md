#MIMMO TASKS

##1. MimmoNamespace e BaseManipulation / ROCCO

 - 1.1. 5 todos in basemanipulation

 - 1.2. Comments 

 - 1.3. line 124 and 149 static_cast unique pointer? if no virtuals are present in the base class, dynamic_cast on reference fails. The only way to do this job, 
   in template context, is using static_cast (checked at compiler runtime) and hoping theuser knows how is doing. But, static_cast works. [&#10004;]
		 

 - 1.4. add mimmonamespace to each mimmo class
 
 - 1.5  m_ndeg/m_displ no more required in base class. 


##2. Apply / ROCCO

 - 2.1. Update execut with new input/output pins  [&#10004;]
 - 2.2. Comments

##3. GenericInput / ROCCO

 - 3.1. Better implementation
 - 3.2. when reading files is not exactly a generic input but is specialized only to dvecarr3E data  
 - 3.3. Comments

##4. Info

 - 4.1. recoding manipulation objects to avoid the use of info (?) DOES NOT EXIST ANYMORE [&#10004;]


##5. FFD Lattice / EDOARDO

- 5.1. modify execute in order to use m_result of base class

- 5.2. use input of base class and avoid the use of displ

- 5.3 	Atomize set/get methods for paramaters & output of classes, according to new pin inout style

- 5.4  Adjusting set/get of primitive class BasicMesh and included Class BasicShape

- 5.5. Comments


##6. Bend / ROCCO

- 6.1. modify execute in order to use input/output of base class (deprecated use of displ and ndeg)  [&#10004;]

- 6.2. useInfo maybe useless

- 6.3. Comments

##7. Mask / ROCCO

- 7.1. modify execute in order to use input/output of base class (deprecated use of displ)  [&#10004;]

- 7.2. useInfo maybe useless

- 7.3. Comments

##8. OutputDoF / ROCCO

- 8.1. modify execute in order to use input/output of base class

- 8.2. modify execute to avoid parent/child relationship (deprecated)

- 8.3. useInfo maybe useless

- 8.4  must become a GenericOutput class

- 8.5. Comments

##9. Rotation and Translation Box / ROCCO

- 9.1. modify execute in order to use input/output of base class or new parameter values for rotation and translation

- 9.2. modify execute to save the result in result and then use input/output pins

- 9.3. useInfo maybe useless

- 9.4. Comments

##10. InOut Pin block/ EDOARDO

- 10.1 Comments


##11. Mimmo Tests/ ROCCO/ EDOARDO

- 10.1. make them coherent with new basemanipulation objects design

- 10.2. Comments

##12. Chain ROCCO

- 11.1	managing total DOF of your manipulation  

- 11.2. Comments
