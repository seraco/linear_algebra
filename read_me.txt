1. upper_triangularisation.m is the main file.
2. In order to run this in your machine, copy all these files in a folder and point your MATLAB directory to this folder and then do the following:
	(i) Define your matrix.
	(ii) Call your upper_triangulisation subroutine by passing your matrix and its size. 
	For example:
	>> clear
	>> A= [1,4,7; 2,5,8; 3,6,10]

	A =
	
	     1     4     7
	     2     5     8
	     3     6    10
	
	>> U = upper_triangulisation(A,3)
	
	U =
	
	     1     4     7
	     0    -3    -6
	     0     0     1
	
	>> 
3. The main function (i.e. upper_triangulisation.m) calls five subroutines to perform the assigned duties. These subroutines are in-line with the theory explained in the class, and the details of these subroutines are clearly written in their respective files.

4. The instructor strongly believes the students would now be able to convert any theory into MATLAB codes after going through this example.

5. Finally, you will be asked to build further on these functions for your coursework and the instructor expects the new scripts meet the standards maintained in the scripts uploaded here (i.e. a neat presentation and proper comment as what you are doing in a particular line of the script). Best wishes!
