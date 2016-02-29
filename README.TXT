To run the Maureparticle example:

1) Compile the Maureparticle and BUILT_TABLES programs. e.g.
     
    gfortran -o maurept.exe Maureparticle_3-27-2015.f
    gfortran -o buildtables.exe BUILD_TABLES.f

2) Move the executables and change directory into the example
   directory.

2) Run the example ADCIRC simulation.
   This is a short run on the ec_95d grid with 
   tides and a constant ENE wind. The FORT.14,FORT.13, 
   and FORT.15 files are in the example directory. You'll
   need a recent ADCIRC executable to do this. I tested with
   the v49_release, but anything later than that should work 
   too. 

3) Run the BUILD_TABLES program to genereate the 
   NODE2EL.TBL and EL2EL.TBL files required by Maureparticle.
   you can do this before ADCIRC has been run since it only 
   needs the FORT.14 file. the EL2EL table can be slow to
   generate for large grids (a faster version is in the works). 

4) Now run Maureparticle. It generates a file called 
   MAUREPT.OUT that contains the particle ID, X,Y location, time,
   and the element the particle is in. 

5) run the ParticlesOut_2_kmz.pl Perl script to genereate a
   kmz file to view the results in Google Earth.

6) Check out the PARTICLES.INP input file. Try some different 
   parameters. It should be self-explanatory. 




