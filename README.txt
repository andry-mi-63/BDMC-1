The code is Bold Diagrammatic Monte Carlo algorithm
for Frohlich and Holstein model. 

Description of the algrithm: 

A. S. Mishchenko, N. Nagaosa, and N. Prokof'ev: Diagrammatic Monte Carlo method for manypolaron problems, Phys. Rev. Lett. 113, 166402 (2014).
  

1 The code is compiled with batch "ifco.bat" using intel Fortran:

.\ifco.bat FHO_25.exe FHO_25.f90

2. Input files:

a) "in.in": model parameters, explained in file  

b) "probab.in": updates probabilities, explained in file

c) "flat.in": flat probabilities of Fenmann diagrams order

d) "gre_gde.in" momenta where imaginary time Green function
                are calculated 