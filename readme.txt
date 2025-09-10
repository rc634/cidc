---------------------------------------------
About the Choptuik Initial Data Code (CIDC)
---------------------------------------------


         -*- -*- -*-

-------------------------------
- ABOUT THE CODE 
-------------------------------

main.cpp -- the main thign we run
- outputs data as csv file 

RK4solver.cpp/hpp - performs rk4 integration

scalar_profile.cpp/hpp - provides scalar field, u and du/dr

rhs.cpp/hpp - provides the Right Hand Side (RHS) of differential equaitons 
- for example psi' = rhs.psi(...)

plot.py - python matplotlib plotting software to visualise the solution

         -*- -*- -*-

-------------------------------
- CHANGE CODE PARAMETERS
-------------------------------
Iniside main.cpp you can modify some parameters:

r1 -- the final radius

dr -- is the stepsize for the initial data

-- total steps is aprox (r1-r0)/dr

-------------------------------
Inside scalar_profile.hpp you can modify 

u0 -- ampitude of scalar field
sigma -- width of scalar field 
rc -- radial coordinate of centre of scalar field
-------------------------------

         -*- -*- -*-

-------------------------------
- RUN THE CODE  
-------------------------------
To compile c++ code :

g++ main.cpp RK4Solver.cpp scalar_profile.cpp rhs.cpp -o cidc -std=c++17

-------------------------------
To run compiled code :

./cidc

-------------------------------
To graph plot the csv file use python : 

python3 plot.py

-------------------------------

         -*- -*- -*-
         
       +-------------+
       | robin croft | 
       +-------------+







