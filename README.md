# HGF Film Simulation

My implementation of ["A Hyperbolic Geometric Flow for Evolving Films and Foams" \[Sadashige Ishida et al. 2017\]](https://sadashigeishida.bitbucket.io/hgf/index.html).
![]("screenshot.png")

## Current Status
Implements only the core algorithms, and output the results on screen with wireframe rendering. You may want to change parameters directly in "src/simulator.cpp" to try different configurations.

## Dependencies
- Eigen
- GLFW
- glut
- LAPACK
- libigl
- OpenGL

## Compile
    
    mkdir build
    cd build
    cmake ..
    make

## Run

    cd build
    ./main [scene]

where the optional argument [scene] can be 1 to 4.

## References
- ["A Hyperbolic Geometric Flow for Evolving Films and Foams" \[Sadashige Ishida et al. 2017\]](https://sadashigeishida.bitbucket.io/hgf/index.html)
