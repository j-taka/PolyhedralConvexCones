# PolyhedralConvexCones
Calculate dual of polyhedral convex cones.
In detail, given the set of vectors a, the algorithm outputs the set of vectors b, vice vasa.

<img src="https://latex.codecogs.com/gif.latex?\bigcap_{i}&space;\bf{a}_i^{T}&space;\bf{x}&space;\geq&space;0&space;\leftrightarrow&space;\sum&space;c_i&space;\bf{b}_i&space;=&space;\bf{x}&space;(\forall&space;i,&space;c_i&space;\geq&space;0)" />

This software is used to classify the contact states in Assembly Plan from Observation:
[1] Jun Takamatsu, Koichi Ogawara, Hiroshi Kimura, and Katsushi ikeuchi, “Recognizing Assembly Tasks Through Human Demonstration,” Int. J. of Robotics Research , Vol. 26, No. 7, pp. 641 – 659, Jul. 2007.

## Requirement
- CMake
- Qhull (http://www.qhull.org/)
- Eigen (http://www.eigen.tuxfamily.org/)

## Compile
1. Downloads all the requirements
2. Compile Qhull and set the environmental variable QHULL_PATH to the directory installed
3. Install Eigen and set the environmental variable Eigen3_DIR to the directory installed
4. Use cmake to compile the program

## Author
Jun Takamatsu <j-taka@is.naist.jp>

## License
[MIT](http://b4b4r07.mit-license.org)

