smallpt4D: simple ray tracer for 4D scenes
==================
* written by S. Kaji
* Based on: smallpt: Global Illumination in 99 lines of C++ by Kevin Beason
* Depends on: Eigen 3 (http://eigen.tuxfamily.org)
             stb_image_write.h by Sean Barrett (included)



# Usage

Edit main.cpp, especially around the end of the file to define a 4D scene.

    > vi main.cpp

Compile with gcc by

    > make

You may want to limit the number of threads used for rendering with OpenMP by

    > export OMP_NUM_THREADS=8

Rendering takes hours!

    > ./smallpt4d

Combine rendered images into a Side-by-side video for HMDs.

    > bash create-video.sh

Play the produced video by a media player such as [Whirligig>http://www.whirligig.xyz/]
(Press "Y" key to enter SBS mode in Whirligig)

