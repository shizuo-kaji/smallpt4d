smallpt4D: simple ray tracer for 4D scenes
==================
* written by S. Kaji
* Based on: smallpt: Global Illumination in 99 lines of C++ by Kevin Beason

# Build

Requirements:
* CMake 3.16+
* A C++17 compiler
* Eigen3 headers
* OpenMP runtime
* OpenCL runtime/SDK (optional, for GPU mode via `-g`)

```bash
cmake --build build --config Release
```

# Usage

Edit `scene_definition.cpp` to define a 4D scene.
Then,

```bash
./build/smallpt4d
```

Command-line options:

```bash
./build/smallpt4d [-m min_bounce] [-s samples] [-w width] [-h height] [-f frames] [-g]
```

* `-m <int>`: minimum bounce count before Russian roulette (default: `5`)
* `-s <int>`: samples per pixel (default: `100`)
* `-w <int>`: output width in pixels (default: `1080`)
* `-h <int>`: output height in pixels (default: `1200`)
* `-f <int>`: number of frames (default: `200`)
* `-g`: enable OpenCL GPU backend when available (sometimes quality degrades due to fp32 precision)

Combine rendered images into a side-by-side video for HMDs:

```bash
python3 create_video.py
```

Play the produced video in a media player such as [Whirligig](http://www.whirligig.xyz/)
(press `Y` to enter SBS mode in Whirligig).
