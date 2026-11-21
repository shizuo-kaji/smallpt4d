# smallpt4D · C++ renderer

The offline path tracer follows light through four-dimensional space and writes a numbered pair of left- and right-eye PNGs for every frame. The [project README](../README.md) introduces the optics and the interactive WebXR app.

## Build

Run the following commands from the `cpp/` directory.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

CMake 3.16+, a C++17 compiler, and Eigen3 headers are required. CMake enables OpenMP and OpenCL when available. If Eigen is installed in a custom location, set `EIGEN3_INCLUDE_DIR` to the directory containing `Eigen/Dense`. The [Makefile](Makefile) runs the same configure/build steps with `make`.

## Render

This example renders one low-resolution stereo pair in a separate output directory. Use a fresh directory for each run to keep earlier frames.

```sh
mkdir -p build/preview
cd build/preview
../smallpt4d -w 320 -h 240 -s 4 -f 1
```

Multi-configuration generators such as Visual Studio put the Release executable in `build/Release/`. From `build/preview/`, use `../Release/smallpt4d.exe` on Windows. The renderer writes `output_L_000.png` and `output_R_000.png` in the working directory.

| Option | Meaning | Default |
| --- | --- | --- |
| `-w`, `-h` | Width and height **per eye** | 1080 × 1200 |
| `-s` | Samples per pixel | 100 |
| `-f` | Number of stereo frames | 200 |
| `-m` | Minimum bounce count before Russian roulette | 5 |
| `-g` | Request OpenCL rendering, falling back to CPU if unavailable. Single-precision GPU results can differ from the CPU. | CPU |

To change a hypersphere's center, radius, color, emission, or material, edit the [scene definition](scene_definition.cpp) and rebuild. The default scene contains 625 lattice hyperspheres and eight enclosing hyperspheres, including walls in the two W directions. The [camera and CPU integrator](main.cpp) and [OpenCL renderer](opencl_renderer.cpp) trace reflection, refraction and diffuse scattering in four dimensions.

## Stereo video

For a moving sequence, use `-f` to render multiple frames into an empty directory, then run the [video helper](create_video.py) from that directory. It requires Python 3 and FFmpeg with the libx264 encoder.

```sh
# From cpp/build/preview/
python3 ../../create_video.py --fps 30 --output sbs.mp4
```

Play `sbs.mp4` in side-by-side mode in a stereo video player. For example, select SBS with `Y` in Whirligig. [Prerecorded examples](../demo/)
