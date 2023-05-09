# traycer

![License](https://img.shields.io/github/license/fonzcastellanos/traycer)

A ray tracer, which generates a 2D image by emitting light rays and simulating the effects of their intersections with virtual 3D objects.

spheres.scene             | table.scene
:------------------------:|:-----------------------------------:
![](spheres.jpg)          | ![](table.jpg)
10 rays/pixel, 10 bounces | 5 rays/pixel, 10 extra lights/light

## Features
- Light rays intersect objects composed of
  - Triangles
  - Spheres
- Phong shading
- Anti-aliasing via jittered supersampling
- Recursive reflection
- Soft shadows

## Build Requirements
- Operating system: macOS or Linux
- C++ compiler supporting at least C++17
- CMake >= 3.2
- OpenGL >= 3.2

My development environment and tooling:
- Hardware: MacBook Pro (Retina, 13-inch, Early 2015)
- Operating system: macOS Mojave == 10.14.16
- C++ compiler: Clang == 10.0.01
- CMake == 3.25.1
- OpenGL == 4.1

## Build Steps

The build system must be generated before the project can be built. From the project directory, generate the build system:
```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE:STRING=Release ..
```
The `-G` option is omitted in `cmake -DCMAKE_BUILD_TYPE:STRING=Release ..`, so CMake will choose a default build system generator type based on your platform. To learn more about generators, see the [CMake docs](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html).

After having generated the build system, go back to the project directory:
```sh
cd ..
```

Build *all* the targets of the project:
```sh
cmake --build build --config Release
```

The built targets are placed in the directory `build`. There is only one executable target: `traycer`.

## Usage
`./build/traycer [options...] <scene-file>`
- `--rays-per-pixel <unsigned integer>`: If > 1, then supersampling is applied. Default is 1.
- `--reflection-bounces <unsigned integer>`: The depth of reflection cursion. Default is 0.
- `--extra-lights-per-light <unsigned integer>`: If > 0, soft shadows are applied. Default is 0. 
- `--render-file <filepath>`: Rendered image will be written to a JPEG file at the path provided.

The [`scenes`](scenes) directory contains example scene files. 