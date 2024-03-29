# traycer

![License](https://img.shields.io/github/license/fonzcastellanos/traycer)

A ray tracer with Phong shading, supersampling, recursive reflections, and soft shadows.

spheres.scene             | table.scene
:------------------------:|:-----------------------------------:
![](spheres.jpg)          | ![](table.jpg)
10 jittered rays/pixel, 10 reflection bounces | 5 jittered rays/pixel, 10 extra lights/light (for producing soft shadows)

## Features
- Ray-surface intersections with the following surfaces:
  - Triangles
  - Spheres
- Anti-aliasing via jittered supersampling
- Phong shading
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

I built `traycer` on a Linux machine as well:
- Operating system: Debian == 11.7
- C++ compiler: gcc == 10.2.1
- CMake == 3.26.3
- OpenGL == 4.5

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

Options:
- `--jitter <number-of-rays-per-pixel>`
  - Enables supersampling via jittering if option argument is greater than 0. 
  - The default option argument is 0.
- `--bounces <number-of-reflection-bounces>`
  - The depth of reflection recursion.
  - The default option argument is 0.
- `--soft-shadows <number-of-extra-lights-per-light>`
  - Soft shadows are applied if option argument is greater than 0. 
  - The default option argument is 0. 
- `--render-to-file <filepath>`
  - Rendered image will be written to a JPEG file at the path provided.

To exit the program, have the window in focus and press `ESC`. You can also terminate the program by pressing `CTRL + C` in the terminal.

The [`scenes`](scenes) directory contains example scene files. 

### Example
```sh
./build/traycer --jitter 10 --bounces 5 --soft-shadows 5 --render-to-file scene.jpg scenes/spheres.scene
```