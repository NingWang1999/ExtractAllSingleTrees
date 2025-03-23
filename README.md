# ExtractAllSingleTrees
Extract all individual trees from a forest/orchard pointclouds.

## Build

The current implementation depends on the [Point Cloud Library (PCL)](https://pointclouds.org) and its dependencies
(e.g., [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and the [Visualization Toolkit (VTK)](https://vtk.org)).
Please install [PCL](https://pointclouds.org/downloads/#cross-platform) (and its dependencies) first.

To build this project, you need [CMake](https://cmake.org/download/) (`>= 3.12`) and a compiler that
supports `>= C++14`.
With CMake, this project can be built on almost all platforms, 
although so far we have only tested it on Windows.


- On Windows with Microsoft Visual Studio, use the `x64 Native Tools Command Prompt for VS XXXX` (**don't** use the
  x86 one), then
  ```
  $ cd path-to-root-dir-of-this project
  $ mkdir Release
  $ cd Release
  $ cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
  $ ninja
  ```

## Usage

The [main.cpp](./makedataset/main.cpp) file demonstrates how to use it. 
To register a "source" forest point cloud to a "target" one, run the built executable like this:
```commandline
./MakeDataset.exe  <source_scan_ply_file>  <target_scan_ply_file>
```
- `<source_scan_ply_file>` and `<target_scan_ply_file>` specify the path to the input point cloud **ply** files (**order matters**: source filename comes first)
- `<output_matrix_txt_file>` specifies the filename to save the estimated 4 by 4 transformation matrix




