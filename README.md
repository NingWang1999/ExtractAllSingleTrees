# ExtractAllSingleTrees
Rough Extraction of All individual Trees from a forest/orchard.

## Build

The current implementation depends on the [Point Cloud Library (PCL)](https://pointclouds.org) and its dependencies
(e.g., [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and the [Visualization Toolkit (VTK)](https://vtk.org)).
Please install [PCL](https://pointclouds.org/downloads/#cross-platform) (and its dependencies) first.

To build this project, you need [CMake](https://cmake.org/download/) (`>= 3.12`) and a compiler that
supports `>= C++17`.
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
To process raw forest or orchard point cloud data and extract individual trees, run the built executable like this:
```commandline
./MakeDataset.exe <raw_data_folder> <rough_data_folder>
```
- `<raw_data_folder>` specifies the directory containing the original PCD files (e.g., LiDAR scans of forests).
- `<rough_data_folder>` specifies the output directory where the processed tree-extracted PCD files will be saved.
- This will process all `.pcd` files in the `raw_data_folder` directory, extract individual trees, and save the results in `rough_data_folder`.

## References & Acknowledgments
This project is based on the original open-source project [GlobalMatch](https://github.com/zexinyang/GlobalMatch). Many thanks to the original author for their valuable contributions!


