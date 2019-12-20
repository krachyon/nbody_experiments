## prerequisites
### C++
- CMake
- Pybind11
- Eigen3
- Boost
- OpenMP

The last two aren't really used, so you could remove them in CMakeLists.txt

### Python
- Jupyter
- ipywidgets


To build the C++ part do something like this:
```
mkdir build_release
cd build_release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

Then there should be a file like `nbody.cpython-38-x86_64-linux-gnu.so` that you can then import into python with `import nbody` or like in the notebook `from build_release import nbody`

Then you should be able to open the notebook (e.g. with `jupyter notebook` and browser) and maybe have it work
