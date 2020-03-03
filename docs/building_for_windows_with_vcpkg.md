# Building osgEarth for Windows using vcpkg

[vcpkg](https://github.com/Microsoft/vcpkg) is an extremly useful C++
package manager. It works on Windows, Linux and MacOS but for this guide
we'll focus on Windows.

### Install vcpkg and the osgEarth dependencies
First, download and bootstrap
[vcpkg](https://github.com/Microsoft/vcpkg) following the instructions on the page.

Next install the dependencies you'll need:
```
vcpkg install osg:x64-windows gdal:x64-windows geos:x64-windows sqlite3:x64-windows protobuf:x64-windows
```
This will take a while the first time you run it. It pulls down lots of dependencies, so go get a cup of coffee.

### Get the osgEarth source
Once all the dependencies are built, you'll need to actually build osgEarth. Clone the source code and make a folder for an out-of-source build (the only good kind of build):
```
git clone <https://github.com/gwaldron/osgearth.git> repo
mkdir build
cd build
```

### Configuring CMake
vcpkg provides a CMake toolchain file that helps osgEarth find all of
its dependencies. You'll need to specify a different build directory
for Release and Debug and specify the build type using
-DCMAKE\_BUILD\_TYPE. This is because some dependencies of osgEarth
don't pick up both debug and release versions without specifying the
build type. This should be fixed in future CMake versions. This is for a
release build. Run from the *build* folder you make earlier:
```
cmake ..\repo -G "Visual Studio 15 2017 Win64"
   -DCMAKE_BUILD_TYPE=Release -DWIN32_USE_MP=ON
   -DCMAKE_TOOLCHAIN_FILE=[vcpkg root]\scripts\buildsystems\vcpkg.cmake
```
Or you can just use the CMake UI.

### Build and install osgEarth

You can build and install osgEarth on the command line using cmake or
you can open up the Visual Studio solution and build it from there:
```
cmake --build . --target INSTALL --config Release
```
### Set up your runtime environment
You'll need to make sure that the vcpkg dependencies and osgEarth are in
your path. So do something like this:
```
set PATH=%PATH%;c:\vcpkg\installed\x64-windows\bin
set PATH=%PATH%;c:\vcpkg\installed\x64-windows\tools\osg
set PATH=%PATH%;c:\Program Files\osgEarth\bin

Note: If you don't want to build osgEarth yourself for your application,
you can actually install it using vcpkg as well. Just use:
```
vcpkg install osgearth:x64-windows
```
