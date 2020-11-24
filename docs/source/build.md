# Building osgEarth

## Building with vcpkg

[vcpkg](https://github.com/Microsoft/vcpkg) is an extremely useful package manager. It works on Windows, Linux and MacOS but for this guide we will focus on Windows.

**Step 1 - Configure vcpkg**

First, download and bootstrap [vcpkg](https://github.com/Microsoft/vcpkg) following the instructions on the page.

Next install the dependencies required to build a fully functional osgEarth. This example assume s 64-bit WIndows build; you can alter that to correspond to your platform/architecture of choice.

Install the required dependencies:

```
vcpkg install osg:x64-windows gdal:x64-windows curl:x64-windows
```

For full functionality, you can install optional dependences as well:

```
vcpkg install sqlite3:x64-windows protobuf:x64-windows geos:x64-windows blend2d:x64-windows webp:x64-windows basisu:x64-windows draco:x64-windows libzip:x64-windows
```

This will take awhile the first time you run it as this pulls down lots of dependencies, so go get a cup of coffee.

Once all the dependencies are built, you’ll need to actually build osgEarth.

**Step 2 - Clone the repository**

Pull down the source from GitHub and create a ```build``` folder for your out-of-source build. We always recommend doing an out-of-source build to avoid problems down the road!

```
git clone https://github.com/gwaldron/osgearth.git
mkdir build
```

**Step 4 - Configure CMake**

vcpkg provides a CMake toolchain file that helps osgEarth find all of its dependencies.

Note: You’ll need to specify a different build directory based on your build configuration (Release, RelWIthDebInfo, Debug) and specify the build type using ```-DCMAKE_BUILD_TYPE```. This is because some dependencies of osgEarth don’t pick up both debug and release versions without specifying the build type. Hopefully this will be fixed in future CMake versions.

Most developers will use a RelWithDebInfo build, like so:

```
cmake \
	-S osgearth \
	-B build \
	-G "Visual Studio 15 2017 Win64" \
	-DCMAKE_BUILD_TYPE=RelWithDebInfo  \
	-DWIN32_USE_MP=ON \
	-DCMAKE_INSTALL_PREFIX=[folder to install osgEarth to] \
	-DCMAKE_TOOLCHAIN_FILE=[vcpkg root]\scripts\buildsystems\vcpkg.cmake
```

**Step 5 - Build and install osgEarth**

You can build and install osgEarth on the command line using CMake or you can open up the Visual Studio solution and build it from there.

```
cmake --build build --target INSTALL --config RelWithDebInfo
```

**Step 6 - Set up your runtime environment**

You’ll need to make sure that the vcpkg dependencies and osgEarth are in your path:

```
set PATH=%PATH%;c:\vcpkg\installed\x64-windows\bin
set PATH=%PATH%;c:\vcpkg\installed\x64-windows\tools\osg
set PATH=%PATH%;[folder where you installed osgEarth]
```

