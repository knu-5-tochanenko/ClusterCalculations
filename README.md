# How to create and run MPI projects in CLion

* Download and install MPI from [Microsoft Website](https://www.microsoft.com/en-us/download/details.aspx?id=100593). You need to download both `msmpisetup.exe` and `msmpisdk.msi`.

> Note that you shouldn't have spacec or cyrillic characters in path to MPI installation.

1. Open Windows Terminal and try running `mpiexec`:

```bash
Microsoft MPI Startup Program [Version 10.1.12498.18]

...
```

If the `mpiexec` command isn't recognised by your system try restarting your PC and check if MPI system variables are installed correctly:

| Name | Value |
| --- | --- |
| MSMPI_BENCHMARKS | `path_to_mpi\Benchmarks\` |
| MSMPI_BIN | `path_to_mpi\Bin\` |
| MSMPI_INC | `path_to_mpi\SDK\Include\` |
| MSMPI_LIB32 | `path_to_mpi\SDK\Lib\x86\` |
| MSMPI_LIB64 | `path_to_mpi\SDK\Lib\x64\` |

And `path_to_mpi\Bin\` path is added to `Path` variable.

> To view system variables open Windows Search -> Input "Variables" -> Click on "Edit the system environment variables" -> Click on "Environemnt Variables on "System Properties" Window.

> `path_to_mpi` is a path to your Microsoft MPI installation.


2. Open CLion and create new `C++ Executable` Project.

3. Open `CMakeLists.txt` where you should see something like this: (`MY_PROJECT` is project name)

```cmake
cmake_minimum_required(VERSION 3.20)
project(MY_PROJECT)

set(CMAKE_CXX_STANDARD 14)

add_executable(MY_PROJECT main.cpp)

```

4. Add MPI Package and MPI Library to project

```cmake
cmake_minimum_required(VERSION 3.20)
project (MY_PROJECT)

set(CMAKE_CXX_STANDARD 14)

# Add MPI Package to Project
find_package(MPI REQUIRED)

add_executable(MY_PROJECT main.cpp)
# Add libraries for code completion and compiling
target_link_libraries(MY_PROJECT PUBLIC MPI::MPI_CXX)
```

5. Add following run configurations:

**1) CMake Application**

| Property | Value |
| --- | --- |
| Name | `Compile` |
| Tagret | `MY_PROJECT` |
| Executable | `MY_PROJECT` |
| Environment Variables | `MPI_CXX_INCLUDE_PATH=path_to_mpi\SDK\Include\;MPI_CXX_LIBRARIES=path_to_mpi\SDK\Lib\x64\;MPI_CXX_BIN_PATH=path_to_mpi\Bin\` |

> IMPORTANT! Replace `path_to_mpi` string with a path to your Microsoft MPI intallation!

**2) Shell Script**

| Property | Value |
| --- | --- |
| Name | `2 Cores` |
| Execute | `Script text` |
| Script text | `mpiexec -n 2 ./cmake-build-debug/MY_PROJECT` |
| Working directory | `project_directory/MY_PROJECT` |
| Execute in the termial | `true` |

> Add more shell scripts for number of cores you would like to use while executing your project. Replace number `2` in `Name` and `Script text` properties with number of cores you would like to use.

> `MY_PROJECT` is a name of `.exe` file in `./cmake-build-debug` folder.

6. Run `Compile` configuration before executing any of `X Cores` configurations.