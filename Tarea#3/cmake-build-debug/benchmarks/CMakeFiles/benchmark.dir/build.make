# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2017.3.3/bin/cmake/bin/cmake

# The command to remove a file.
RM = /opt/clion-2017.3.3/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug"

# Include any dependencies generated for this target.
include benchmarks/CMakeFiles/benchmark.dir/depend.make

# Include the progress variables for this target.
include benchmarks/CMakeFiles/benchmark.dir/progress.make

# Include the compile flags for this target's objects.
include benchmarks/CMakeFiles/benchmark.dir/flags.make

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o: benchmarks/CMakeFiles/benchmark.dir/flags.make
benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o: ../benchmarks/benchmarkMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmark.dir/benchmarkMain.cpp.o -c "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMain.cpp"

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark.dir/benchmarkMain.cpp.i"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMain.cpp" > CMakeFiles/benchmark.dir/benchmarkMain.cpp.i

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark.dir/benchmarkMain.cpp.s"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMain.cpp" -o CMakeFiles/benchmark.dir/benchmarkMain.cpp.s

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.requires:

.PHONY : benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.requires

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.provides: benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.requires
	$(MAKE) -f benchmarks/CMakeFiles/benchmark.dir/build.make benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.provides.build
.PHONY : benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.provides

benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.provides.build: benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o


benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o: benchmarks/CMakeFiles/benchmark.dir/flags.make
benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o: ../benchmarks/benchmarkMatrixAdd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o -c "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMatrixAdd.cpp"

benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.i"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMatrixAdd.cpp" > CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.i

benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.s"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks/benchmarkMatrixAdd.cpp" -o CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.s

benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.requires:

.PHONY : benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.requires

benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.provides: benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.requires
	$(MAKE) -f benchmarks/CMakeFiles/benchmark.dir/build.make benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.provides.build
.PHONY : benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.provides

benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.provides.build: benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o


# Object files for target benchmark
benchmark_OBJECTS = \
"CMakeFiles/benchmark.dir/benchmarkMain.cpp.o" \
"CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o"

# External object files for target benchmark
benchmark_EXTERNAL_OBJECTS =

benchmarks/benchmark: benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o
benchmarks/benchmark: benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o
benchmarks/benchmark: benchmarks/CMakeFiles/benchmark.dir/build.make
benchmarks/benchmark: src/libanpi.a
benchmarks/benchmark: /usr/local/lib/libboost_filesystem.so
benchmarks/benchmark: /usr/local/lib/libboost_system.so
benchmarks/benchmark: /usr/local/lib/libboost_unit_test_framework.so
benchmarks/benchmark: benchmarks/CMakeFiles/benchmark.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable benchmark"
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmark.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
benchmarks/CMakeFiles/benchmark.dir/build: benchmarks/benchmark

.PHONY : benchmarks/CMakeFiles/benchmark.dir/build

benchmarks/CMakeFiles/benchmark.dir/requires: benchmarks/CMakeFiles/benchmark.dir/benchmarkMain.cpp.o.requires
benchmarks/CMakeFiles/benchmark.dir/requires: benchmarks/CMakeFiles/benchmark.dir/benchmarkMatrixAdd.cpp.o.requires

.PHONY : benchmarks/CMakeFiles/benchmark.dir/requires

benchmarks/CMakeFiles/benchmark.dir/clean:
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" && $(CMAKE_COMMAND) -P CMakeFiles/benchmark.dir/cmake_clean.cmake
.PHONY : benchmarks/CMakeFiles/benchmark.dir/clean

benchmarks/CMakeFiles/benchmark.dir/depend:
	cd "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3" "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/benchmarks" "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug" "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks" "/home/nabetse28/Documentos/ANPI-TEC/Tarea#3/cmake-build-debug/benchmarks/CMakeFiles/benchmark.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : benchmarks/CMakeFiles/benchmark.dir/depend

