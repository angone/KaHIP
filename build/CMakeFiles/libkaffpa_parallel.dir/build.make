# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /software/spackages/linux-centos8-x86_64/gcc-8.3.1/cmake-3.17.3-3bmddqa2rihg4wv3x55bwjyhu6g4pciv/bin/cmake

# The command to remove a file.
RM = /software/spackages/linux-centos8-x86_64/gcc-8.3.1/cmake-3.17.3-3bmddqa2rihg4wv3x55bwjyhu6g4pciv/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aangone/vsp/KaHIP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aangone/vsp/KaHIP/build

# Include any dependencies generated for this target.
include CMakeFiles/libkaffpa_parallel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/libkaffpa_parallel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/libkaffpa_parallel.dir/flags.make

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.o: ../lib/parallel_mh/parallel_mh_async.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.o -c /home/aangone/vsp/KaHIP/lib/parallel_mh/parallel_mh_async.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/parallel_mh/parallel_mh_async.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/parallel_mh/parallel_mh_async.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.o: ../lib/parallel_mh/population.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.o -c /home/aangone/vsp/KaHIP/lib/parallel_mh/population.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/parallel_mh/population.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/parallel_mh/population.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.o: ../lib/parallel_mh/galinier_combine/gal_combine.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.o -c /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/gal_combine.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/gal_combine.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/gal_combine.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.o: ../lib/parallel_mh/galinier_combine/construct_partition.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.o -c /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/construct_partition.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/construct_partition.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/parallel_mh/galinier_combine/construct_partition.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.o: ../lib/parallel_mh/exchange/exchanger.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.o -c /home/aangone/vsp/KaHIP/lib/parallel_mh/exchange/exchanger.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/parallel_mh/exchange/exchanger.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/parallel_mh/exchange/exchanger.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.o: ../lib/tools/graph_communication.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.o -c /home/aangone/vsp/KaHIP/lib/tools/graph_communication.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/tools/graph_communication.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/tools/graph_communication.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.s

CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.o: CMakeFiles/libkaffpa_parallel.dir/flags.make
CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.o: ../lib/tools/mpi_tools.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.o"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.o -c /home/aangone/vsp/KaHIP/lib/tools/mpi_tools.cpp

CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.i"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/lib/tools/mpi_tools.cpp > CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.i

CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.s"
	/software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/lib/tools/mpi_tools.cpp -o CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.s

libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/parallel_mh_async.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/population.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/gal_combine.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/galinier_combine/construct_partition.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/parallel_mh/exchange/exchanger.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/tools/graph_communication.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/lib/tools/mpi_tools.cpp.o
libkaffpa_parallel: CMakeFiles/libkaffpa_parallel.dir/build.make

.PHONY : libkaffpa_parallel

# Rule to build all files generated by this target.
CMakeFiles/libkaffpa_parallel.dir/build: libkaffpa_parallel

.PHONY : CMakeFiles/libkaffpa_parallel.dir/build

CMakeFiles/libkaffpa_parallel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/libkaffpa_parallel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/libkaffpa_parallel.dir/clean

CMakeFiles/libkaffpa_parallel.dir/depend:
	cd /home/aangone/vsp/KaHIP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aangone/vsp/KaHIP /home/aangone/vsp/KaHIP /home/aangone/vsp/KaHIP/build /home/aangone/vsp/KaHIP/build /home/aangone/vsp/KaHIP/build/CMakeFiles/libkaffpa_parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/libkaffpa_parallel.dir/depend
