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
include parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/depend.make

# Include the progress variables for this target.
include parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/progress.make

# Include the compile flags for this target's objects.
include parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/flags.make

parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o: parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/flags.make
parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o: ../parallel/parallel_src/app/edge_list_to_metis_graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o"
	cd /home/aangone/vsp/KaHIP/build/parallel/parallel_src && /software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o -c /home/aangone/vsp/KaHIP/parallel/parallel_src/app/edge_list_to_metis_graph.cpp

parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.i"
	cd /home/aangone/vsp/KaHIP/build/parallel/parallel_src && /software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aangone/vsp/KaHIP/parallel/parallel_src/app/edge_list_to_metis_graph.cpp > CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.i

parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.s"
	cd /home/aangone/vsp/KaHIP/build/parallel/parallel_src && /software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aangone/vsp/KaHIP/parallel/parallel_src/app/edge_list_to_metis_graph.cpp -o CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.s

# Object files for target edge_list_to_metis_graph
edge_list_to_metis_graph_OBJECTS = \
"CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o"

# External object files for target edge_list_to_metis_graph
edge_list_to_metis_graph_EXTERNAL_OBJECTS = \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/parallel_graph_access.cpp.o" \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/io/parallel_graph_io.cpp.o" \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management.cpp.o" \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management_refinement.cpp.o" \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management_coarsening.cpp.o" \
"/home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/libedgelist.dir/extern/argtable3-3.0.3/argtable3.c.o"

parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/app/edge_list_to_metis_graph.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/parallel_graph_access.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/io/parallel_graph_io.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management_refinement.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/lib/data_structure/balance_management_coarsening.cpp.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/libedgelist.dir/extern/argtable3-3.0.3/argtable3.c.o
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/build.make
parallel/parallel_src/edge_list_to_metis_graph: parallel/modified_kahip/liblibmodified_kahip_interface.a
parallel/parallel_src/edge_list_to_metis_graph: /software/spackages/linux-centos8-x86_64/gcc-8.3.1/openmpi-4.0.3-pncrqo3qrs323qe3lwg6pohfvix3xdha/lib/libmpi_cxx.so
parallel/parallel_src/edge_list_to_metis_graph: /software/spackages/linux-centos8-x86_64/gcc-8.3.1/openmpi-4.0.3-pncrqo3qrs323qe3lwg6pohfvix3xdha/lib/libmpi.so
parallel/parallel_src/edge_list_to_metis_graph: /software/spackages/linux-centos8-x86_64/gcc-8.3.1/gcc-9.3.0-s4u74h7rdafctph4bcbuaemterx5zhzc/lib64/libgomp.so
parallel/parallel_src/edge_list_to_metis_graph: /lib64/libpthread.so
parallel/parallel_src/edge_list_to_metis_graph: parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aangone/vsp/KaHIP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable edge_list_to_metis_graph"
	cd /home/aangone/vsp/KaHIP/build/parallel/parallel_src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/edge_list_to_metis_graph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/build: parallel/parallel_src/edge_list_to_metis_graph

.PHONY : parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/build

parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/clean:
	cd /home/aangone/vsp/KaHIP/build/parallel/parallel_src && $(CMAKE_COMMAND) -P CMakeFiles/edge_list_to_metis_graph.dir/cmake_clean.cmake
.PHONY : parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/clean

parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/depend:
	cd /home/aangone/vsp/KaHIP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aangone/vsp/KaHIP /home/aangone/vsp/KaHIP/parallel/parallel_src /home/aangone/vsp/KaHIP/build /home/aangone/vsp/KaHIP/build/parallel/parallel_src /home/aangone/vsp/KaHIP/build/parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : parallel/parallel_src/CMakeFiles/edge_list_to_metis_graph.dir/depend

