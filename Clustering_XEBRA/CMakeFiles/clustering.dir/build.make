# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA

# Include any dependencies generated for this target.
include CMakeFiles/clustering.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/clustering.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/clustering.dir/flags.make

CMakeFiles/clustering.dir/Clustering.cpp.o: CMakeFiles/clustering.dir/flags.make
CMakeFiles/clustering.dir/Clustering.cpp.o: Clustering.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/clustering.dir/Clustering.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/clustering.dir/Clustering.cpp.o -c /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/Clustering.cpp

CMakeFiles/clustering.dir/Clustering.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/clustering.dir/Clustering.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/Clustering.cpp > CMakeFiles/clustering.dir/Clustering.cpp.i

CMakeFiles/clustering.dir/Clustering.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/clustering.dir/Clustering.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/Clustering.cpp -o CMakeFiles/clustering.dir/Clustering.cpp.s

CMakeFiles/clustering.dir/Clustering.cpp.o.requires:

.PHONY : CMakeFiles/clustering.dir/Clustering.cpp.o.requires

CMakeFiles/clustering.dir/Clustering.cpp.o.provides: CMakeFiles/clustering.dir/Clustering.cpp.o.requires
	$(MAKE) -f CMakeFiles/clustering.dir/build.make CMakeFiles/clustering.dir/Clustering.cpp.o.provides.build
.PHONY : CMakeFiles/clustering.dir/Clustering.cpp.o.provides

CMakeFiles/clustering.dir/Clustering.cpp.o.provides.build: CMakeFiles/clustering.dir/Clustering.cpp.o


CMakeFiles/clustering.dir/AnalysisSort.C.o: CMakeFiles/clustering.dir/flags.make
CMakeFiles/clustering.dir/AnalysisSort.C.o: AnalysisSort.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/clustering.dir/AnalysisSort.C.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/clustering.dir/AnalysisSort.C.o -c /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/AnalysisSort.C

CMakeFiles/clustering.dir/AnalysisSort.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/clustering.dir/AnalysisSort.C.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/AnalysisSort.C > CMakeFiles/clustering.dir/AnalysisSort.C.i

CMakeFiles/clustering.dir/AnalysisSort.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/clustering.dir/AnalysisSort.C.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/AnalysisSort.C -o CMakeFiles/clustering.dir/AnalysisSort.C.s

CMakeFiles/clustering.dir/AnalysisSort.C.o.requires:

.PHONY : CMakeFiles/clustering.dir/AnalysisSort.C.o.requires

CMakeFiles/clustering.dir/AnalysisSort.C.o.provides: CMakeFiles/clustering.dir/AnalysisSort.C.o.requires
	$(MAKE) -f CMakeFiles/clustering.dir/build.make CMakeFiles/clustering.dir/AnalysisSort.C.o.provides.build
.PHONY : CMakeFiles/clustering.dir/AnalysisSort.C.o.provides

CMakeFiles/clustering.dir/AnalysisSort.C.o.provides.build: CMakeFiles/clustering.dir/AnalysisSort.C.o


# Object files for target clustering
clustering_OBJECTS = \
"CMakeFiles/clustering.dir/Clustering.cpp.o" \
"CMakeFiles/clustering.dir/AnalysisSort.C.o"

# External object files for target clustering
clustering_EXTERNAL_OBJECTS =

clustering: CMakeFiles/clustering.dir/Clustering.cpp.o
clustering: CMakeFiles/clustering.dir/AnalysisSort.C.o
clustering: CMakeFiles/clustering.dir/build.make
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libCore.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libRIO.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libNet.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libHist.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libGraf.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libGraf3d.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libGpad.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libTree.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libRint.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libPostscript.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libMatrix.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libPhysics.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libMathCore.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libThread.so
clustering: /usr/local/root_v6.12.04/root-6.12.04_build/lib/libMultiProc.so
clustering: CMakeFiles/clustering.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable clustering"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/clustering.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/clustering.dir/build: clustering

.PHONY : CMakeFiles/clustering.dir/build

CMakeFiles/clustering.dir/requires: CMakeFiles/clustering.dir/Clustering.cpp.o.requires
CMakeFiles/clustering.dir/requires: CMakeFiles/clustering.dir/AnalysisSort.C.o.requires

.PHONY : CMakeFiles/clustering.dir/requires

CMakeFiles/clustering.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/clustering.dir/cmake_clean.cmake
.PHONY : CMakeFiles/clustering.dir/clean

CMakeFiles/clustering.dir/depend:
	cd /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA /home/ab602/Thesis/Xebra_G4_Analysis/Clustering_XEBRA/CMakeFiles/clustering.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/clustering.dir/depend

