# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build

# Include any dependencies generated for this target.
include CMakeFiles/MakeMatchingHits.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/MakeMatchingHits.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/MakeMatchingHits.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MakeMatchingHits.dir/flags.make

CMakeFiles/MakeMatchingHits.dir/main.cpp.o: CMakeFiles/MakeMatchingHits.dir/flags.make
CMakeFiles/MakeMatchingHits.dir/main.cpp.o: /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/main.cpp
CMakeFiles/MakeMatchingHits.dir/main.cpp.o: CMakeFiles/MakeMatchingHits.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MakeMatchingHits.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MakeMatchingHits.dir/main.cpp.o -MF CMakeFiles/MakeMatchingHits.dir/main.cpp.o.d -o CMakeFiles/MakeMatchingHits.dir/main.cpp.o -c /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/main.cpp

CMakeFiles/MakeMatchingHits.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MakeMatchingHits.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/main.cpp > CMakeFiles/MakeMatchingHits.dir/main.cpp.i

CMakeFiles/MakeMatchingHits.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MakeMatchingHits.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/main.cpp -o CMakeFiles/MakeMatchingHits.dir/main.cpp.s

CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o: CMakeFiles/MakeMatchingHits.dir/flags.make
CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o: /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/src/FakeSFG.cxx
CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o: CMakeFiles/MakeMatchingHits.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o -MF CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o.d -o CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o -c /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/src/FakeSFG.cxx

CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/src/FakeSFG.cxx > CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.i

CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/src/FakeSFG.cxx -o CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.s

# Object files for target MakeMatchingHits
MakeMatchingHits_OBJECTS = \
"CMakeFiles/MakeMatchingHits.dir/main.cpp.o" \
"CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o"

# External object files for target MakeMatchingHits
MakeMatchingHits_EXTERNAL_OBJECTS =

MakeMatchingHits: CMakeFiles/MakeMatchingHits.dir/main.cpp.o
MakeMatchingHits: CMakeFiles/MakeMatchingHits.dir/src/FakeSFG.cxx.o
MakeMatchingHits: CMakeFiles/MakeMatchingHits.dir/build.make
MakeMatchingHits: libFakeSFG.dylib
MakeMatchingHits: /opt/root_install/lib/libCore.so
MakeMatchingHits: /opt/root_install/lib/libImt.so
MakeMatchingHits: /opt/root_install/lib/libRIO.so
MakeMatchingHits: /opt/root_install/lib/libNet.so
MakeMatchingHits: /opt/root_install/lib/libHist.so
MakeMatchingHits: /opt/root_install/lib/libGraf.so
MakeMatchingHits: /opt/root_install/lib/libGraf3d.so
MakeMatchingHits: /opt/root_install/lib/libGpad.so
MakeMatchingHits: /opt/root_install/lib/libROOTDataFrame.so
MakeMatchingHits: /opt/root_install/lib/libTree.so
MakeMatchingHits: /opt/root_install/lib/libTreePlayer.so
MakeMatchingHits: /opt/root_install/lib/libRint.so
MakeMatchingHits: /opt/root_install/lib/libPostscript.so
MakeMatchingHits: /opt/root_install/lib/libMatrix.so
MakeMatchingHits: /opt/root_install/lib/libPhysics.so
MakeMatchingHits: /opt/root_install/lib/libMathCore.so
MakeMatchingHits: /opt/root_install/lib/libThread.so
MakeMatchingHits: /opt/root_install/lib/libMultiProc.so
MakeMatchingHits: /opt/root_install/lib/libROOTVecOps.so
MakeMatchingHits: CMakeFiles/MakeMatchingHits.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable MakeMatchingHits"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MakeMatchingHits.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MakeMatchingHits.dir/build: MakeMatchingHits
.PHONY : CMakeFiles/MakeMatchingHits.dir/build

CMakeFiles/MakeMatchingHits.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MakeMatchingHits.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MakeMatchingHits.dir/clean

CMakeFiles/MakeMatchingHits.dir/depend:
	cd /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build /Users/lorenzo/T2K-uniGe/ND280Software/FakeSFG/build/CMakeFiles/MakeMatchingHits.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/MakeMatchingHits.dir/depend
