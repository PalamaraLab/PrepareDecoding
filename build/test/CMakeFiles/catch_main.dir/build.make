# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/abhidg/src/PrepareDecoding

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/abhidg/src/PrepareDecoding/build

# Include any dependencies generated for this target.
include test/CMakeFiles/catch_main.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/catch_main.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/catch_main.dir/flags.make

test/CMakeFiles/catch_main.dir/catch_main.cpp.o: test/CMakeFiles/catch_main.dir/flags.make
test/CMakeFiles/catch_main.dir/catch_main.cpp.o: ../test/catch_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhidg/src/PrepareDecoding/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/catch_main.dir/catch_main.cpp.o"
	cd /home/abhidg/src/PrepareDecoding/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/catch_main.dir/catch_main.cpp.o -c /home/abhidg/src/PrepareDecoding/test/catch_main.cpp

test/CMakeFiles/catch_main.dir/catch_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/catch_main.dir/catch_main.cpp.i"
	cd /home/abhidg/src/PrepareDecoding/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhidg/src/PrepareDecoding/test/catch_main.cpp > CMakeFiles/catch_main.dir/catch_main.cpp.i

test/CMakeFiles/catch_main.dir/catch_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/catch_main.dir/catch_main.cpp.s"
	cd /home/abhidg/src/PrepareDecoding/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhidg/src/PrepareDecoding/test/catch_main.cpp -o CMakeFiles/catch_main.dir/catch_main.cpp.s

# Object files for target catch_main
catch_main_OBJECTS = \
"CMakeFiles/catch_main.dir/catch_main.cpp.o"

# External object files for target catch_main
catch_main_EXTERNAL_OBJECTS =

test/libcatch_main.a: test/CMakeFiles/catch_main.dir/catch_main.cpp.o
test/libcatch_main.a: test/CMakeFiles/catch_main.dir/build.make
test/libcatch_main.a: test/CMakeFiles/catch_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/abhidg/src/PrepareDecoding/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcatch_main.a"
	cd /home/abhidg/src/PrepareDecoding/build/test && $(CMAKE_COMMAND) -P CMakeFiles/catch_main.dir/cmake_clean_target.cmake
	cd /home/abhidg/src/PrepareDecoding/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/catch_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/catch_main.dir/build: test/libcatch_main.a

.PHONY : test/CMakeFiles/catch_main.dir/build

test/CMakeFiles/catch_main.dir/clean:
	cd /home/abhidg/src/PrepareDecoding/build/test && $(CMAKE_COMMAND) -P CMakeFiles/catch_main.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/catch_main.dir/clean

test/CMakeFiles/catch_main.dir/depend:
	cd /home/abhidg/src/PrepareDecoding/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abhidg/src/PrepareDecoding /home/abhidg/src/PrepareDecoding/test /home/abhidg/src/PrepareDecoding/build /home/abhidg/src/PrepareDecoding/build/test /home/abhidg/src/PrepareDecoding/build/test/CMakeFiles/catch_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/catch_main.dir/depend

