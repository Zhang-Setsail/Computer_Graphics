# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = "/home/cs18/Desktop/ASS34/Assignment3&4"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/cs18/Desktop/ASS34/Assignment3&4"

# Include any dependencies generated for this target.
include CGL/CMakeFiles/CGL.dir/depend.make

# Include the progress variables for this target.
include CGL/CMakeFiles/CGL.dir/progress.make

# Include the compile flags for this target's objects.
include CGL/CMakeFiles/CGL.dir/flags.make

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o: CGL/src/vector2D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/vector2D.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector2D.cpp"

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/vector2D.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector2D.cpp" > CMakeFiles/CGL.dir/src/vector2D.cpp.i

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/vector2D.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector2D.cpp" -o CMakeFiles/CGL.dir/src/vector2D.cpp.s

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o


CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o: CGL/src/vector3D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/vector3D.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector3D.cpp"

CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/vector3D.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector3D.cpp" > CMakeFiles/CGL.dir/src/vector3D.cpp.i

CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/vector3D.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector3D.cpp" -o CMakeFiles/CGL.dir/src/vector3D.cpp.s

CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o


CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o: CGL/src/vector4D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/vector4D.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector4D.cpp"

CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/vector4D.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector4D.cpp" > CMakeFiles/CGL.dir/src/vector4D.cpp.i

CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/vector4D.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/vector4D.cpp" -o CMakeFiles/CGL.dir/src/vector4D.cpp.s

CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o


CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o: CGL/src/matrix3x3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/matrix3x3.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix3x3.cpp"

CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/matrix3x3.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix3x3.cpp" > CMakeFiles/CGL.dir/src/matrix3x3.cpp.i

CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/matrix3x3.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix3x3.cpp" -o CMakeFiles/CGL.dir/src/matrix3x3.cpp.s

CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o


CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o: CGL/src/matrix4x4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/matrix4x4.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix4x4.cpp"

CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/matrix4x4.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix4x4.cpp" > CMakeFiles/CGL.dir/src/matrix4x4.cpp.i

CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/matrix4x4.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/matrix4x4.cpp" -o CMakeFiles/CGL.dir/src/matrix4x4.cpp.s

CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o


CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o: CGL/src/quaternion.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/quaternion.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/quaternion.cpp"

CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/quaternion.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/quaternion.cpp" > CMakeFiles/CGL.dir/src/quaternion.cpp.i

CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/quaternion.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/quaternion.cpp" -o CMakeFiles/CGL.dir/src/quaternion.cpp.s

CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o


CGL/CMakeFiles/CGL.dir/src/complex.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/complex.cpp.o: CGL/src/complex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CGL/CMakeFiles/CGL.dir/src/complex.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/complex.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/complex.cpp"

CGL/CMakeFiles/CGL.dir/src/complex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/complex.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/complex.cpp" > CMakeFiles/CGL.dir/src/complex.cpp.i

CGL/CMakeFiles/CGL.dir/src/complex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/complex.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/complex.cpp" -o CMakeFiles/CGL.dir/src/complex.cpp.s

CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/complex.cpp.o


CGL/CMakeFiles/CGL.dir/src/color.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/color.cpp.o: CGL/src/color.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CGL/CMakeFiles/CGL.dir/src/color.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/color.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/color.cpp"

CGL/CMakeFiles/CGL.dir/src/color.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/color.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/color.cpp" > CMakeFiles/CGL.dir/src/color.cpp.i

CGL/CMakeFiles/CGL.dir/src/color.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/color.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/color.cpp" -o CMakeFiles/CGL.dir/src/color.cpp.s

CGL/CMakeFiles/CGL.dir/src/color.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/color.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/color.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/color.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/color.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/color.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/color.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/color.cpp.o


CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o: CGL/src/osdtext.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/osdtext.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdtext.cpp"

CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/osdtext.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdtext.cpp" > CMakeFiles/CGL.dir/src/osdtext.cpp.i

CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/osdtext.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdtext.cpp" -o CMakeFiles/CGL.dir/src/osdtext.cpp.s

CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o


CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o: CGL/src/osdfont.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/osdfont.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdfont.cpp"

CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/osdfont.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdfont.cpp" > CMakeFiles/CGL.dir/src/osdfont.cpp.i

CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/osdfont.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/osdfont.cpp" -o CMakeFiles/CGL.dir/src/osdfont.cpp.s

CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o


CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o: CGL/src/viewer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/viewer.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/viewer.cpp"

CGL/CMakeFiles/CGL.dir/src/viewer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/viewer.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/viewer.cpp" > CMakeFiles/CGL.dir/src/viewer.cpp.i

CGL/CMakeFiles/CGL.dir/src/viewer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/viewer.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/viewer.cpp" -o CMakeFiles/CGL.dir/src/viewer.cpp.s

CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o


CGL/CMakeFiles/CGL.dir/src/base64.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/base64.cpp.o: CGL/src/base64.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CGL/CMakeFiles/CGL.dir/src/base64.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/base64.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/base64.cpp"

CGL/CMakeFiles/CGL.dir/src/base64.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/base64.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/base64.cpp" > CMakeFiles/CGL.dir/src/base64.cpp.i

CGL/CMakeFiles/CGL.dir/src/base64.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/base64.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/base64.cpp" -o CMakeFiles/CGL.dir/src/base64.cpp.s

CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/base64.cpp.o


CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o: CGL/src/lodepng.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/lodepng.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/lodepng.cpp"

CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/lodepng.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/lodepng.cpp" > CMakeFiles/CGL.dir/src/lodepng.cpp.i

CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/lodepng.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/lodepng.cpp" -o CMakeFiles/CGL.dir/src/lodepng.cpp.s

CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o


CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o: CGL/src/tinyxml2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/tinyxml2.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/tinyxml2.cpp"

CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/tinyxml2.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/tinyxml2.cpp" > CMakeFiles/CGL.dir/src/tinyxml2.cpp.i

CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/tinyxml2.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/tinyxml2.cpp" -o CMakeFiles/CGL.dir/src/tinyxml2.cpp.s

CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o


CGL/CMakeFiles/CGL.dir/src/path.cpp.o: CGL/CMakeFiles/CGL.dir/flags.make
CGL/CMakeFiles/CGL.dir/src/path.cpp.o: CGL/src/path.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CGL/CMakeFiles/CGL.dir/src/path.cpp.o"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CGL.dir/src/path.cpp.o -c "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/path.cpp"

CGL/CMakeFiles/CGL.dir/src/path.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGL.dir/src/path.cpp.i"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/path.cpp" > CMakeFiles/CGL.dir/src/path.cpp.i

CGL/CMakeFiles/CGL.dir/src/path.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGL.dir/src/path.cpp.s"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/src/path.cpp" -o CMakeFiles/CGL.dir/src/path.cpp.s

CGL/CMakeFiles/CGL.dir/src/path.cpp.o.requires:

.PHONY : CGL/CMakeFiles/CGL.dir/src/path.cpp.o.requires

CGL/CMakeFiles/CGL.dir/src/path.cpp.o.provides: CGL/CMakeFiles/CGL.dir/src/path.cpp.o.requires
	$(MAKE) -f CGL/CMakeFiles/CGL.dir/build.make CGL/CMakeFiles/CGL.dir/src/path.cpp.o.provides.build
.PHONY : CGL/CMakeFiles/CGL.dir/src/path.cpp.o.provides

CGL/CMakeFiles/CGL.dir/src/path.cpp.o.provides.build: CGL/CMakeFiles/CGL.dir/src/path.cpp.o


# Object files for target CGL
CGL_OBJECTS = \
"CMakeFiles/CGL.dir/src/vector2D.cpp.o" \
"CMakeFiles/CGL.dir/src/vector3D.cpp.o" \
"CMakeFiles/CGL.dir/src/vector4D.cpp.o" \
"CMakeFiles/CGL.dir/src/matrix3x3.cpp.o" \
"CMakeFiles/CGL.dir/src/matrix4x4.cpp.o" \
"CMakeFiles/CGL.dir/src/quaternion.cpp.o" \
"CMakeFiles/CGL.dir/src/complex.cpp.o" \
"CMakeFiles/CGL.dir/src/color.cpp.o" \
"CMakeFiles/CGL.dir/src/osdtext.cpp.o" \
"CMakeFiles/CGL.dir/src/osdfont.cpp.o" \
"CMakeFiles/CGL.dir/src/viewer.cpp.o" \
"CMakeFiles/CGL.dir/src/base64.cpp.o" \
"CMakeFiles/CGL.dir/src/lodepng.cpp.o" \
"CMakeFiles/CGL.dir/src/tinyxml2.cpp.o" \
"CMakeFiles/CGL.dir/src/path.cpp.o"

# External object files for target CGL
CGL_EXTERNAL_OBJECTS =

CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/complex.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/color.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/base64.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/src/path.cpp.o
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/build.make
CGL/libCGL.a: CGL/CMakeFiles/CGL.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/cs18/Desktop/ASS34/Assignment3&4/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_16) "Linking CXX static library libCGL.a"
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && $(CMAKE_COMMAND) -P CMakeFiles/CGL.dir/cmake_clean_target.cmake
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CGL.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CGL/CMakeFiles/CGL.dir/build: CGL/libCGL.a

.PHONY : CGL/CMakeFiles/CGL.dir/build

CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/vector2D.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/vector3D.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/vector4D.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/matrix3x3.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/matrix4x4.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/quaternion.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/complex.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/color.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/osdtext.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/osdfont.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/viewer.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/base64.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/lodepng.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/tinyxml2.cpp.o.requires
CGL/CMakeFiles/CGL.dir/requires: CGL/CMakeFiles/CGL.dir/src/path.cpp.o.requires

.PHONY : CGL/CMakeFiles/CGL.dir/requires

CGL/CMakeFiles/CGL.dir/clean:
	cd "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" && $(CMAKE_COMMAND) -P CMakeFiles/CGL.dir/cmake_clean.cmake
.PHONY : CGL/CMakeFiles/CGL.dir/clean

CGL/CMakeFiles/CGL.dir/depend:
	cd "/home/cs18/Desktop/ASS34/Assignment3&4" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/cs18/Desktop/ASS34/Assignment3&4" "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" "/home/cs18/Desktop/ASS34/Assignment3&4" "/home/cs18/Desktop/ASS34/Assignment3&4/CGL" "/home/cs18/Desktop/ASS34/Assignment3&4/CGL/CMakeFiles/CGL.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CGL/CMakeFiles/CGL.dir/depend

