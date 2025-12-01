# Compiler, MODIFIED
CXX = g++
CXXFLAGS = -std=c++17 -O2 -fPIC -Iinclude -I/usr/local/include

# Source and object files
SRC = src/drag.cpp
OBJ = $(SRC:.cpp=.o)

# Output shared library
LIB = libdrag.so

# Default target
all: $(LIB)

# Compile source to object
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link shared library, MODIFIED
$(LIB): $(OBJ)
	$(CXX) -shared -o $@ $^ -lGeographicLib -L/usr/local/lib -lnrlmsise00

# Clean build files
clean:
	rm -f $(OBJ) $(LIB)

# Install target (flat layout)
install: $(LIB)
# Create include and lib directories
	mkdir -p /opt/drag/include
	mkdir -p /opt/drag/lib

# Copy the shared library
	cp $(LIB) /opt/drag/lib/

# Copy public header(s) to flat include
	cp include/drag.H /opt/drag/include/

# Copy Eigen headers
	cp -r include/Eigen /opt/drag/include/
