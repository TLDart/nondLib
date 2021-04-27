CXX = g++
CXX_FLAGS = -std=c++14 -O3
CXX_WARN_FLAGS = -Wall -Wextra

.PHONY: all
all: example2d example3d example5d

%: %.cpp ../nondlib.hpp
	$(CXX) $(CXX_FLAGS) $(CXX_WARN_FLAGS) $< -o $@
