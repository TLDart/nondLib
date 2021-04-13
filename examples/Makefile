CXX = g++
CXX_FLAGS = -std=c++14 -O3
CXX_WARN_FLAGS = -Wall -Wextra

.PHONY: all
all: example2d example3d example5d

example5d: example5d.cpp ../nondlib.hpp
	$(CXX) $(CXX_FLAGS) $(CXX_WARN_FLAGS) $< -o $@

example3d: example3d.cpp ../nondlib.hpp
	$(CXX) $(CXX_FLAGS) $(CXX_WARN_FLAGS) $< -o $@

example2d: example2d.cpp ../nondlib.hpp
	$(CXX) $(CXX_FLAGS) $(CXX_WARN_FLAGS) $< -o $@
