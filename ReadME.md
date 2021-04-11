# A Non Dominated Points Library

This Repository is a dedicated to keeping a Library for NonDominated Points.

## Important Remarks

* All of the algorithms listed in this library modify the initial set of Points. Meaning that any dominated point will be lost forever after the algorithm runs.
* All of the algorithms are non-stable. Meaning that after running the algorithms in order in relation to the output is maintained.
* In the current version of the algorithm the original vector is modified to a vector without non-dominated points.

## General Usage

To use this library you copy the "nondlib.hpp" file into your project and include it where appropriate. Note that, this library requires at least C++14.

```bash
#Makefile Example
example: myfile.cpp myfile2.cpp myheader2.h maxima.hpp
    g++ -O3 myfile.cpp  myfile2.cpp myheader2.h maxima.hpp -o example

```

## Function Details

The current version of the library provides the following functions:

```cpp
template <typename C, typename M>
void filterQuadD(std::vector<C> &v, M &maxima);
/* Complexity O(n^2)
    Parameters:
    v - Specifies which pointset the algorithm is going to be applied to 
    maxima - Specifies the minimization, maximization behaviour
*/

template <typename C, typename M>
void filterDimSweep2D(std::vector<C> &v, M const &maxima);
/* Complexity O(nlog(n))
    Parameters:
    v - Specifies which pointset the algorithm is going to be applied to 
    maxima - Specifies the minimization, maximization behaviour
*/

template <typename C, typename M>
void filterDimSweep3D(std::vector<C> &v, M const &maxima, size_t obj);
/* Complexity O(nlog(n))
    Parameters:
    v - Specifies which pointset the algorithm is going to be applied to 
    maxima - Specifies the minimization, maximization behaviour
    obj - Specifies which parameter is gonna be used for the sort, Example: 0, would sort in x, and use y and z to solve a 2d problem to find the nonDominated Points
*/

template <typename C, typename M>
void filterDivConqDC(std::vector<C> &v, M const &maxima, int const &dims);
/* Complexity O(log^k-1).n) for k >= 3 
    Parameters:
    v - Specifies which pointset the algorithm is going to be applied to 
    maxima - Specifies the minimization, maximization behaviour
    dims - Number of dimensions of v
*/

template <typename C, typename M>
bool updateMaximaND(std::vector<C> &v, M &maxima, int obj, C &p) {;
/* Complexity O(nd)
    Parameters:
    v - Specifies which pointset the algorithm is going to be applied to 
    maxima - Specifies the minimization, maximization behaviour
    obj - Specifies which parameter is gonna be used for the sort, Example: 0, would sort in x, and use y and z to solve a 2d problem to find the nonDominated Points
    point - Point to update the set with\
    Notes:
        This algorithm is not stable
*/
```

## Example

```cpp
std::vector<std::vector<double>> example = {{1,1}, {2,4}, {1.5,1.6}, {3,2}, {3,3}};
std::vector<double> mx = {1,1,1};
nondlib::inplace::filterDimSweep2D<std::vector<double>,std::vector<double>>(example,mx);

//example is now {{2,4}, {3,3}};
```