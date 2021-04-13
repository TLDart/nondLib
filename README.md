# nondLib: A C++ library for archiving operations of nondominated points

[GitHub](https://github.com/TLDart/nondLib)

**Maintainer:** [Duarte M. Dias](duartedias@student.dei.uc.pt)

**Contributors:**
    [Duarte M. Dias](duartedias@student.dei.uc.pt),
    [Alexandre D. Jesus](https://adbjesus.com).
    [Luís Paquete](https://www.uc.pt/go/paquete/).

----------------------------​

This library implements several operations to filter and maintain a set of
nondominated points. In particular, the current version supports the following
operations:

- `filter` - filter a set of nondominated points from a larger set of points
- `update` - try to insert a new point into a nondominated point set

Both in-place and not-in-place versions exist for all currently implemented
functions, such that the former implies altering the container that is passed to
the function, and the latter simply returns a new container.

The details behind some of the functions are explained in the following article:

Duarte M. Dias, Alexandre D. Jesus, Luís Paquete, A software library for
archiving nondominated points, 2021 (submitted).

**Note:** this is an early version of the library, and as such there may
be breaking changes.

## Usage

The library is provided as a single header for ease of portability. As
such, all there is to do is copy the `nondlib.hpp` file into your
project and include it in your `.cpp` files. Note that the use of this
library currently requires at least C++ 14.

For a usage example, see the [examples](/examples) folder.

## Library Functions Details

Currently, most functions require an archive to be passed as a reference
to an `std::vector<C>`, such that type `C` denotes a point in the
objective space implemented as a container. Exactly which containers are
supported for `C` depends on the particular library function but random
access containers are always valid, e.g `std::vector`, `std::array`, and
`std::deque`.

Moreover, most functions also require a random access container
parameter named `maxima`, that denotes whether minimization (`maxima[i]
= -1`) or maximization (`maxima[i] = 1`) should be considered for a
given objective `i`.

### In-place functions

The following in-place functions are implemented

#### Quadratic-time algorithm for filter operation and any number of dimensions

```cpp
// Parameters:
//   v - Objective point set to be filtered
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - abs(maxima[i]) != 1
template <typename C, typename M>
void nondlib::inplace::filterQuadD(std::vector<C> &v, M const& maxima);
```

#### Dimension-sweep algorithm for filter operation in 2D

```cpp
// Parameters:
//   v - Objective point set to be filtered
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - abs(maxima[i]) != 1
template <typename C, typename M>
void nondlib::inplace::filterDimSweep2D(std::vector<C> &v, M const& maxima);
```

#### Dimension-sweep algorithm for filter operation in 3D

```cpp
// Parameters:
//   v - Objective point set to be filtered
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//   obj - Specifies which objective is gonna be used for sorting, defaults to
//         0, which would sort the points by objective 0, and use objectives 1
//         and 2 to solve a 2d problem and find the nondomianted points
//
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - abs(maxima[i]) != 1
//   - obj > 2
template <typename C, typename M>
void nondlib::inplace::filterDimSweep3D(std::vector<C> &v, M const& maxima, size_t obj = 0); // TODO default obj = 0
```

#### Multidimensional divide-and-conquer algorithm for filter operation in any number of dimensions

```cpp
// Parameters:
//   v - Objective point set to be filtered
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - abs(maxima[i]) != 1
template <typename C, typename M>
void nondlib::inplace::filterDivConqDC(std::vector<C> &v, M const& maxima); // TODO remove dims
```

#### Dimension-sweep algorithm for update operation and any number of dimensions 

```cpp
// Parameters:
//   v - Nondominated point set to be updated
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//   p - point to be inserted into set v (if not dominated) with a type that is
//       random access container and that can be used to construct a type C
//       container
//
// Returns: a boolean denoting whether or not the nondominated set v was updated
//       
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - p.size() != maxima.size()
//   - abs(maxima[i]) != 1
template <typename C, typename M, typename P>
bool nondlib::inplace::updateMaximaND(std::vector<C> &v, M const& maxima, P &&p);
    // TODO static_assert(std::is_same_v<std::remove_cvref_t<C>, std::remove_cvref_t<P>>);
```

### Not-in-place functions

The not-in-place versions of the library are implemented under the
`nondlib::notinplace` namespace using the same function names and
parameters. The main difference is that point set `v` is passed as a
const reference, and that the filtered/updated set is returned from the
function. Example for the `filterQuadD` function:

```cpp
// Parameters:
//   v - Objective point set to be filtered
//   maxima - A random access container denoting whether minimization or
//            maximization should be considered for each objective
//
// Returns: The filtered nondominated set. 
//
// Undefined behavior if:
//   - v[i].size() != maxima.size() for any i
//   - abs(maxima[i]) != 1
template <typename C, typename M>
std::vector<C> nondlib::inplace::filterQuadD(std::vector<C> &v, M const& maxima);
```
