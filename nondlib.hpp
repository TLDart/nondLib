#ifndef __NONDLIB_HPP_
#define __NONDLIB_HPP_

#include <algorithm>
#include <functional>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace nondlib {

// "Private" namespace with helper functions and classes. Warning: may break
// between minor or even patch versions since it is not expected to be used by a
// user.
namespace priv {

template <typename T>
void maximadc_filter2(std::vector<T> &U, std::vector<T> const &V) {
  auto r1 = std::numeric_limits<typename T::value_type>::min();
  auto end = U.begin();
  for (auto itu = U.begin(), itv = V.begin();;) {
    if (itu == U.end()) {
      break;
    }
    if (itv == V.end()) {
      for (; itu != U.end(); ++itu) {
        if ((*itu)[1] >= r1) {
          if (itu != end)
            *end = std::move(*itu);
          ++end;
        }
      }
      break;
    }
    if ((*itv)[0] >= (*itu)[0]) {
      r1 = std::max(r1, (*itv)[1]);
      ++itv;
    } else {
      if ((*itu)[1] >= r1) {
        if (itu != end)
          *end = std::move(*itu);
        ++end;
      }
      ++itu;
    }
  }
  U.erase(end, U.end());
}

template <typename T>
void maximadc_filter3(std::vector<T> &U, std::vector<T> const &V) {
  using D = typename T::value_type;
  std::set<std::array<D, 2>> aux;

  auto dominated = [&aux](auto const &p) {
    auto tmp = std::array<D, 2>{p[1], p[2]};
    auto it = aux.lower_bound(tmp);
    if (it == aux.end()) {
      return false;
    }
    return tmp[1] <= (*it)[1];
  };

  auto insert = [&aux](auto const &p) {
    auto tmp = std::array<D, 2>{p[1], p[2]};
    auto it = aux.lower_bound(tmp);
    if (it != aux.end() && tmp[1] <= (*it)[1]) {  // dominated
      return;
    }
    if (it != aux.end() && tmp[0] == (*it)[0]) {
      it = aux.erase(it);
    }
    while (it != aux.begin()) {
      it = std::prev(it);
      if (tmp[1] >= (*it)[1]) {
        it = aux.erase(it);
      } else {
        it = std::next(it);
        break;
      }
    }
    aux.insert(it, tmp);
  };

  auto end = U.begin();
  for (auto itu = U.begin(), itv = V.begin();;) {
    if (itu == U.end()) {
      break;
    }
    if (itv == V.end()) {
      for (; itu != U.end(); ++itu) {
        if (!dominated(*itu)) {
          if (itu != end)
            *end = std::move(*itu);
          ++end;
        }
      }
      break;
    }
    if ((*itv)[0] >= (*itu)[0]) {
      insert((*itv));
      ++itv;
    } else {
      if (!dominated((*itu))) {
        if (itu != end)
          *end = std::move(*itu);
        ++end;
      }
      ++itu;
    }
  }
  U.erase(end, U.end());
}

// Filters all points in U that are not dominated by V
template <typename T>
void maximadc_filterk(std::vector<T> &U, std::vector<T> &V, size_t k, size_t base) {
  if (V.empty() || U.empty()) {
    return;
  }

  if (V.size() == 1 || U.size() == 1) {  // Filter naively
    auto end = std::remove_if(U.begin(), U.end(), [&V, &k](auto const &u) {
      for (auto const &v : V) {
        size_t i = 0;
        for (i = 0; i < k; ++i) {
          if (v[i] < u[i]) {
            break;
          }
        }
        if (i == k) {
          return true;
        }
      }
      return false;
    });
    U.erase(end, U.end());
    return;
  }

  if (k == base) {
    if (base == 2) {
      maximadc_filter2(U, V);
    } else {
      maximadc_filter3(U, V);
    }
    return;
  }

  // Equipartition V
  std::vector<size_t> inds;
  inds.reserve(V.size());
  for (size_t i = 0; i < V.size(); i++)
    inds.emplace_back(i);

  auto cmp = [&V, &k](auto const &lhs, auto const &rhs) {
    if (V[lhs][k - 1] > V[rhs][k - 1]) {
      return true;
    } else if (V[lhs][k - 1] < V[rhs][k - 1]) {
      return false;
    } else {
      return lhs > rhs;
    }
  };

  auto mid = inds.begin() + inds.size() / 2;
  std::nth_element(inds.begin(), mid, inds.end(), cmp);

  std::vector<bool> isV1(V.size(), false);
  for (auto it = inds.begin(); it != mid; ++it) {
    isV1[*it] = true;
  }

  auto ud = V[*mid][k - 1];

  std::vector<T> V1, V2;
  V1.reserve(V.size() / 2);
  V2.reserve(V.size() / 2 + 1);
  for (size_t i = 0; i < V.size(); ++i) {
    if (isV1[i]) {
      V1.push_back(std::move(V[i]));
    } else {
      V2.push_back(std::move(V[i]));
    }
  }

  // Partition U by V2[0]
  std::vector<T> U1, U2;
  U1.reserve(U.size());
  U2.reserve(U.size());
  for (size_t i = 0; i < U.size(); ++i) {
    if (U[i][k - 1] > ud) {
      U1.push_back(std::move(U[i]));
    } else {
      U2.push_back(std::move(U[i]));
    }
  }

  // Step 2
  maximadc_filterk(U2, V2, k, base);
  maximadc_filterk(U2, V1, k - 1, base);
  maximadc_filterk(U1, V1, k, base);

  // Merge U1 and U2 but in original order
  size_t i = 0;
  for (auto ita = U1.begin(), itb = U2.begin();;) {
    if (ita == U1.end()) {
      for (; itb != U2.end(); ++itb) {
        U[i++] = std::move(*itb);
      }
      break;
    }
    if (itb == U2.end()) {
      for (; ita != U1.end(); ++ita) {
        U[i++] = std::move(*ita);
      }
      break;
    }
    if ((*ita)[0] >= (*itb)[0]) {
      U[i++] = std::move(*ita++);
    } else {
      U[i++] = std::move(*itb++);
    }
  }
  U.erase(U.begin() + i, U.end());

  // Merge V1 and V2 but in original order
  i = 0;
  for (auto ita = V1.begin(), itb = V2.begin();;) {
    if (ita == V1.end()) {
      for (; itb != V2.end(); ++itb) {
        V[i++] = std::move(*itb);
      }
      break;
    }
    if (itb == V2.end()) {
      for (; ita != V1.end(); ++ita) {
        V[i++] = std::move(*ita);
      }
      break;
    }
    if ((*ita)[0] >= (*itb)[0]) {
      V[i++] = std::move(*ita++);
    } else {
      V[i++] = std::move(*itb++);
    }
  }
}

template <typename T>
void maximadc_maximak(std::vector<T> &S, size_t k, size_t base) {
  if (S.size() <= 1) {
    return;
  }

  // Find the median
  std::vector<size_t> inds;
  inds.reserve(S.size());
  for (size_t i = 0; i < S.size(); i++)
    inds.emplace_back(i);

  auto cmp = [&S, &k](auto const &lhs, auto const &rhs) {
    if (S[lhs][k - 1] > S[rhs][k - 1]) {
      return true;
    } else if (S[lhs][k - 1] < S[rhs][k - 1]) {
      return false;
    } else {
      return lhs > rhs;
    }
  };

  auto mid = inds.begin() + inds.size() / 2;
  std::nth_element(inds.begin(), mid, inds.end(), cmp);

  std::vector<bool> isA(S.size(), false);
  for (auto it = inds.begin(); it != mid; ++it) {
    isA[*it] = true;
  }

  std::vector<T> A, B;
  A.reserve(S.size() / 2);
  B.reserve(S.size() / 2 + 1);
  for (size_t i = 0; i < S.size(); ++i) {
    if (isA[i]) {
      A.push_back(std::move(S[i]));
    } else {
      B.push_back(std::move(S[i]));
    }
  }

  // Step 2
  maximadc_maximak(A, k, base);
  maximadc_maximak(B, k, base);
  maximadc_filterk(B, A, k - 1, base);

  // Merge A \cup B but in original order
  size_t i = 0;
  for (auto ita = A.begin(), itb = B.begin();;) {
    if (ita == A.end()) {
      while (itb != B.end()) {
        S[i++] = std::move(*itb++);
      }
      break;
    }
    if (itb == B.end()) {
      while (ita != A.end()) {
        S[i++] = std::move(*ita++);
      }
      break;
    }
    if ((*ita)[0] >= (*itb)[0]) {
      S[i++] = std::move(*ita++);
    } else {
      S[i++] = std::move(*itb++);
    }
  }
  S.erase(S.begin() + i, S.end());
}

template <typename T>
bool cmp2(T &lhs, T &rhs) {
  if (lhs.c[1] < rhs.c[1])
    return true;
  if (lhs.c[1] == rhs.c[1])
    return lhs.c[0] < rhs.c[0];
  return false;
}

template <typename T>
bool cmpLex(T const &lhs, T const &rhs) {
  for (int i = 0; i < lhs.c.size(); i++) {
    if (lhs.c[i] < rhs.c[i])
      return true;
    if (lhs.c[i] > rhs.c[i])
      return false;
  }
  return false;
}

// Sorts in lexicographical order, using the n + 1
// dimension as a tie breaker
struct CompleteColSortReverse {
  template <typename T, typename U>
  bool operator()(T const &lhs, U const &rhs) const {
    for (size_t i = 0; i < lhs.size(); i++) {
      if (lhs[i] < rhs[i])
        return false;
      if (lhs[i] > rhs[i])
        return true;
    }
    return true;
  }
};

class AnyColSortReverse {
  size_t param;

 public:
  explicit AnyColSortReverse(size_t p)
      : param(p) {}

  template <typename T, typename U>
  bool operator()(T const &lhs, U const &rhs) const {
    if (lhs[param] < rhs[param]) {
      return false;
    } else if (lhs[param] == rhs[param]) {
      for (size_t i = 0; i < lhs.size(); i++) {
        if (lhs[i] < rhs[i])
          return false;
      }
    }
    return true;
  }
};

template <typename T, typename M>
void multiplyMaxima(T &v, M const &maxima) {
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v[i][j] *= maxima[j];
    }
  }
}

template <typename T, typename U>
bool dominates(T &p1, U &p2) {
  for (size_t i = 0; i < p1.size(); i++) {
    if (p2[i] > p1[i])
      return false;
  }
  return true;
}

template <typename T>
bool _dominated2d(std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> &aux,
                  std::vector<T> const &p, size_t lo, size_t hi) {
  /*
      This function is an integrating part of the 3d maxima, it used to solve a
     2d problem, uses a auxiliary set<vector<double>> to keep track, of the
     dominated points Beware that the trick here is that we use a reverse array
     in order not to use reverse iterators
  */
  auto it = aux.upper_bound({p[lo], p[hi]});  // Gets the iterator to the element after which the
                                              // element should be placed, the array is guaranteed
                                              // to have atleast 1 element;
  if (it != aux.begin()) {                    // Here we want to test wether the previous elemnt
                                              // dominates current. In case there is no such
                                              // element (aka if the element will be the first in
                                              // the list) we cant do that
    auto prev = it;
    prev--;
    if ((*prev)[1] >= p[hi]) {  // Condition evaluates if the point is already
                                // dominated by some other point
      return true;
    }
  }

  auto start = it;
  while (it != aux.end()) {  // Here we loop until the end and check if there are
                             // dominated elements, then we remove them all at once
    if ((*it)[1] > p[hi]) {  // What happens here is, we use a aux var to store the deleted
                             // iterator, then there is no need to increment rit since one
                             // element is remove and it point to the next position
      break;
    } else  // Because the points are dominated this means that once itY > pY,
            // all it other itY > pY are therefore dont need checking
      it++;
  }
  aux.erase(start, it);
  aux.insert({p[lo], p[hi]});
  return false;  // add the point to the set, meaning that the point is not
                 // dominated
}
}  // namespace priv

namespace inplace {
using namespace nondlib::priv;

template <typename C, typename M>
void filterQuadD(std::vector<C> &v, M const &maxima) {
  /*
  Algorithm for calculating the maxima of a point set. Complexity O(n^2 * d),
  where n is the numbers of points in the set and d is the number of dimensions
  pointset defines a structure that containing all N-Dimensional points.
  maxima defines a d-sized array containing either 1 or -1 that defines if the d
  is to be maximized (1) or minimized (-1) Particular to this function Points
  are removed using a swap-pop scheme
  */
  multiplyMaxima(v, maxima);
  for (size_t point = 0; point < v.size(); point++) {
    for (size_t rempoint = point + 1; rempoint < v.size(); ++rempoint) {
      bool dominated = true, dominator = true;
      for (size_t dim = 0; dim < v[point].size() && (dominated || dominator);
           dim++) {                              // And we verify if it respects the conditions
        if (v[point][dim] > v[rempoint][dim]) {  // If the points contains atleast 1 dimension
                                                 // greater than some other point it is not
                                                 // dominated and cant be removed
          dominated = false;
        }
        if (v[point][dim] < v[rempoint][dim]) {  // Neither Dominated nor dominator
          dominator = false;
        }
      }
      if (dominated) {                                    // Remove the point if it is dominated
        std::swap(*(v.begin() + point), *(v.end() - 1));  // Swap Point with last
        v.pop_back();                                     // Remove last
        point--;               // Decrement the loop because of the removed point
        break;                 // Skip to the next point
      } else if (dominator) {  // Same idea but for the other point
        std::swap(*(v.begin() + rempoint), *(v.end() - 1));
        v.pop_back();
        rempoint--;
      }
    }
  }
  multiplyMaxima(v, maxima);
}

template <typename C, typename M>
void filterDimSweep2D(std::vector<C> &v, M const &maxima) {
  if (v.empty()) {
    return;
  }

  /*Calculates the non dominated points on a 2d set using a matrix of vectors as
  its core data structure Complexity : NlogN
  */
  multiplyMaxima(v, maxima);

  std::sort(v.begin(), v.end(),
            CompleteColSortReverse());  // Sort in x in reverse order
  double maxy = v[0][1];                // Y coord of the last point in x
  size_t idx = 1;
  for (size_t point = 1; point < v.size(); point++) {
    if (v[point][1] > maxy) {
      if (idx != point) {
        v[idx][0] = v[point][0];
        v[idx][1] = v[point][1];
      }
      maxy = v[point][1];
      idx++;
    }
  }
  v.erase(v.begin() + idx,
          v.end());  // May remove this, but taking in consideration it is in the
                     // end, this operation should be done in constant time
  multiplyMaxima(v, maxima);
}

/*
 * Note:
 *   - Undefined behavior if v[i].size() != maxima.size() for all i
 *   - Undefined behavior if v.size() == 0
 */
template <typename C, typename M>
void filterDimSweep3D(std::vector<C> &v, M const &maxima, size_t obj = 0) {
  if (v.empty()) {
    return;
  }

  // std::cout << v.size() << std::endl;
  multiplyMaxima(v, maxima);

  std::sort(v.begin(), v.end(), AnyColSortReverse(obj));
  size_t hi, lo;
  hi = (obj + 1) % 3 > (obj + 2) % 3 ? (obj + 1) % 3 : (obj + 2) % 3;
  lo = 3 - obj - hi;
  size_t first = 1, point;

  using T = typename C::value_type;
  std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> tmp = {{v[0][lo], v[0][hi]}};

  for (point = 1; point < v.size(); point++) {
    if (!(_dominated2d<T>(tmp, v[point], lo, hi))) {
      if (first != point) {
        v[first][0] = v[point][0];
        v[first][1] = v[point][1];
        v[first][2] = v[point][2];
      }
      first++;
    }
  }
  v.erase(v.begin() + first, v.end());

  multiplyMaxima(v, maxima);
}

template <typename C, typename M>
void filterDivConqD(std::vector<C> &v, M const &maxima, size_t base = 3) {
  if (base != 2 && base != 3) {
    throw("Base case for divide and conquer must be 2 or 3!");
  }

  size_t const dims = maxima.size();
  if (dims < 2) {
    throw("There must be at least 2 objectives!");
  } else if (dims == 2) {
    return filterDimSweep2D(v, maxima);
  } else if (dims == 3) {
    return filterDimSweep3D(v, maxima);
  }

  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v[i][j] *= maxima[j];
    }
  }

  auto lexsort = [](auto const &lhs, auto const &rhs) {
    for (size_t i = 0; i < lhs.size(); i++) {
      if (lhs[i] > rhs[i]) {
        return true;
      }
      if (lhs[i] < rhs[i]) {
        return false;
      }
    }
    return true;
  };

  sort(v.begin(), v.end(), lexsort);
  maximadc_maximak(v, dims, base);

  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v[i][j] *= maxima[j];
    }
  }
}

// Retuns true if point was inserted in v, i.e. if it is non-dominated, or false
// otherwise.
template <typename C, typename M, typename P>
bool updateMaximaND(std::vector<C> &v, M const &maxima, P &&p) {
  auto dominates = [&maxima](auto const &lhs, auto const &rhs) {
    for (size_t i = 0; i < maxima.size(); ++i) {
      if (lhs[i] * maxima[i] < rhs[i] * maxima[i]) {
        return false;
      }
    }
    return true;
  };

  size_t i = 0;
  for (; i < v.size(); ++i) {
    if (dominates(v[i], p)) {
      return false;  // return if v[i] dominates p
    }
    if (dominates(p, v[i])) {
      break;
    }
  }
  size_t end = i++;
  for (; i < v.size(); ++i) {
    if (!dominates(p, v[i])) {
      v[end++] = std::move(v[i]);
    }
  }
  v.erase(v.begin() + end, v.end());
  v.emplace_back(std::forward<P>(p));
  return true;
}
}  // namespace inplace

namespace notinplace {
using namespace nondlib::priv;

template <typename C, typename M>
std::vector<C> filterQuadD(std::vector<C> const &v, M const &maxima) {
  auto res = std::vector<C>();
  res.reserve(v.size());
  for (auto it = v.begin(); it != v.end(); ++it) {
    nondlib::inplace::updateMaximaND(res, maxima, *it);
  }
  return res;
}

template <typename C, typename M>
std::vector<C> filterDimSweep3D(std::vector<C> const &v, M const &maxima, size_t obj = 0) {
  if (v.empty()) {
    return {};
  }
  std::vector<C> pointset;
  pointset.reserve(v.size());

  auto aux = v;
  multiplyMaxima(aux, maxima);
  std::sort(aux.begin(), aux.end(), AnyColSortReverse(obj));
  size_t hi, lo;
  hi = (obj + 1) % 3 > (obj + 2) % 3 ? (obj + 1) % 3 : (obj + 2) % 3;
  lo = 3 - obj - hi;

  using T = typename C::value_type;
  std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> tmp = {{aux[0][lo], aux[0][hi]}};
  pointset.push_back(std::move(aux[0]));
  for (size_t i = 1; i < aux.size(); i++) {
    if (!(_dominated2d(tmp, aux[i], lo, hi))) {
      pointset.push_back(std::move(aux[i]));
    }
  }
  multiplyMaxima(pointset, maxima);
  return pointset;
}

template <typename C, typename M>
std::vector<C> filterDimSweep2D(std::vector<C> const &v, M const &maxima) {
  std::vector<C> res = v;
  nondlib::inplace::filterDimSweep2D(res, maxima);
  return res;
}

template <typename C, typename M>
std::vector<C> filterDivConqD(std::vector<C> const &v, M const &maxima, size_t base = 3) {
  const size_t dims = maxima.size();
  if (dims < 2) {
    throw("There must be at least 2 objectives!");
  } else if (dims == 2) {
    return filterDimSweep2D(v, maxima);
  } else if (dims == 3) {
    return filterDimSweep3D(v, maxima);
  }

  std::vector<C> res = v;
  nondlib::inplace::filterDivConqD(res, maxima, base);
  return res;
}

// Retuns true if point was inserted in v, i.e. if it is non-dominated, or false
// otherwise.
template <typename C, typename M, typename P>
std::vector<C> updateMaximaND(std::vector<C> const &v, M const &maxima, P &&p) {
  auto dominates = [&maxima](auto const &lhs, auto const &rhs) {
    for (size_t i = 0; i < maxima.size(); ++i) {
      if (lhs[i] * maxima[i] < rhs[i] * maxima[i]) {
        return false;
      }
    }
    return true;
  };

  auto it = v.begin();
  for (; it != v.end(); ++it) {
    if (dominates(*it, p))
      return v;
    if (dominates(p, *it))
      break;
  }

  std::vector<C> res;
  res.reserve(v.size());
  std::copy(v.begin(), it, std::back_inserter(res));
  for (it = it != v.end() ? ++it : it; it != v.end(); ++it) {
    if (!dominates(p, *it)) {
      res.push_back(*it);
    }
  }
  res.emplace_back(std::forward<P>(p));
  return res;
}

}  // namespace notinplace
}  // namespace nondlib

#endif
