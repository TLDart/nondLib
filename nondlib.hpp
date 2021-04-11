#ifndef MAXIMA
#define MAXIMA
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <list>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nondlib {

template <typename T>
class Point {
   public:
    T c;
    int rank;
    bool isA;
    bool isB;

    Point(T c, int rank = 0, bool isA = false, bool isB = false) : c(c), rank(0), isA(isA), isB(isB) {}

    bool operator==(const Point &other) const { return this->c == other.c; };

    bool operator<(const Point &other) const { return this->c < other.c; };

    std::string printdims() const {
        std::string s;
        s = std::to_string(this->c[0]);
        for (int i = 1; i < c.size(); i++) {
            s += "," + std::to_string(this->c[i]);
        }
        return s;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &p) {
        return os << "Point(" << p.printdims() << "," << p.rank << ")"
                  << "rank " << p.rank << std::endl;
    }
};

template <typename T>
class cmpN {
    int param;

   public:
    cmpN(size_t p) : param(p) {}

    bool operator()(T const &lhs, T const &rhs) {
        if (lhs.c.get()[param] < rhs.c.get()[param])
            return true;
        else if (lhs.c.get()[param] == rhs.c.get()[param]) {
            for (size_t i = 0; i < lhs.c.get().size(); ++i) {
                if (lhs.c.get()[i] < rhs.c.get()[i]) {
                    return true;
                }
                if (lhs.c.get()[i] > rhs.c.get()[i]) {
                    return false;
                }
            }
        }
        return false;
    }
};

template <typename T>
bool cmp2(T &lhs, T &rhs) {
    if (lhs.c[1] < rhs.c[1]) return true;
    if (lhs.c[1] == rhs.c[1]) return lhs.c[0] < rhs.c[0];
    return false;
}

template <typename T>
bool cmpLex(T const &lhs, T const &rhs) {
    for (int i = 0; i < lhs.c.size(); i++) {
        if (lhs.c[i] < rhs.c[i]) return true;
        if (lhs.c[i] > rhs.c[i]) return false;
    }
    return false;
}

struct CompleteColSortReverse  // Sorts in lexicographical order, using the n + 1 dimension as a tie breaker
{
    bool operator()(const std::vector<double> &lhs, const std::vector<double> &rhs) const {
        for (int i = 0; i < lhs.size(); i++) {
            if (lhs[i] < rhs[i]) return false;
            if (lhs[i] > rhs[i]) return true;
        }
        return true;
    }
};

class AnyColSortReverse {
   public:
    size_t param;
    AnyColSortReverse(size_t p) : param(p) {}

    bool operator()(const std::vector<double> &lhs, const std::vector<double> &rhs) const {
        if (lhs[param] < rhs[param]) {
            return false;
        } else if (lhs[param] == rhs[param]) {
            for (size_t i = 0; i < lhs.size(); i++) {
                if (lhs[i] < rhs[i]) return false;
            }
        }
        return true;
    }
};

template <typename T, typename M>
void multiplyMaxima(T v, M maxima) {
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[0].size(); ++j) {  // Initial multiplication required if because of maximization/ minimization criteria
            v[i][j] *= maxima[j];
        }
    }
}
template <typename T>
bool compareReverse(const std::vector<double> &lhs, const std::vector<double> &rhs, int param) {
    if (lhs[param] < rhs[param]) {
        return false;
    } else if (lhs[param] == rhs[param]) {
        for (size_t i = 0; i < lhs.size(); i++) {
            if (lhs[i] < rhs[i]) return false;
        }
    }
    return true;
}
template <typename C>
bool _dominates(C &p1, C &p2) {
    for (size_t i = 0; i < p1.size(); i++) {
        if (p2[i] > p1[i]) return false;
    }
    return true;
}
// Retuns true if point was inserted in v, i.e. if it is non-dominated, or false otherwise.
template <typename C, typename M>
bool updateMaximaND(std::vector<C> &v, M &maxima, int obj, C &p) {
    size_t i = 0;
    for (; i < v.size(); ++i) {
        if (_dominates(v[i], p)) return false;  // return if v[i] dominates p
        if (_dominates(p, v[i])) break;
    }
    size_t end = i++;
    for (; i < v.size(); ++i) {
        bool ok = !_dominates(p, v[i]);
        if (ok) v[end++] = std::move(v[i]);
    }
    v.erase(v.begin() + end, v.end());
    v.push_back(p);
    return true;
}

template <typename C, typename M>
void updateMaxima2D(std::vector<C> &v, M &maxima, C &point) {
    // Complexity, Search in O(logN), insert in O(n)
    // Vector is reverse sorted

    size_t idx, last;
    size_t nxt, lst;

    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);
    point[0] *= maxima[0];
    point[1] *= maxima[1];

    // Update to binary search
    for (idx = 0; idx < v.size(); idx++) {  // get the point location
        if (v[idx][0] <= point[0]) {
            break;  // find the index of the point to remove
        }
    }

    nxt = idx + 1;
    lst = idx + 1;

    if (idx == 0 || v[idx - 1][1] < point[1]) {
        v.insert(v.begin() + idx, point);           // Insert the point
        for (int i = nxt; i < v.size() - 1; i++) {  // Find all elements are to be removed with this set update
            if (v[i][1] > point[0]) {
                lst = i;
                break;
            }
        }
        if (lst != nxt) {
            for (size_t i = lst; i < v.size() - 1; i++) {  // Copy the elements using 2 pointers to the correct location
                v[nxt][0] = v[lst][0];
                v[nxt][1] = v[lst][1];
                nxt++;
                break;
            }
        }
        v.erase(v.begin() + nxt, v.end());  // Remove the elements from the end of the vector
    }
    // else do nothing, the point is dominated
    multiplyMaxima<std::vector<C>, M>(v, maxima);
}

template <typename C, typename M>
void updateMaxima3D(std::vector<C> &v, M &maxima, int obj, C &point) {
    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);

    point[0] *= maxima[0];
    point[1] *= maxima[1];
    point[2] *= maxima[2];

    size_t hi, lo, pos, first;
    bool valid = false;
    hi = (obj + 1) % 3 > (obj + 2) % 3 ? (obj + 1) % 3 : (obj + 2) % 3;
    lo = 3 - obj - hi;

    using T = typename C::value_type;
    std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> aux = {};

    for (pos = 0; pos < (int)v.size(); pos++) {  // Finding the insert Pos
        if (compareReverse<int>(point, v[pos], obj)) {
            break;
        }
        aux.insert({v[pos][lo], v[pos][hi]});  // No need to check insertion since all points in the all are already confirmed non dominated
    }
    if (pos == 0) {
        aux.insert({point[lo], point[hi]});
        v.insert(v.begin() + pos, point);
        first = 1;
        valid = true;
    } else if (!(_dominated2d(aux, point, lo, hi))) {  // If the point is not dominated by the current points
        v.insert(v.begin() + pos, point);
        first = pos + 1;
        valid = true;
    }

    if (valid) {
        for (size_t i = first; i < v.size(); i++) {
            if (!(_dominated2d(aux, v[i], lo, hi))) {
                if (first != i) {
                    v[first][0] = v[i][0];
                    v[first][1] = v[i][1];
                    v[first][2] = v[i][2];
                }
                first++;
            }
        }
        v.erase(v.begin() + first, v.end());  // Erase a range of points from the end of the vector.
    }

    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);
}

template <typename T>
void ecdf2(std::vector<std::reference_wrapper<Point<T>>> &S) {
    if (S.size() == 0) return;
    auto last = std::numeric_limits<typename T::type::value_type>::max();
    for (size_t i = 0; i < S.size(); ++i) {
        if (S[i].get().c.get()[1] >= last) {
            S[i].get().rank = 1;
        } else {
            last = std::min(last, S[i].get().c.get()[1]);
        }
    }
    return;
}
template <typename T>
void ecdf2_modified(std::vector<std::reference_wrapper<Point<T>>> &S) {
    if (S.size() == 0) return;
    auto last = std::numeric_limits<typename T::type::value_type>::max();
    for (size_t i = 0; i < S.size(); ++i) {
        if (S[i].get().isB && S[i].get().c.get()[1] >= last) {
            S[i].get().rank = 1;
        } else if (S[i].get().isA) {
            last = std::min(last, S[i].get().c.get()[1]);
        }
    }
    return;
}

template <typename T>
void ecdfk_modified(std::vector<std::reference_wrapper<Point<T>>> &S, int k) {
    // Base case
    if (S.size() == 1) {
        S[0].get().rank = 0;
        return;
    }
    if (S.size() == 0) {
        return;
    }
    if (k == 2) {
        ecdf2_modified(S);
        return;
    }

    // Find the median
    std::vector<size_t> sA;
    sA.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        sA.emplace_back(i);
    }

    auto cmp = cmpN<Point<T>>(k - 1);
    auto cmp2 = [&cmp, &S](auto const &lhs, auto const &rhs) { return cmp(S[lhs], S[rhs]); };
    auto mid = sA.begin() + sA.size() / 2;
    std::nth_element(sA.begin(), mid, sA.end(), cmp2);

    std::vector<bool> isA(S.size(), false);
    for (auto it = sA.begin(); it != mid; ++it) {
        isA[*it] = true;
    }

    std::vector<std::reference_wrapper<Point<T>>> A, B;
    A.reserve(S.size());
    B.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        if (isA[i]) {
            A.push_back(S[i]);
        } else {
            B.push_back(S[i]);
        }
    }

    // Step 2
    ecdfk_modified(A, k);
    ecdfk_modified(B, k);

    // Build aux
    std::vector<std::pair<Point<T>, std::reference_wrapper<Point<T>>>> aux;
    aux.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        if (S[i].get().rank > 0) continue;
        if (isA[i]) {
            aux.emplace_back(Point<T>(S[i].get().c, 0, S[i].get().isA, false), S[i]);
        } else {
            aux.emplace_back(Point<T>(S[i].get().c, 0, false, S[i].get().isB), S[i]);
        }
    }

    // Merge
    std::vector<std::reference_wrapper<Point<T>>> v;
    v.reserve(aux.size());
    for (auto &p : aux) {
        v.push_back(p.first);
    };
    ecdfk_modified(v, k - 1);

    for (auto &p : aux) {
        if (p.first.isB) {
            p.second.get().rank += p.first.rank;
        }
    }
}

template <typename T>
void ecdfk(std::vector<std::reference_wrapper<Point<T>>> &S, int k) {
    // Base case
    if (S.size() == 1) {
        S[0].get().rank = 0;
        return;
    }
    if (S.size() == 0) {
        return;
    }
    if (k == 2) {
        ecdf2(S);
        return;
    }

    // Find the median
    std::vector<size_t> sA;
    sA.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        sA.emplace_back(i);
    }

    auto cmp = cmpN<Point<T>>(k - 1);
    auto cmp2 = [&cmp, &S](auto const &lhs, auto const &rhs) { return cmp(S[lhs], S[rhs]); };
    auto mid = sA.begin() + sA.size() / 2;
    std::nth_element(sA.begin(), mid, sA.end(), cmp2);

    std::vector<bool> isA(S.size(), false);
    for (auto it = sA.begin(); it != mid; ++it) {
        isA[*it] = true;
    }

    std::vector<std::reference_wrapper<Point<T>>> A;
    A.reserve(S.size());
    std::vector<std::reference_wrapper<Point<T>>> B;
    B.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        if (isA[i]) {
            A.push_back(S[i]);
        } else {
            B.push_back(S[i]);
        }
    }

    // Step 2
    ecdfk(A, k);
    ecdfk(B, k);

    // Build aux
    std::vector<std::pair<Point<T>, std::reference_wrapper<Point<T>>>> aux;
    aux.reserve(S.size());
    for (size_t i = 0; i < S.size(); i++) {
        if (S[i].get().rank != 0) continue;
        if (isA[i]) {
            aux.emplace_back(Point<T>(S[i].get().c, 0, true, false), S[i]);
        } else {
            aux.emplace_back(Point<T>(S[i].get().c, 0, false, true), S[i]);
        }
    }

    // Merge
    std::vector<std::reference_wrapper<Point<T>>> v;
    v.reserve(aux.size());
    for (auto &p : aux) {
        v.push_back(p.first);
    };
    ecdfk_modified(v, k - 1);

    for (auto &p : aux) {
        if (p.first.isB) {
            p.second.get().rank += p.first.rank;
        }
    }
}

template <typename T>
std::vector<Point<T>> setMaximaK(std::vector<Point<T>> &S, int dims) {
    for (auto &p : S) {
        for (int k = 0; k < dims; k++) {
            p.c[k] = -p.c[k];
        }
    }

    // printPointsP(S, "prev sort");
    sort(S.begin(), S.end(), cmpLex<Point<T>>);
    std::vector<Point<std::reference_wrapper<T>>> aux;
    aux.reserve(S.size());
    for (auto &s : S) {
        aux.emplace_back(std::ref(s.c));
    }
    std::vector<std::reference_wrapper<Point<std::reference_wrapper<T>>>> aux2(aux.begin(), aux.end());
    ecdfk(aux2, dims);

    for (auto &p : S) {
        for (int k = 0; k < dims; k++) {
            p.c[k] = -p.c[k];
        }
    }

    std::vector<Point<T>> nondom;
    for (int i = 0; i < aux.size(); i++) {
        if (aux[i].rank == 0) {
            nondom.push_back(S[i]);
        }
    }
    return nondom;
}

template <typename T>
bool _dominated2d(std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> &aux, std::vector<double> p, size_t lo, size_t hi) {
    /*
        This function is an integrating part of the 3d maxima, it used to solve a 2d problem, uses a auxiliary set<vector<double>> to keep track, of the dominated points
        Beware that the trick here is that we use a reverse array in order not to use reverse iterators
    */
    auto it = aux.upper_bound({p[lo], p[hi]});  // Gets the iterator to the element after which the element should be placed, the array is guaranteed to have atleast 1 element;
    if (it != aux.begin()) {                    // Here we want to test wether the previous elemnt dominates current. In case there is no such element (aka if the element will be the first in the list) we cant do that
        auto prev = it;
        prev--;
        if ((*prev)[1] >= p[hi]) {  // Condition evaluates if the point is already dominated by some other point
            return true;
        }
    }

    auto start = it;
    while (it != aux.end()) {    // Here we loop until the end and check if there are dominated elements, then we remove them all at once
        if ((*it)[1] > p[hi]) {  // What happens here is, we use a aux var to store the deleted iterator, then there is no need to increment rit since one element is remove and it point to the next position
            break;
        } else  // Because the points are dominated this means that once itY > pY, all it other itY > pY are therefore dont need checking
            it++;
    }
    aux.erase(start, it);
    aux.insert({p[lo], p[hi]});
    return false;  // add the point to the set, meaning that the point is not dominated
}

namespace inplace {
template <typename C, typename M>
void filterQuadD(std::vector<C> &v, M &maxima) {
    /*
    Algorithm for calculating the maxima of a point set. Complexity O(n^2 * d), where n is the numbers of points in the set and d is the number of dimensions
    pointset defines a structure that containing all N-Dimensional points.
    maxima defines a d-sized array containing either 1 or -1 that defines if the d is to be maximized (1) or minimized (-1)
    Particular to this function Points are removed using a swap-pop scheme
    */
    multiplyMaxima<std::vector<C>, M>(v, maxima)

        std::vector<std::vector<double>>
            result;
    for (int point = 0; point <= v.size() - 1; point++) {
        for (int rempoint = point + 1; rempoint != v.size(); rempoint++) {  // Since we do a comparison on each side there is no need to test the points before the point
            bool dominated = true, dominator = true;
            for (int dim = 0; dim != v[point].size() && (dominated || dominator); dim++) {  // And we verify if it respects the conditions
                if (v[point][dim] > v[rempoint][dim]) {                                     // If the points contains atleast 1 dimension greater than some other point it is not dominated and cant be removed
                    dominated = false;
                }
                if (v[point][dim] < v[rempoint][dim]) {  // Neither Dominated nor dominator
                    dominator = false;
                }
            }
            if (dominated) {                                      // Remove the point if it is dominated
                std::swap(*(v.begin() + point), *(v.end() - 1));  // Swap Point with last
                v.pop_back();                                     // Remove last
                point--;                                          // Decrement the loop because of the removed point
                break;                                            // Skip to the next point
            } else if (dominator) {                               // Same idea but for the other point
                std::swap(*(v.begin() + rempoint), *(v.end() - 1));
                v.pop_back();
                rempoint--;
            }
        }
    }
    multiplyMaxima<std::vector<C>, M>(v, maxima)
}

template <typename C, typename M>
void filterDimSweep2D(std::vector<C> &v, M const &maxima) {
    /*Calculates the non dominated points on a 2d set using a matrix of vectors as its core data structure
    Complexity : NlogN
    */
    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);

    std::sort(v.begin(), v.end(), CompleteColSortReverse());  // Sort in x in reverse order
    double maxy = v[0][1];                                    // Y coord of the last point in x
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
    v.erase(v.begin() + idx, v.end());  // May remove this, but taking in consideration it is in the end, this operation should be done in constant time
    multiplyMaxima<std::vector<C>, M>(v, maxima);
}

/*
 * Note:
 *   - Undefined behavior if v[i].size() != maxima.size() for all i
 *   - Undefined behavior if v.size() == 0
 */

template <typename C, typename M>
void filterDimSweep3D(std::vector<C> &v, M const &maxima, size_t obj) {
    // std::cout << v.size() << std::endl;
    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);

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

    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);
}

template <typename C, typename M>
void filterDivConqDCC(std::vector<C> &v, M const &maxima, int const &dims) {
    using RefC = std::reference_wrapper<C>;
    using PointRefC = Point<RefC>;
    using RefPointRefC = std::reference_wrapper<PointRefC>;
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[0].size(); ++j) {
            v[i][j] *= -maxima[j];
        }
    }
    auto lexsort = [](auto const &lhs, auto const &rhs) {
        for (size_t i = 0; i < lhs.size(); i++) {
            if (lhs[i] > rhs[i]) return false;
            if (lhs[i] < rhs[i]) return true;
        }
        return true;
    };
    sort(v.begin(), v.end(), lexsort);
    std::vector<PointRefC> a(v.begin(), v.end());
    std::vector<RefPointRefC> b(a.begin(), a.end());
    ecdfk(b, dims);
    size_t end = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].rank == 0) {
            v[end++] = std::move(v[i]);
        }
    }
    v.erase(v.begin() + end, v.end());
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[0].size(); ++j) {
            v[i][j] *= -maxima[j];
        }
    }
}

}  // namespace inplace

namespace notinplace {
template <typename C, typename M>
std::vector<C> filterDivConqDC(std::vector<C> &v, M const &maxima, int const &dims) {
    using RefC = std::reference_wrapper<C>;
    using PointRefC = Point<RefC>;
    using RefPointRefC = std::reference_wrapper<PointRefC>;
    std::vector<C> pointset;
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[0].size(); ++j) {
            v[i][j] *= -maxima[j];
        }
    }
    auto lexsort = [](auto const &lhs, auto const &rhs) {
        for (size_t i = 0; i < lhs.size(); i++) {
            if (lhs[i] > rhs[i]) return false;
            if (lhs[i] < rhs[i]) return true;
        }
        return true;
    };
    sort(v.begin(), v.end(), lexsort);
    std::vector<PointRefC> a(v.begin(), v.end());
    std::vector<RefPointRefC> b(a.begin(), a.end());
    ecdfk(b, dims);
    size_t end = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].rank == 0) {
            pointset.push_back(a);
        }
    }
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[0].size(); ++j) {
            v[i][j] *= -maxima[j];
        }
    }
    for (size_t i = 0; i < pointset.size(); ++i) {
        for (size_t j = 0; j < poinset[0].size(); ++j) {
            pointset[i][j] *= -maxima[j];
        }
    }
    return pointset;
}

template <typename C, typename M>
std::vector<C> filterDimSweep3D(std::vector<C> &v, M const &maxima, size_t obj) {
    std::vector<C> pointset;
    // std::cout << v.size() << std::endl;
    multiplyMaxima<std::vector<C> &, M const &>(v, maxima);

    std::sort(v.begin(), v.end(), AnyColSortReverse(obj));
    size_t hi, lo;
    hi = (obj + 1) % 3 > (obj + 2) % 3 ? (obj + 1) % 3 : (obj + 2) % 3;
    lo = 3 - obj - hi;
    size_t first = 1, point;

    using T = typename C::value_type;
    std::set<std::array<T, 2>, std::greater<std::array<T, 2>>> tmp = {{v[0][lo], v[0][hi]}};

    for (point = 1; point < v.size(); point++) {
        if (!(_dominated2d<T, C>(tmp, v[point], lo, hi))) {
            pointset.push_back(v[point]);
        }
    }
    multiplyMaxima(pointset, maxima);
    multiplyMaxima(v, maxima);
    return pointset;
}
template <typename C, typename M>
std::vector<C> filterDimSweep2D(std::vector<C> &v, M const &maxima) {
    /*Calculates the non dominated points on a 2d set using a matrix of vectors
    as its core data structure Complexity : NlogN
    */
    std::vector<C> pointset;
    multiplyMaxima(v, maxima);

    std::sort(v.begin(), v.end(),
              CompleteColSortReverse());  // Sort in x in reverse order
    double maxy = v[0][1];                // Y coord of the last point in x
    size_t idx = 1;
    for (size_t point = 1; point < v.size(); point++) {
        if (v[point][1] > maxy) {
            pointset.push_back(v[point]);
            maxy = v[point][1];
        }
    }
    multiplyMaxima(pointset, maxima);
    multiplyMaxima(v, maxima);
    return pointset;
}
}  // namespace notinplace
}  // namespace nondlib
#endif