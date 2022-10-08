// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <limits>
#include "Point.hpp"



template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;


  KDTreeNode<N, ElemType>* root;
  int depth; //profundidad
  Point<N>* guess = nullptr; //punto potencialmente mas cercano
  float bestDist = std::numeric_limits<float>::max();

 private:
  size_t dimension_;//dimension de los puntos
  size_t size_; //# de elementos almacenados
  struct KDTreeNode {
      Point<N> punto;
      KDTreeNode* childs[2]{};
      KDTreeNode() {
          this->childs[0] = nullptr;
          this->childs[1] = nullptr;
      }
      KDTreeNode(Point<N> pt){
          this->punto = pt;
          this->childs[0] = nullptr;
          this->childs[1] = nullptr;
      }
      KDTreeNode(Point<N> pt, KDTreeNode* left, KDTreeNode* right) {
          this->punto = pt;
          this->childs[0] = left;
          this->childs[1] = right;
      }

      KDTreeNode(const value_type& value);

      const Point<N>& getPoint() const {
          return punto;
      }
      void SetPoint(const Point<N>& pt) const {
          punto = pt;
      }
      KDTreeNode* const* getChild() const {
          return childs;
      }

  }

};

template <size_t N, typename ElemType>
//typename  KDTree<N, ElemType>::KDTreeNode(const value_type& value);

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  // TODO(me): Fill this in.
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  // TODO(me): Fill this in.
  return 0;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  // TODO(me): Fill this in.
  return 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  // TODO(me): Fill this in.
  return true;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  // TODO(me): Fill this in.
  return true;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
