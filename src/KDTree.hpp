// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <algorithm>
#include <map>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <limits>
#include <vector>
#include <queue>
#include <stdlib.h>
#include "Point.hpp"


/*PRIORITY QUEUE*/
template <typename T>
class PQueueM {
public:
   
    explicit PQueueM(size_t maxSize);
    void enqueue(const T& value, double priority);
    T dequeueMin();
    size_t size() const;
    bool empty() const;
    size_t maxSize() const;
    double best()  const;
    double worst() const;

private:
    std::multimap<double, T> elems;
    size_t SizeMax;
};

template <typename T>
PQueueM<T>::PQueueM(size_t maxSize) {
    SizeMax = maxSize;
}

template <typename T>
void PQueueM<T>::enqueue(const T& value, double priority) {
    // Add the element
    elems.insert(std::make_pair(priority, value));

    // If there are too many elements in the queue, drop off the last one.
    if (size() > maxSize()) {
        typename std::multimap<double, T>::iterator last = elems.end();
        --last;
        elems.erase(last);
    }
}

template <typename T>
T PQueueM<T>::dequeueMin() {
    // Copy the best value.
    T result = elems.begin()->second;
    elems.erase(elems.begin());

    return result;
}

template <typename T>
size_t PQueueM<T>::size() const {
    return elems.size();
}

template <typename T>
bool PQueueM<T>::empty() const {
    return elems.empty();
}

template <typename T>
size_t PQueueM<T>::maxSize() const {
    return SizeMax;
}

template <typename T>
double PQueueM<T>::best() const {
    return empty() ? numeric_limits<double>::infinity() : elems.begin()->first;
}

template <typename T>
double PQueueM<T>::worst() const {
    return empty() ? std::numeric_limits<double>::infinity() : elems.rbegin()->first;
}

/*End priority queue details*/


template <size_t N, typename ElemType>
class KDTree
{

public:
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree& rhs);

    KDTree& operator=(const KDTree& rhs);

    size_t dimension() const;

    size_t size() const;
    bool empty() const;

    bool contains(const Point<N>& pt) const;

    void insert(const Point<N>& pt, const ElemType& value);

    ElemType& operator[](const Point<N>& pt);

    ElemType& at(const Point<N>& pt);

    const ElemType& at(const Point<N>& pt) const;

    ElemType knn_value(const Point<N>& key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

private:
    ElemType def;
    size_t dimension_;//point dimension 
    size_t size_;     //number of stored points
    Point<N>* key = nullptr; //black point
    Point<N>* guess = nullptr; //potentially nearest point
    //double bestDist = std::numeric_limits<double>::infinity();
    int k = 1; //k nearest points

    struct KDTreeNode
    {
        Point<N> pt_;
        ElemType value_;
        KDTreeNode* children[2];

        KDTreeNode(const Point<N>& point_, const ElemType& value)
        {

            children[0] = nullptr;
            children[1] = nullptr;
            pt_ = point_;
            value_ = value;
        }

        KDTreeNode(Point<N> pt, ElemType value, KDTreeNode* c_left, KDTreeNode* c_right)
        {
            pt_ = pt;
            value_ = value;
            children[0] = c_left;
            children[1] = c_right;
        }
        ~KDTreeNode(){
            //delete children[0];
            //delete children[1];
        }
    };
    
    bool search(const Point<N>& point_, KDTreeNode**& result)
    {
        int dimension{ 0 };
        result = &root;

        while (*result && (*result)->pt_ != point_)
        {
            dimension = dimension > dimension_ - 1 ? 0 : dimension;
            result = &((*result)->children[point_[dimension] > (*result)->pt_[dimension]]);
            dimension++;
        }

        return (*result) != nullptr;
    }
    bool search(const Point<N>& point_, KDTreeNode**& result) const
    {
        int dimension{ 0 };

        result = &root;
        while (*result && (*result)->pt_ != point_)
        {
            dimension = dimension > dimension_ - 1 ? 0 : dimension;
            result = &((*result)->children[point_[dimension] > (*result)->pt_[dimension]]);
            dimension++;
        }
        return (*result) != nullptr;
    }
    mutable KDTreeNode* root; //variable
    void destroy(KDTreeNode* node)
    {
        if (node != nullptr)
        {
            destroy(node->children[0]);
            destroy(node->children[1]);
            delete node;
        }
    }
    static KDTreeNode* copyNodes(const KDTreeNode* node)
    {
        if (node != nullptr)
        {
            KDTreeNode* nodeCopy = new KDTreeNode(node->pt_, node->value_, copyNodes(node->children[0]), copyNodes(node->children[1]));
            return nodeCopy;
        }
    }
    
    void nearest_neighbor(Point<N> key, KDTreeNode* currentNode, ElemType*& guest, double& bestDist, int depth) const;
    void nearest_neighborPQ(Point<N> key, KDTreeNode* currentNode, int depth, PQueueM<ElemType> &pqP ) const;
    ElemType nneighbor(const Point<N> &key) const {
        double bestDist = std::numeric_limits<double>::infinity();
        ElemType* val = nullptr;
        nearest_neighbor(key, root, val, bestDist, 0);
        //return val;
        return val == nullptr ? ElemType() : *val;
    }
};

template <size_t N, typename ElemType>//  KDTreeNode*& guest
void KDTree<N, ElemType>::nearest_neighbor(Point<N> key, KDTreeNode* currentNode, ElemType*& guest, double& bestDist, int depth)
const
{
    if (!currentNode)
        return;
    double currentDistance = distance(currentNode->pt_, key);

    if (currentDistance < bestDist)
    {
        bestDist = currentDistance;
        guest = &currentNode->value_;
    }
    int axis = depth % dimension_;
    bool child = key[axis] < currentNode->pt_[axis];

    nearest_neighbor(key, currentNode->children[!child], guest, bestDist, ++depth);
    //if (fabs(currentNode->pt_[axis] - key[axis]) < guest)
    if (fabs(currentNode->pt_[axis] - key[axis]) < bestDist)
        nearest_neighbor(key, currentNode->children[child], guest, bestDist, ++depth);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::nearest_neighborPQ(Point<N> key, KDTreeNode* currentNode, int depth, PQueueM<ElemType>& pqP) const
{
    if (!currentNode)
        return;

    pqP.enqueue(currentNode->value_, distance(currentNode->pt_, key));

    int axis = depth % dimension_;
    bool child = key[axis] < currentNode->pt_[axis];

    nearest_neighborPQ(key, currentNode->children[!child], ++depth,pqP);

    if (pqP.size() < pqP.maxSize() || fabs(currentNode->pt_[axis] - key[axis]) < pqP.worst())
        nearest_neighborPQ(key, currentNode->children[child], ++depth, pqP);
    
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree()
{
    size_ = 0;
    root = nullptr;
    dimension_ = N;
    def = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree()
{
    destroy(root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs)
{
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
    def = rhs.def;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs)
{
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
    def = rhs.def;
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const
{
    return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const
{
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const
{
    return size_ == 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const
{
    KDTreeNode** tmp;
    bool result = search(pt, tmp);
    return result;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value)
{
    KDTreeNode** tmp;
    if (empty())
    {
        ++size_;
        this->root = new KDTreeNode(pt, value);
    }
    else if (search(pt, tmp))
    {
        at(pt) = value;
    }
    else
    {
        ++size_;
        (*tmp) = new KDTreeNode(pt, value);
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt)
{
    KDTreeNode** tmp;
    if (search(pt, tmp))
        return (*tmp)->value_;
    insert(pt, def);
    search(pt, tmp);
    return (*tmp)->value_;
    ;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt)
{
    KDTreeNode** tmp;
    if (!search(pt, tmp))
    {
        throw std::out_of_range("");
    }
    return ((*tmp)->value_);
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const
{
    KDTreeNode** tmp;
    if (!search(pt, tmp))
    {
        throw std::out_of_range("The point is out of range");
    }
    return ((*tmp)->value_);
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const
{
    if (k == 1)
        //return nearest_neighbor(key, root, guess, bestDist, 0);
        return nneighbor(key);

    PQueueM<ElemType> pqPoints(k);
    nearest_neighborPQ(key, root, 0, pqPoints);

    std::multimap<int, ElemType, std::greater<int>> MapF;

    while (!pqPoints.empty()) {
        ElemType element = pqPoints.dequeueMin();
        for (typename std::multimap<int, ElemType>::iterator it = MapF.begin(); it != MapF.end(); ++it)
        {
            if (it->second == element) {
                MapF.insert(std::make_pair(it->first + 1, it->second));
                MapF.erase(it);
                goto salir;//break;
            }
        }
        MapF.insert(std::make_pair(1,element));
    salir: ;
    }
    return MapF.begin()->second;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key, size_t k) const
{
    std::vector<ElemType> values;
    if (k == 1) {
        values.push_back(nneighbor(key));
    }
    else {
        ElemType* resp = knn_value(key, k);
        for (auto ite = resp.begin(); ite != resp.end(); ite++) {
            values.push_back(ite);
        }
    }
    return values;
}

#endif // SRC_KDTREE_HPP_