#ifndef MAX_HEAP
#define MAX_HEAP

#include <structDist.hpp>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include <stdexcept>


class MaxHeap {
public:
    void insert(const structDist& element) {
        heap.push(element);
        id_map[element.id] = element;
    }
    structDist getMax() {
        while (!heap.empty() && id_map.find(heap.top().id) == id_map.end()) {
            heap.pop();
        }
        if (heap.empty()) {
            throw std::runtime_error("Heap is empty");
        }
        return heap.top();
    }

    std::vector<structDist> getTwoMax() {
        std::vector<structDist> result;
        std::priority_queue<structDist> tempHeap = heap; // Copiar el heap

        // Obtener el primer máximo
        if (!tempHeap.empty()) {
            result.push_back(tempHeap.top());
            tempHeap.pop();
        }

        while (!tempHeap.empty() && id_map.find(tempHeap.top().id) == id_map.end()) {
            tempHeap.pop();
        }
        if (!tempHeap.empty()) {
            result.push_back(tempHeap.top());
        }

        return result;
    }

    void extractMax() {
        if (heap.empty()) {
            throw std::runtime_error("Heap is empty");
        }
        structDist maxElement = heap.top();
        heap.pop();
        id_map.erase(maxElement.id);
    }
    void removeById(int id) {
        auto it = id_map.find(id);
        if (it != id_map.end()) {
            id_map.erase(it);
        }
    }

private:
    std::priority_queue<structDist> heap;
    std::unordered_map<int, structDist> id_map;
};


#endif