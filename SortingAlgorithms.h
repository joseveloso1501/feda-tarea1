#include <vector>
#include <string>

class SortingAlgorithms
{
public:
    void bubbleSort(std::vector<int> &array);
    void mergeSort(std::vector<int> &array, int left, int right);
    void quickSort(std::vector<int> &array, int low, int high);
    void insertionSort(std::vector<int> &array);

private:
    void merge(std::vector<int> &array, int left, int middle, int right);
    int partition(std::vector<int> &array, int low, int high);
};