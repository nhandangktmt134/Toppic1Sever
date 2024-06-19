#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    int label;
    kDTreeNode *left;
    kDTreeNode *right;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->label = 0;
        this->left = left;
        this->right = right;
    }
    kDTreeNode(vector<int> data, int label, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->label = label;
        this->left = nullptr;
        this->right = nullptr;
    }

    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (int i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    int count;

public:
    kDTree(int k = 2) : k(k), root(nullptr), count(0) {}
    kDTreeNode *copyTree(kDTreeNode *root)
    {
        if (root == nullptr)
            return nullptr;
        kDTreeNode *newNode = new kDTreeNode(root->data, root->label);
        newNode->left = copyTree(root->left);
        newNode->right = copyTree(root->right);

        return newNode;
    }
    void deleteTree(kDTreeNode *node)
    {
        if (node == nullptr)
            return;
        deleteTree(node->left);
        deleteTree(node->right);
        delete node;
    }
    ~kDTree()
    {
        deleteTree(root);
    }

    const kDTree &operator=(const kDTree &other)
    {
        this->k = other.k;
        deleteTree(root);
        this->root = copyTree(other.root);
        this->count = other.count;
        return *this;
    }
    kDTree(const kDTree &other)
    {
        this->k = other.k;
        deleteTree(root);
        this->root = copyTree(other.root);
        this->count = other.count;
    }
    int nodeCount(kDTreeNode *root) const
    {
        if (root == nullptr)
            return 0;

        return nodeCount(root->left) + nodeCount(root->right) + 1;
    }
    int nodeCount() const
    {
        return nodeCount(root);
    }

    int recHeight(kDTreeNode *root) const
    {
        if (root == nullptr)
            return 0;

        int leftHeight = recHeight(root->left);
        int rightHeight = recHeight(root->right);
        return max(leftHeight, rightHeight) + 1;
    }
    int height() const
    {
        return recHeight(root);
    }
    int recLeafCount(kDTreeNode *root) const
    {
        if (root == nullptr)
            return 0;
        if (root->left == nullptr && root->right == nullptr)
            return 1;

        int leftLeaves = recLeafCount(root->left);
        int rightLeaves = recLeafCount(root->right);
        return leftLeaves + rightLeaves;
    }
    int leafCount() const
    {
        return recLeafCount(root);
    }
    void recInorderTraversal(kDTreeNode *root) const
    {
        // Implement inorderTraversal logic here
        if (root == nullptr)
        {
            return;
        }
        recInorderTraversal(root->left);
        cout << *root << " ";
        recInorderTraversal(root->right);
    }
    void inorderTraversal() const
    {
        recInorderTraversal(root);
    }
    void recPreorderTraversal(kDTreeNode *root) const
    {
        // Implement preorderTraversal logic here
        if (root == nullptr)
        {
            return;
        }
        cout << *root << " ";
        recPreorderTraversal(root->left);
        recPreorderTraversal(root->right);
    }
    void preorderTraversal() const
    {
        recPreorderTraversal(root);
    }
    void recPostorderTraversal(kDTreeNode *root) const
    {
        // Implement postorderTraversal logic here
        if (root == nullptr)
        {
            return;
        }
        recPostorderTraversal(root->left);
        recPostorderTraversal(root->right);
        cout << *root << " ";
    }
    void postorderTraversal() const
    {
        recPostorderTraversal(root);
    }
    kDTreeNode *recInsert(kDTreeNode *temp, const vector<int> &point, int level) const
    {
        // Bước 1: Điều kiện dừng
        if (temp == nullptr)
        {
            return new kDTreeNode(point);
        }

        // Bước 2: Tính alpha
        int alpha = level % k;

        // Bước 3 và Bước 4: So sánh và gọi đệ quy
        if (point[alpha] < temp->data[alpha])
        {
            temp->left = recInsert(temp->left, point, level + 1);
        }
        else
        {
            temp->right = recInsert(temp->right, point, level + 1);
        }

        return temp;
    }

    void insert(const vector<int> &point)
    {
        this->root = this->recInsert(root, point, 0);
        this->count++;
    }
    kDTreeNode *minNode(kDTreeNode *a, kDTreeNode *b, kDTreeNode *c, int d) const
    {
        kDTreeNode *temp = a;
        if (b != nullptr && b->data[d] < temp->data[d])
        {
            temp = b;
        }
        if (c != nullptr && c->data[d] < temp->data[d])
        {
            temp = c;
        }
        return temp;
    }
    kDTreeNode *findMin(kDTreeNode *root, int d, int level) const
    {
        if (root == nullptr)
            return nullptr;

        int alpha = level % k;

        if (alpha == d)
        {
            if (root->left == nullptr)
            {
                return root; // Không có cây con bên trái, trả về node hiện tại
            }
            else
            {
                return findMin(root->left, d, level + 1);
            }
        }
        else
        {
            return minNode(root, findMin(root->left, d, level + 1), findMin(root->right, d, level + 1), d);
        }
    }
    kDTreeNode *recRemove(kDTreeNode *node, const vector<int> &point, int level) const
    {
        if (node == nullptr)
            return nullptr;
        int alpha = level % k;
        if (node->data == point)
        {
            if (node->right != nullptr)
            {
                // Tìm node nho nhat trong cây con bên phải cho chieu hiện tại
                kDTreeNode *min = findMin(node->right, alpha, level + 1);
                node->data = min->data; // Replace node data with min
                node->right = recRemove(node->right, min->data, level + 1);
            }
            else if (node->left != nullptr)
            {
                kDTreeNode *min = findMin(node->left, alpha, level + 1);
                node->data = min->data; // Replace node data with min
                node->right = node->left;
                node->left = nullptr;
                node->right = recRemove(node->right, min->data, level + 1);
            }
            else
            {
                // Leaf node
                delete node;
                return nullptr;
            }
            return node;
        }
        if (point[alpha] < node->data[alpha])
        {
            node->left = recRemove(node->left, point, level + 1);
        }
        else
        {
            node->right = recRemove(node->right, point, level + 1);
        }

        return node;
    }
    void remove(const vector<int> &point)
    {
        this->root = this->recRemove(root, point, 0);
        this->count--;
    }
    bool searchRec(kDTreeNode *node, const vector<int> &point, int level) const
    {
        if (node == nullptr)
            return false;

        int alpha = level % k;
        if (point[alpha] < node->data[alpha])
        {
            return searchRec(node->left, point, level + 1);
        }
        else if (point[alpha] > node->data[alpha])
        {
            return searchRec(node->right, point, level + 1);
        }
        else
        {
            if (node->data == point)
                return true; // Point found
            else
            {
                bool leftNode = searchRec(node->left, point, level + 1);
                bool rightNode = searchRec(node->right, point, level + 1);
                return leftNode || rightNode;
            }
        }
    }
    bool search(const vector<int> &point)
    {
        return searchRec(root, point, 0);
    }
    void merge(vector<vector<int>> &points, int left, int mid, int right, int dim) const
    {
        int n1 = mid - left + 1;
        int n2 = right - mid;
        vector<vector<int>> leftPoints(n1), rightPoints(n2);

        for (int i = 0; i < n1; i++)
        {
            leftPoints[i] = points[left + i];
        }
        for (int j = 0; j < n2; j++)
        {
            rightPoints[j] = points[mid + 1 + j];
        }

        int i = 0, j = 0, k = left;
        while (i < n1 && j < n2)
        {
            if (leftPoints[i][dim] <= rightPoints[j][dim])
            {
                points[k] = leftPoints[i];
                i++;
            }
            else
            {
                points[k] = rightPoints[j];
                j++;
            }
            k++;
        }

        while (i < n1)
        {
            points[k] = leftPoints[i];
            i++;
            k++;
        }
        while (j < n2)
        {
            points[k] = rightPoints[j];
            j++;
            k++;
        }
    }
    void mergeSort(vector<vector<int>> &points, int l, int r, int dim) const
    {
        if (l >= r)
            return;

        // Tính toán chỉ số trung gian để chia mảng
        int m = l + (r - l) / 2;

        // Gọi đệ quy cho nửa trái và nửa phải của mảng
        mergeSort(points, l, m, dim);
        mergeSort(points, m + 1, r, dim);

        // Gộp hai nửa đã được sắp xếp
        merge(points, l, m, r, dim);
    }

    kDTreeNode *recBuildTree(const vector<vector<int>> &pointList, int level)
    {
        // Bước 1: kiem tra dieu kien dung
        if (pointList.empty())
        {
            return nullptr;
        }

        // Buoc 2: sap xep mang theo chieu alpha
        int alpha = level % k;
        vector<vector<int>> sortedPoints = pointList;
        mergeSort(sortedPoints, 0, sortedPoints.size() - 1, alpha);

        // B3: tim trung vi
        int medianIndex = (sortedPoints.size() - 1) / 2;
        while (medianIndex > 0 && sortedPoints[medianIndex][alpha] == sortedPoints[medianIndex - 1][alpha])
        {
            medianIndex--;
        }
        // Bước 4: Chia cay
        vector<vector<int>> leftPoints(sortedPoints.begin(), sortedPoints.begin() + medianIndex);
        vector<vector<int>> rightPoints(sortedPoints.begin() + medianIndex + 1, sortedPoints.end());

        // Bước 5: Khoi tao cay
        kDTreeNode *node = new kDTreeNode(sortedPoints[medianIndex]);
        this->count++;
        node->left = recBuildTree(leftPoints, level + 1);
        node->right = recBuildTree(rightPoints, level + 1);

        return node;
    }
    void buildTree(const vector<vector<int>> &pointList)
    {
        this->root = this->recBuildTree(pointList, 0);
    }
    double distance(const vector<int> &a, const vector<int> &b)
    {
        double dist = 0.0;
        for (size_t i = 0; i < a.size(); i++)
        {
            dist += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return sqrt(dist);
    }
    int distanceSquare(const vector<int> &a, const vector<int> &b)
    {
        int dist = 0;
        for (size_t i = 0; i < a.size(); i++)
        {
            dist += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return dist;
    }
    /////////
    void recNearestNeighbour(kDTreeNode *node, const vector<int> &target, kDTreeNode *&best, int level, int k)
    {
        // Bước 1: Dừng lại nếu node hiện tại là null
        if (node == nullptr)
            return;

        // Bước 2: Tính alpha
        int alpha = level % k;
        kDTreeNode *nextBranch = (target[alpha] < node->data[alpha]) ? node->left : node->right;
        kDTreeNode *oppositeBranch = (target[alpha] < node->data[alpha]) ? node->right : node->left;

        // Di về phía có giá trị theo chiều alpha gần với target hơn
        recNearestNeighbour(nextBranch, target, best, level + 1, k);

        // Bước 3: Tính khoảng cách giữa target và node
        double tempDistance = distance(node->data, target);
        double bestDistance = (best != nullptr) ? distance(best->data, target) : 1.79769e+308;

        // Bước 4: Nếu khoảng cách giữa target và node nhỏ hơn khoảng cách gần nhất
        if (tempDistance < bestDistance)
        {
            best = node;
            bestDistance = tempDistance;
        }

        // Kiểm tra xem có cần khám phá nhánh đối diện không
        double d = abs(target[alpha] - node->data[alpha]);
        if (d <= bestDistance)
        {
            recNearestNeighbour(oppositeBranch, target, best, level + 1, k);
        }
    }

    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best)
    {
        best = nullptr;
        this->recNearestNeighbour(root, target, best, 0, k);
    }
    void sortBestList(vector<kDTreeNode *> &bestList, const vector<int> &target)
    {
        int n = bestList.size();
        bool swapped;
        for (int i = 0; i < n - 1; i++)
        {
            swapped = false;
            for (int j = 0; j < n - i - 1; j++)
            {
                int dist1 = distanceSquare(bestList[j]->data, target);
                int dist2 = distanceSquare(bestList[j + 1]->data, target);
                if (dist1 > dist2)
                {
                    swap(bestList[j], bestList[j + 1]);
                    swapped = true;
                }
            }
            if (!swapped)
            {
                break;
            }
        }
    }
    void kNearestNeighbourREC(kDTreeNode *temp, const vector<int> &target, int k, vector<kDTreeNode *> &bestList, int level)
    {
        if (temp == nullptr)
            return;
        int alpha = level % this->k;
        kDTreeNode *nextBranch = (target[alpha] < temp->data[alpha]) ? temp->left : temp->right;
        kDTreeNode *oppositeBranch = (target[alpha] < temp->data[alpha]) ? temp->right : temp->left;
        kNearestNeighbourREC(nextBranch, target, k, bestList, level + 1);
        int d = (target[alpha] - temp->data[alpha]) * (target[alpha] - temp->data[alpha]);
        // cout<<"Target "<<target[alpha] <<"+ "<<"Temp "<<temp->data[alpha]<<endl;
        int dist = distanceSquare(temp->data, target);
        if (bestList.size() < k)
        {
            bestList.push_back(temp);
            sortBestList(bestList, target);
        }
        else if (dist < distanceSquare(bestList.back()->data, target))
        {
            bestList.pop_back();
            bestList.push_back(temp);
            sortBestList(bestList, target);
        }
        if (d <= dist || bestList.size() < k)
        {
            kNearestNeighbourREC(oppositeBranch, target, k, bestList, level + 1);
        }
    }
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList)
    {
        kNearestNeighbourREC(root, target, k, bestList, 0);
    }
    void mergePointLabel(vector<vector<int>> &points, vector<int> &labels, int left, int mid, int right, int dim) const
    {
        int n1 = mid - left + 1;
        int n2 = right - mid;
        vector<vector<int>> leftPoints(n1), rightPoints(n2);
        vector<int> leftLabels(n1), rightLabels(n2);

        for (int i = 0; i < n1; i++)
        {
            leftPoints[i] = points[left + i];
            leftLabels[i] = labels[left + i]; // Sao chép label tương ứng với leftPoints
        }
        for (int j = 0; j < n2; j++)
        {
            rightPoints[j] = points[mid + 1 + j];
            rightLabels[j] = labels[mid + 1 + j]; // Sao chép label tương ứng với rightPoints
        }

        int i = 0, j = 0, k = left;
        while (i < n1 && j < n2)
        {
            if (leftPoints[i][dim] <= rightPoints[j][dim])
            {
                points[k] = leftPoints[i];
                labels[k] = leftLabels[i]; // Cập nhật label tương ứng với leftPoints
                i++;
            }
            else
            {
                points[k] = rightPoints[j];
                labels[k] = rightLabels[j]; // Cập nhật label tương ứng với rightPoints
                j++;
            }
            k++;
        }

        while (i < n1)
        {
            points[k] = leftPoints[i];
            labels[k] = leftLabels[i]; // Cập nhật label tương ứng với leftPoints
            i++;
            k++;
        }
        while (j < n2)
        {
            points[k] = rightPoints[j];
            labels[k] = rightLabels[j]; // Cập nhật label tương ứng với rightPoints
            j++;
            k++;
        }
    }

    void mergeSortPointLabel(vector<vector<int>> &points, vector<int> &labels, int l, int r, int dim) const
    {
        if (l >= r)
            return;

        // Tính toán chỉ số trung gian để chia mảng
        int m = l + (r - l) / 2;

        // Gọi đệ quy cho nửa trái và nửa phải của mảng
        mergeSortPointLabel(points, labels, l, m, dim);
        mergeSortPointLabel(points, labels, m + 1, r, dim);

        // Gộp hai nửa đã được sắp xếp
        mergePointLabel(points, labels, l, m, r, dim);
    }
    kDTreeNode *recBuildTreeLabel(const vector<vector<int>> &pointList, const vector<int> &label, int level)
    {
        // Bước 1: kiem tra dieu kien dung
        if (pointList.empty())
        {
            return nullptr;
        }

        // Buoc 2: sap xep mang theo chieu alpha
        int alpha = level % k;
        vector<vector<int>> sortedPoints = pointList;
        vector<int> sortedLabels = label;
        mergeSortPointLabel(sortedPoints, sortedLabels, 0, sortedPoints.size() - 1, alpha);

        // B3: tim trung vi
        int medianIndex = (sortedPoints.size() - 1) / 2;
        while (medianIndex > 0 && sortedPoints[medianIndex][alpha] == sortedPoints[medianIndex - 1][alpha])
        {
            medianIndex--;
        }
        // Bước 4: Chia cay
        vector<vector<int>> leftPoints(sortedPoints.begin(), sortedPoints.begin() + medianIndex);
        vector<vector<int>> rightPoints(sortedPoints.begin() + medianIndex + 1, sortedPoints.end());
        vector<int> leftLabel(sortedLabels.begin(), sortedLabels.begin() + medianIndex);
        vector<int> rightLabel(sortedLabels.begin() + medianIndex + 1, sortedLabels.end());
        // Bước 5: Khoi tao cay
        kDTreeNode *node = new kDTreeNode(sortedPoints[medianIndex], sortedLabels[medianIndex]);
        this->count++;
        node->left = recBuildTreeLabel(leftPoints, leftLabel, level + 1);
        node->right = recBuildTreeLabel(rightPoints, rightLabel, level + 1);

        return node;
    }
    void buildTreeLabel(const vector<vector<int>> &pointList, const vector<int> &label)
    {
        this->count = pointList.size();
        this->root = this->recBuildTreeLabel(pointList, label, 0);
    };
    friend class kNN;
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    kDTree kdtree;

public:
    kNN(int k = 5) : k(k){};
    vector<vector<int>> listListToVectorVector(const list<list<int>> &data)
    {
        vector<vector<int>> result;
        // Iterate over each inner list in the outer list
        for (const auto &innerList : data)
        {
            vector<int> innerVec;
            // Iterate over each integer in the inner list
            for (const auto &num : innerList)
            {
                // Add the integer to the inner vector
                innerVec.push_back(num);
            }
            // Add the inner vector to the result vector of vectors
            result.push_back(innerVec);
        }
        return result;
    }
    vector<int> listListToVector(const list<list<int>> &data)
    {
        vector<int> result;
        // Iterate over each inner list in the outer list
        for (const auto &innerList : data)
        {
            // Iterate over each integer in the inner list
            for (const auto &num : innerList)
            {
                // Add the integer to the result vector
                result.push_back(num);
            }
        }
        return result;
    }
    void fit(Dataset &X_train, Dataset &y_train)
    {
        this->X_train = &X_train;
        this->y_train = &y_train;
        if (X_train.data.size())
        {
            int dismension = X_train.data.front().size();
            kdtree.k = dismension;
            vector<vector<int>> pointList;
            vector<int> label;
            //
            pointList = listListToVectorVector(X_train.data);
            label = listListToVector(y_train.data);
            kdtree.buildTreeLabel(pointList, label);
        }
    };
    Dataset predict(Dataset &X_test)
    {
        Dataset result;
        result.columnName.push_back("label");
        vector<vector<int>> target;
        target = listListToVectorVector(X_test.data);
        // get
        for (const auto &targetVec : target)
        {
            vector<kDTreeNode *> best;
            kdtree.kNearestNeighbour(targetVec, this->k, best);
            vector<int> labelCounts(10, 0); // Initialize with zeros for labels 0 to 9
            for (auto node : best)
            {
                labelCounts[node->label]++;
            }
            // Find the most frequent label
            int mostFrequentLabel = 0; // Default to label 0
            int maxCount = 0;
            for (int i = 0; i < 10; ++i)
            { // Iterate through labels 0 to 9
                if (labelCounts[i] > maxCount || (labelCounts[i] == maxCount && i < mostFrequentLabel))
                {
                    mostFrequentLabel = i;
                    maxCount = labelCounts[i];
                }
            }
            result.data.push_back({mostFrequentLabel});
        }

        return result;
    };
    double score(const Dataset &y_test, const Dataset &y_pred)
    {
        if (y_test.data.empty() || y_pred.data.empty())
        {
            return -1; // Return 0 if any dataset is empty
        }

        auto it_test = y_test.data.begin();
        auto it_pred = y_pred.data.begin();
        int count = 0;
        size_t numSamples = y_test.data.size();

        while (it_test != y_test.data.end() && it_pred != y_pred.data.end())
        {
            int labelTest = it_test->front(); // Get label from y_test
            int labelPred = it_pred->front(); // Get label from y_pred
            if (labelTest == labelPred)
            {
                count++; // Increment count if labels match
            }
            ++it_test;
            ++it_pred;
        }
        return ((double)count / numSamples);
    };

    //* test case thôi
    // void print_Y(const Dataset &y)
    // {
    //     cout << y.columnName[0] << ": ";
    //     for (auto it : y.data)
    //     {
    //         cout << it.front() << " ";
    //     }
    //     cout << endl;
    // }
};

// Please add more or modify as needed
