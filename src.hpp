#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************
    
    friend class resistive_network;

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************
    
    int get_rows() const { return m; }
    int get_cols() const { return n; }

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // TODO: 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) {
        m = m_;
        n = n_;
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    // TODO: 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) {
        m = obj.m;
        n = obj.n;
        if (obj.data == nullptr) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; i++) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; j++) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
    }

    // TODO: 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // TODO: 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // TODO: 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        
        // 先释放原有内存
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
        
        m = obj.m;
        n = obj.n;
        if (obj.data == nullptr) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; i++) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; j++) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
        return *this;
    }

    // TODO: 重载括号，返回矩阵的第i行(1-based)、第j列(0-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n) {
            throw matrix_error();
        }
        return data[i-1][j];
    }

    // TODO: 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                result.data[i][j] = fraction(0);
                for (int k = 0; k < lhs.n; k++) {
                    result.data[i][j] = result.data[i][j] + lhs.data[i][k] * rhs.data[k][j];
                }
            }
        }
        return result;
    }

    // TODO: 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (data == nullptr || m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // TODO: 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (data == nullptr || m == 0 || n == 0 || m != n) {
            throw matrix_error();
        }
        
        // 创建临时矩阵用于高斯消元
        matrix temp(*this);
        fraction det(1);
        
        // 高斯消元
        for (int i = 0; i < m; i++) {
            // 找到第i列的主元
            int pivot = i;
            for (int j = i + 1; j < m; j++) {
                // 这里需要比较绝对值，但fraction类没有提供比较运算符
                // 我们只需要找到非零元素
                if (!(temp.data[pivot][i] == fraction(0)) || !(temp.data[j][i] == fraction(0))) {
                    if (temp.data[pivot][i] == fraction(0)) {
                        pivot = j;
                    }
                }
            }
            
            // 如果主元为0，行列式为0
            if (temp.data[pivot][i] == fraction(0)) {
                return fraction(0);
            }
            
            // 交换行
            if (pivot != i) {
                for (int j = 0; j < n; j++) {
                    fraction t = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot][j];
                    temp.data[pivot][j] = t;
                }
                det = det * fraction(-1);
            }
            
            // 消元
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0))) {
                    fraction factor = temp.data[j][i] / temp.data[i][i];
                    for (int k = i; k < n; k++) {
                        temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                    }
                }
            }
        }
        
        // 计算对角线元素的乘积
        for (int i = 0; i < m; i++) {
            det = det * temp.data[i][i];
        }
        
        return det;
    }
};

class resistive_network {
private:

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 矩阵A 和 矩阵C
    matrix adjacency, conduction;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************
    
    // 解线性方程组 Ax = b，返回x
    matrix solve_linear_system(const matrix &A, const matrix &b) {
        int n = A.get_rows();
        // 创建增广矩阵 [A|b]
        matrix aug(n, n + 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 0; j < n; j++) {
                aug(i, j) = A.data[i-1][j];
            }
            aug(i, n) = b.data[i-1][0];
        }
        
        // 高斯-约旦消元
        for (int i = 0; i < n; i++) {
            // 找主元
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (!(aug.data[j][i] == fraction(0))) {
                    if (aug.data[pivot][i] == fraction(0)) {
                        pivot = j;
                    }
                }
            }
            
            // 交换行
            if (pivot != i) {
                for (int j = 0; j <= n; j++) {
                    fraction t = aug.data[i][j];
                    aug.data[i][j] = aug.data[pivot][j];
                    aug.data[pivot][j] = t;
                }
            }
            
            // 归一化当前行
            fraction div = aug.data[i][i];
            for (int j = 0; j <= n; j++) {
                aug.data[i][j] = aug.data[i][j] / div;
            }
            
            // 消元
            for (int j = 0; j < n; j++) {
                if (j != i && !(aug.data[j][i] == fraction(0))) {
                    fraction factor = aug.data[j][i];
                    for (int k = 0; k <= n; k++) {
                        aug.data[j][k] = aug.data[j][k] - factor * aug.data[i][k];
                    }
                }
            }
        }
        
        // 提取解
        matrix x(n, 1);
        for (int i = 0; i < n; i++) {
            x.data[i][0] = aug.data[i][n];
        }
        return x;
    }

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

    // TODO: 设置电阻网络，构建矩阵A和C。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        
        // 构建关联矩阵 A (n x m)
        adjacency = matrix(interface_size, connection_size);
        for (int i = 0; i < connection_size; i++) {
            adjacency(from[i], i) = fraction(1);
            adjacency(to[i], i) = fraction(-1);
        }
        
        // 构建电导矩阵 C (m x m)，对角矩阵
        conduction = matrix(connection_size, connection_size);
        for (int i = 0; i < connection_size; i++) {
            conduction(i+1, i) = fraction(1) / resistance[i];
        }
    }

    ~resistive_network() = default;

    // TODO: 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) {
            return fraction(0);
        }
        
        // 设置电流：在interface_id1注入1A，在interface_id2流出1A
        matrix I(interface_size, 1);
        I(interface_id1, 0) = fraction(1);
        I(interface_id2, 0) = fraction(-1);
        
        // 构建导纳矩阵 G = A * C * A^T
        matrix AT = adjacency.transposition();
        matrix G = adjacency * conduction * AT;
        
        // 去掉最后一行和最后一列（因为u_n = 0）
        int n = interface_size - 1;
        matrix G_reduced(n, n);
        matrix I_reduced(n, 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                G_reduced(i, j-1) = G(i, j-1);
            }
            I_reduced(i, 0) = I(i, 0);
        }
        
        // 解方程 G_reduced * U = I_reduced
        matrix U = solve_linear_system(G_reduced, I_reduced);
        
        // 计算电压差
        fraction u1 = (interface_id1 <= n) ? U(interface_id1, 0) : fraction(0);
        fraction u2 = (interface_id2 <= n) ? U(interface_id2, 0) : fraction(0);
        
        // 等效电阻 = 电压差 / 电流
        return u1 - u2;
    }

    // TODO: 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]) {
        // 构建导纳矩阵 G = A * C * A^T
        matrix AT = adjacency.transposition();
        matrix G = adjacency * conduction * AT;
        
        // 构建电流向量（去掉最后一个节点）
        int n = interface_size - 1;
        matrix G_reduced(n, n);
        matrix I_reduced(n, 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                G_reduced(i, j-1) = G(i, j-1);
            }
            I_reduced(i, 0) = current[i-1];
        }
        
        // 解方程 G_reduced * U = I_reduced
        matrix U = solve_linear_system(G_reduced, I_reduced);
        
        return U(id, 0);
    }


    // TODO: 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]) {
        // 构建电压向量
        matrix U(interface_size, 1);
        for (int i = 0; i < interface_size; i++) {
            U(i+1, 0) = voltage[i];
        }
        
        // 计算每条边的电压差: V = A^T * U
        matrix AT = adjacency.transposition();
        matrix V = AT * U;
        
        // 计算功率: P = V^T * C * V = sum(C_i * V_i^2)
        fraction power(0);
        for (int i = 0; i < connection_size; i++) {
            fraction v = V(i+1, 0);
            fraction c = conduction(i+1, i);
            power = power + c * v * v;
        }
        
        return power;
    }
};


#endif //SRC_HPP