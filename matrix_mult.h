#include <iostream>
#include <vector>
#include <functional>

using namespace std;

template <typename Func>
long long execution_time_ms(Func function, const vector<vector<int>> &A, const vector<vector<int>> &B, string alg)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  function(A, B, alg);
  auto end_time = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
}

void read_matrix(vector<vector<int>> &M)
{
  for (auto &row : M)
    for (auto &element : row)
      cin >> element;
}

void print_matrix(const vector<vector<int>> &M)
{
  for (int i = 0; i < M.size(); i++)
  {
    for (int j = 0; j < M[i].size(); j++)
      cout << M[i][j] << " ";
    cout << endl;
  }
}

vector<vector<int>> standard_mm(const vector<vector<int>> &A, const vector<vector<int>> &B)
{
  int n = A.size();
  int m = A[0].size();
  int k = B[0].size();

  vector<vector<int>> C(n, vector<int>(k, 0));

  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++)
      for (int k = 0; k < m; k++)
        C[i][j] += A[i][k] * B[k][j];

  return C;
}

vector<vector<int>> transpose_mm(const vector<vector<int>> &A, const vector<vector<int>> &B)
{
  
}

void strassen(int A[][MAX], int B[][MAX], int C[][MAX], int n)
{
  // Caso base para matrices de 1x1
  if (n == 1)
  {
    C[0][0] = A[0][0] * B[0][0];
    return;
  }

  // Particionado de las matrices en bloques
  int m = n / 2;
  int A11[MAX][MAX], A12[MAX][MAX], A21[MAX][MAX], A22[MAX][MAX];
  int B11[MAX][MAX], B12[MAX][MAX], B21[MAX][MAX], B22[MAX][MAX];
  int C11[MAX][MAX], C12[MAX][MAX], C21[MAX][MAX], C22[MAX][MAX];
  int P1[MAX][MAX], P2[MAX][MAX], P3[MAX][MAX], P4[MAX][MAX], P5[MAX][MAX], P6[MAX][MAX], P7[MAX][MAX];

  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      A11[i][j] = A[i][j];
      A12[i][j] = A[i][j + m];
      A21[i][j] = A[i + m][j];
      A22[i][j] = A[i + m][j + m];

      B11[i][j] = B[i][j];
      B12[i][j] = B[i][j + m];
      B21[i][j] = B[i + m][j];
      B22[i][j] = B[i + m][j + m];
    }
  }

  // Cálculo de los productos de Strassen
  strassen(sum(A11, A22, m), sum(B11, B22, m), P1, m);
  strassen(sum(A21, A22, m), B11, P2, m);
  strassen(A11, sub(B12, B22, m), P3, m);
  strassen(A22, sub(B21, B11, m), P4, m);
  strassen(sum(A11, A12, m), B22, P5, m);
  strassen(sub(A21, A11, m), sum(B11, B12, m), P6, m);
  strassen(sub(A12, A22, m), sum(B21, B22, m), P7, m);

  // Cálculo de los bloques de la matriz resultante
  sum(sub(sum(P1, P4, m), P5, m), P7, m, C11);
  sum(P3, P5, m, C12);
  sum(P2, P4, m, C21);
  sum(sub(sum(P1, P3, m), P2, m), P6, m, C22);

  // Copia de los bloques a la matriz resultante
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      C[i][j] = C11[i][j];
      C[i][j + m] = C12[i][j];
      C[i + m][j] = C21[i][j];
      C[i + m][j + m] = C22[i][j];
    }
  }
}

vector<vector<int>> matrix_multiplication(const vector<vector<int>> &A, const vector<vector<int>> &B, string alg)
{
  if (alg == "transpose_mm")
    return transpose_mm(A, B);
  return standard_mm(A, B);
}
