#include <iostream>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <string.h>
#include "matrix_mult.h"
#include "SortingAlgorithms.h"
using namespace std;

// bubblesort

void bubbleSort(int arr[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)
    {
        // Últimos i elementos ya están en su lugar correcto
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j] > arr[j + 1])
            {
                // Intercambia arr[j] y arr[j+1]
                int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

// mergesort

void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    n1 = 0;
    int n2 = r - m;

    // Crea dos arreglos temporales para guardar los elementos
    int L[n1], R[n2];

    // Copia los datos a los arreglos temporales L[] y R[]
    for (i = 0; i < n1; i++)
    {
        L[i] = arr[l + i];
    }
    for (j = 0; j < n2; j++)
    {
        R[j] = arr[m + 1 + j];
    }

    // Combina los arreglos temporales de vuelta en arr[l..r]
    i = 0; // Índice inicial del subarreglo izquierdo
    j = 0; // Índice inicial del subarreglo derecho
    k = l; // Índice inicial del subarreglo mezclado
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    // Copia los elementos restantes de L[], si los hay
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    // Copia los elementos restantes de R[], si los hay
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(int arr[], int l, int r)
{
    if (l < r)
    {
        // Encuentra el punto medio
        int m = l + (r - l) / 2;

        // Ordena la primera y segunda mitad del arreglo
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        // Mezcla las dos mitades ordenadas
        merge(arr, l, m, r);
    }
}

// quicksort

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(int arr[], int low, int high)
{
    int pivot = arr[high]; // Toma el último elemento como pivote
    int i = (low - 1);     // Índice del elemento más pequeño

    for (int j = low; j <= high - 1; j++)
    {
        // Si el elemento actual es menor o igual que el pivote
        if (arr[j] <= pivot)
        {
            i++;                    // Incrementa el índice del elemento más pequeño
            swap(&arr[i], &arr[j]); // Intercambia el elemento actual con el elemento más pequeño
        }
    }
    swap(&arr[i + 1], &arr[high]); // Intercambia el pivote con el elemento más pequeño
    return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        // Encuentra el índice de partición, arr[p] es ahora el pivote
        int p = partition(arr, low, high);

        // Ordena los elementos antes y después de la partición
        quickSort(arr, low, p - 1);
        quickSort(arr, p + 1, high);
    }
}

// multiplicacion de matrices iterativo cubico

void multiplyMatrix(int mat1[][MAX], int mat2[][MAX], int res[][MAX], int n)
{
    // Inicializa la matriz de resultado con 0
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = 0;
        }
    }

    // Multiplica las matrices mat1 y mat2 y almacena el resultado en res
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

// multiplicacion de matrices con algoritmo iterativo cúbico optimizado para mantener la localidad de los datos

void matrixMultiplication(int A[][MAX], int B[][MAX], int C[][MAX], int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            int sum = 0;
            for (k = 0; k < n; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

// multiplicacion de matrices con algoritmo de Strassen

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

int main(int argv, char *argc[])
{
    srand(time(NULL));

    int n, m, k, N, n_1, i, num_of_experiments;
    string selected_algorithm;

    switch (atoi(argc[1]))
    //1, 2, 3 = algoritmos de ordenamiento con dataset tipo 1, 2 o 3
    //4, 5, 6 multiplicacion de matrices con dataset aleatorios 
    {
    case 2:
        selected_algorithm = "transpose_mm";
        break;
    case 3:
        selected_algorithm = "optimized_mm";
        break;
    case 4:
        selected_algorithm = "optimized_mm";
        break;
    case 5:
        selected_algorithm = "optimized_mm";
        break;
    case 6:
        selected_algorithm = "optimized_mm";
        break;
    default:
        selected_algorithm = "standard_mm";
        break;
    }

    if ((argv > 2) && (strcmp(argc[2], "--test") == 0))
    {
        cin >> n >> m >> k;
        vector<vector<int>> M_A(n, vector<int>(m, 0)), M_B(m, vector<int>(k, 0));
        read_matrix(M_A);
        read_matrix(M_B);

        vector<vector<int>> result;
        result = matrix_multiplication(M_A, M_B, selected_algorithm);
        print_matrix(result);
        return 0;
    }

    string outfile_name = selected_algorithm + "_results.csv";
    ofstream outfile(outfile_name);
    string column_names = "n,time[ms]\n";
    outfile << column_names;

    i = 100, n_1 = 1, N = 1000, num_of_experiments = 10;

    for (int n = n_1; n <= N; n += i)
    {
        cout << n << endl;
        double mm_total_time = 0;
        vector<vector<int>> M_A(n, vector<int>(n, 0)), M_B(n, vector<int>(n, 0));

        for (int j = 0; j < num_of_experiments; j++)
        {
            long long single_execution_time = execution_time_ms(M_A, M_B, selected_algorithm);
            mm_total_time += single_execution_time;
        }
        double mm_avg_time = mm_total_time / num_of_experiments;
        outfile << n << "," << mm_avg_time << endl;
    }
    outfile.close();
    return 0;
}