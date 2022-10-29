#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>

using namespace std;

#define TESTS_COUNT 1000

int** generate1(int m, int n) {
	int** a = new int* [m];

	for (size_t i = 0; i < m; ++i) {
		a[i] = new int[n];
	}

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] = ((n/m) * i + j) * 2;
		}
	}

	return a;
}

int** generate2(int m, int n) {
    int** a = new int* [m];

    for (size_t i = 0; i < m; ++i) {
        a[i] = new int[n];
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = ((n/m) * i * j) * 2;
        }
    }

    return a;
}

void print(int** a, int m, int n) {
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			std::cout << a[i][j] << "  ";
		}
		std::cout << "\n";
	}
}

void deleteMatrix(int** a, int m) {
    for (size_t i = 0; i < m; ++i) {
        delete[] a[i];
    }
    delete[] a;
}

bool linearLadderSearch(int** a, int m, int n, int target) {
	int i = 0, j = n-1;
	while (i < m and j >= 0) {
		if (a[i][j] == target) {
			return true;
		}
		else if (a[i][j] > target) {
			j--;
		}
		else {
			i++;
		}
	}
    return false;
}

bool binaryRowSearch(int** a, int m, int n, int target) {
    for (int i = 0; i < m; ++i) {
        int l = 0;
        int r = n - 1;
        int mid = 0;
        while (l < r) {
            mid = (l + r) / 2;
            if (a[i][mid] >= target) {
                r = mid;
            } else {
                l = mid + 1;
            }
        }
        if (a[i][l] == target) {
            return true;
        }
    }
    return false;
}

int exponentialColumnSearch(int** a, int m, int i, int j,  int target) {
    int delta = 1, previousI = i;
    while (i < m and a[i][j] < target) {
        previousI = i;
        i += delta;
        delta *= 2;
    }
    if (i >= m) {
        i = m - 1;
    }
    int l = previousI;
    int r = i;
    int mid;
    while (l < r) {
        mid = (l + r) / 2;
        if (a[mid][j] >= target) {
            r = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}

bool exponentialLadderSearch2(int** a, int m, int n, int target) {
    int i = 0, j = n - 1;
    while (i < m and j >= 0) {
        if (a[i][j] == target) {
            return true;
        }
        else if (a[i][j] > target) {
            j--;
        }
        else {
            if (i == m - 1) {
                return false;
            }
            i = exponentialColumnSearch(a, m, i, j, target);
        }
    }
    return false;
}

int exponentialRowSearch(int** a, int n, int i, int j,  int target) {
    int delta = 1, previousJ = j;
    while (j >= 0 and a[i][j] > target) {
        previousJ = j;
        j -= delta;
        delta *= 2;
    }
    if (j < 0) {
        j = 0;
    }
    int l = j;
    int r = previousJ;
    int mid;
    while (l < r) {
        mid = (l + r) / 2;
        if (a[i][mid] <= target) {
            r = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}


bool exponentialLadderSearch(int** a, int m, int n, int target) {
    int i = 0, j = n - 1;
    while (i < m and j >= 0) {
        if (a[i][j] == target) {
            return true;
        }
        else if (a[i][j] > target) {
            if (j == 0) {
                return false;
            }
            j = exponentialRowSearch(a, n, i, j, target);
        }
        else {
            i++;
        }
    }
    return false;
}

int main()
{
    int target;
    int** a;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();

    ofstream fout("C:\\Users\\ageev\\CLionProjects\\algolab1\\results.txt");
    fout << "m\tLinearLadderSearch\tBinaryRowSearch\tExponentialLadderSearch\n";
    // Фиксируем N = 2^13
    int n = pow(2, 13);
    // Тестируем 3 алгоритма на данных 1
    // Меняем M в промежутке 2^k, k принадлежит [0, 13]
    for (int m = 1; m <= n; m *= 2) {
        cout << "Size: " << m << endl;
        // Генерация данных 1
        a = generate1(m, n);
        target = 2 * n + 1;

        long long sumLin = 0, sumBin = 0, sumExp = 0;

        // Делаем 1000 запусков для точности значений времени
        for (int i = 0; i < TESTS_COUNT; ++i) {
            // холостой ход, чтобы нивелировать погрешности при первом запуске
            linearLadderSearch(a, m, n, target);
            binaryRowSearch(a, m, n, target);
            exponentialLadderSearch(a, m, n, target);

            t0 = std::chrono::high_resolution_clock::now();
            linearLadderSearch(a, m, n, target);
            t1 = std::chrono::high_resolution_clock::now();
            sumLin += chrono::duration_cast<chrono::nanoseconds >(t1-t0).count();

            t0 = std::chrono::high_resolution_clock::now();
            binaryRowSearch(a, m, n, target);
            t1 = std::chrono::high_resolution_clock::now();
            sumBin += chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();

            t0 = std::chrono::high_resolution_clock::now();
            exponentialLadderSearch(a, m, n, target);
            t1 = std::chrono::high_resolution_clock::now();
            sumExp += chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();
        }
        // Вычисляем среднее время за 1000 запусков
        long long timeLin = sumLin / TESTS_COUNT, timeBin = sumBin / TESTS_COUNT, timeExp = sumExp / TESTS_COUNT;
        // Записываем данные в текстовый файл
        fout << m << "\t" << timeLin << "\t" << timeBin << "\t" << timeExp << endl;

        deleteMatrix(a, m);
    }

    ofstream fout2("C:\\Users\\ageev\\CLionProjects\\algolab1\\results2.txt");
    fout2 << "m\tData1\tData2\n";
    // Тестируем алгоритм O(log(n) - log (n) + 1) на двух данных
    // Меняем M в промежутке 2^k, k принадлежит [0, 13]
    for (int m = 1; m <= n; m *= 2) {
        cout << "Size: " << m << endl;
        // Генерация данных 1
        int** a1 = generate1(m, n);
        int target1 = 2 * n + 1;
        // Генерация данных 2
        int** a2 = generate2(m, n);
        int target2 = 16 * n + 1;

        long long sum1 = 0, sum2 = 0;
        // Делаем 1000 запусков для точности значений времени
        for (int i = 0; i < TESTS_COUNT; ++i) {
            // холостой ход, чтобы нивелировать погрешности при первом запуске
            exponentialLadderSearch(a1, m, n, target1);
            exponentialLadderSearch(a1, m, n, target1);

            auto t0 = std::chrono::high_resolution_clock::now();
            exponentialLadderSearch(a1, m, n, target1);
            auto t1 = std::chrono::high_resolution_clock::now();
            sum1 += chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

            // холостой ход, чтобы нивелировать погрешности при первом запуске
            exponentialLadderSearch(a2, m, n, target2);
            exponentialLadderSearch(a2, m, n, target2);

            t0 = std::chrono::high_resolution_clock::now();
            exponentialLadderSearch(a2, m, n, target2);
            t1 = std::chrono::high_resolution_clock::now();
            sum2 += chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        }
        // Вычисляем среднее время за 1000 запусков
        long long time1 = sum1 / TESTS_COUNT, time2 = sum2 / TESTS_COUNT;
        // Записываем данные в текстовый файл
        fout2 << m << "\t" << time1 << "\t" << time2 << endl;

        deleteMatrix(a1, m);
        deleteMatrix(a2, m);
    }
}
