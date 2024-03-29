# Алгоритмы. Лабораторная работа №1

Работу выполнил [Максим Агеев](https://t.me/maxveega)
[Ссылка на репозиторий в GitHub](https://github.com/TobbyMax/algolab1)

**Цель**: сравнить работу трех алгоритмов поиска числа в матрице, значения в стобцах и строках которой упорядочены по возрастанию, оценить время выполнения алгоритмов при двух способах генерации данных, сделать выводы об эффективности каждого из алгоритмов.

## Алгоритм 1. Линейный поиск лесенкой.

Первый алгоритм (`linearLadderSearch`) реализует следующую логику: 
1. Начинаем поиск с правого верхнего угла матрицы (`a[0][n-1]`);
2. Пока мы не вышли за пределы матрицы (`i < m and j >= 0`), проверяем:
	a) если мы нашли искомое число, возвращаем `true`;
	b) если число, находящееся по текущему индексу больше искомого, сдвигаемся на 1 влево (уменьшаем j на единицу);
	c)  если число, находящееся по текущему индексу меньше искомого, сдвигаемся на 1 вниз (увеличиваем i на единицу);
3. Если мы вышли за пределы матрицы, значит, искомого числа в матрице нет - возвращаем `false`.

#### Код

```c++
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
```


#### Сложность алгоритма: **O(N+M)**



## Алгоритм 2. Бинарный поиск по строкам.

Второй  алгоритм (`binaryRowSearch`) реализует следующую логику: 
1. Игнорируем факт того, что значения в столбцах матрицы упорядочены;
2. Для каждой строки матрицы запускаем бинарный поиск;
3. Если после поиска во всех строках искомое число не найдено, выводим `false`.

#### Код

```c++
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
```


#### Сложность алгоритма: **O(M * log(N))**



## Алгоритм 3. Поиск лесенкой + экспоненциальный поиск.

Третий алгоритм (`exponentialLadderSearch`) реализует следующую логику: 
1. Начинаем поиск с правого верхнего угла матрицы (`a[0][n-1]`);
2. Пока мы не вышли за пределы матрицы (`i < m and j >= 0`), проверяем:
	a) если мы нашли искомое число, возвращаем `true`;
	b) если число, находящееся по текущему индексу больше искомого, запускаем экспоненциальный поиск по строке влево (`exponentialRowSearch`), после опредения границ части строки, где должно находиться искомое число, запускаем на этом отрезке бинарный поиск;
	c)  если число, находящееся по текущему индексу меньше искомого, сдвигаемся на 1 вниз (увеличиваем i на единицу);
3. Если мы вышли за пределы матрицы, значит, искомого числа в матрице нет - возвращаем `false`.

*Примечание: Исходя из того, что по техническому заданию число строк в тестовых данных не превышает числа столбцов, я использовал алгоритм, который использует линейный поиск для столбцов и  экспоненциальный поиск для строк, возиможны также модификации алгоритма, при которых используется экспоненциальный поиск для столбцов и  линейный поиск для строк, а также модификации, использующие экспоненциальный поиск и для столбцов, и для строк.*

#### Код

```c++
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
    // бинарный поиск
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
```


#### Сложность алгоритма: **O(M * (log(N) - log(M) + 1))**


## Генерация данных

### Способ 1

При первом способе генерации каждый элемент массива задается формулой `a[i][j] = ((n/m) * i + j) * 2`, а искомое число - формулой `target = 2 * n + 1`

#### Код

```c++
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
```

```c++
int** a = generate1(m, n);
int target = 2 * n + 1;
```


### Способ 2

При втором способе генерации каждый элемент массива задается формулой `a[i][j] = ((n/m) * i * j) * 2`, а искомое число - формулой `target = 16 * n + 1`

#### Код

```c++
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
```

```c++
int** a = generate2(m, n);
int target = 16 * n + 1;
```



## Запуск

Фиксируем значение N = 2^13, меняем значение M = 2^k, увеличивая k в промежутке от 0 до 13. Для каждого значения M производим 1000 тестов для каждого алгоритма на первых данных и замеряем время выполнения, затем высчитываем среднее время, чтобы добиться точности и нивелировать влияние на время выпонения таких факторов, как другие запущенные в текущий момент на ПК процессы. Результаты запусков (среднее время выполнения каждого алгоритма для конкретного M) записываем в файл results.txt. Аналогично, делаем по 1000 тестов на обоих данных для третьего алгоритма и высчитываем среднее время выполнения. Результаты запусков (среднее время выполнения алгоритма на каждом типе данных для конкретного M) записываем в файл results2.txt.
Запуск производится в IDE Сlion в компиляторе GNU 11.2.0 для C++ (стандарт 23).

#### Код

```c++
#define TESTS_COUNT 1000

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
        // Генерация данных 1 (аналогично для данных 2)  
        a = generate1(m, n);  // для данных 2: a = generate2(m, n);
        target = 2 * n + 1;  // для данных 2:  target = 16 * n + 1; 
  
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
```


## Графики

После запуска копируем полученные данные из текстового файла в файл PI1_Ageev_algolab1.xlsx  и строим по ним графики.

Графики и данные приложены в файле:
![[PI1_Ageev_algolab1.xlsx]]



## Выводы

### По тестам алгоритмов на первых данных:

1. Бинарный поиск работает быстро при небольшом количестве строк (так как по строке поиск осуществляется за O(log(N))). Так, при **M = 2^5** бинарный алгоритм выполняется за  `ConstBin * 2^5 * log( 2^13) = 13 * ConstBin * 2^5`, в то время как линейный поиск лесенкой тратит на поиск `ConstLin * (2^13 + 2^5) = ConstLin * 2^5 (2^8 + 1) = 257 * ConstLin * 2^5 > ConstLin * 2^13`.   Однако при количестве строк равном 512 бинарный поиск начинает проигрывать линейному и экпонециальному поискам примерно на *7000 наносекунд*. Это связано с тем, что бинарный алгоритм при **M = 512** требует  `ConstBin * 2^9 * log( 2^13) = 13 * ConstBin * 2^9** операций, линейному алгоритму требуется **ConstLin * (2^13 + 2^9) = ConstLin * 2^9 (2^4 + 1) = 17 * ConstLin * 2^9` операций, а экспоненциальному - `ConstExp * 2^9 * (log(2^13) - log(2^9)  + 1 ) = 5 * ConstExp * 2^9`, то есть количество операций за счет констант примерно равно, и  все три алгоритма на этих входных данных работают почти одинаково.

2. При максимальном количестве строк **M = 8192** линейный алгоритм быстрее бинарного в 14 раз, так как при **M = N = 2^13** линейный алгоритм требует `ConstLin * (2^13 + 2^13) = ConstLin * 2^14` операций, а бинарный алгоритм выполняется уже за  `ConstBin * 2^13 * log( 2^13) = 13 * ConstBin * 2^13 > ConstBin * 2^16` (если предположить, что **СonstLin = ConstBin**, то отношение количеств операций равно 13, что приблизительно равно отношению времени выполнения алгоритмов).

3. Экпоненциальный поиск лесенкой при небольшом количестве строк  (до 256) лишь немного (*до 3000 наносекунд*) проигрывает бинарному поиску по строкам, так как количество операций при **M <= 2^8** отличается, по сути,только на константу, к примеру, при **M = 2^7**, бинарный поиск выполняется за `ConstBin * 2^7 * log( 2^13) = 13 * ConstBin * 2^7`, а экспоненциальный тратит `ConstExp * 2^7 * (log(2^13) - log(2^7)  + 1 ) = 7 * ConstExp * 2^7`.

4. Однако при росте M бинарный поиск теряет свое преимущество. Как было уже отмечено выше, при **M >= 512** экспоненциальный поиск работает быстрее бинарного на *7000 наносекунд*, а при **M = 8192** экспоненциальный поиск превосходит бинарный почти в 7 раз, потому что при **M = N = 2^13** экспоненциальный алгоритм вырождается до линейного  и тратит `ConstExp * 2^13 * (log(2^13) - log(2^13)  + 1 ) = ConstExp * 2^13` операций, а бинарный алгоритм выполняется уже за  `ConstBin * 2^13 * log( 2^13) = 13 * ConstBin * 2^13 > ConstBin * 2^16`. 

5. Так как экспоненциальный алгоритм поиска при **M < 512** растет идентично бинарному, то при небольшом количестве строк экпоненциальный поиск работает быстрее линейного. При **M = 512**  алгоритмы работают с одинаковой скоростью (*t = 40000 наносекунд*). При дальнейшем увеличении M линейный поиск лесенкой показывает себя примерно в 1,5 раза лучше экспоненциального, к примеру, при M = 2^11 линейному алгоритму требуется `ConstLin * (2^13 + 2^11) = ConstLin * 2^11 (2^2 + 1) = 5 * ConstLin * 2^11` операций, а экспоненциальному - `ConstExp * 2^11 * (log(2^13) - log(2^11)  + 1 ) = 3 * ConstExp * 2^11`, а при **M = N = 2^13** линейный поиск по времени превосходит экспоненциальный в 2 раза (только из-за большей константы у экспоненциального).

6. При **M = 1**, все алгоритмы за счет оптимизации компилятора работают одинаково быстро.

7. Таким образом, бинарный поиск по строкам хорошо показывает себя в случаях, когда количество столбцов намного больше числа строк, линейный алгоритм эффективен при **N ≈ M**, а экспоненциальный поиск лесенкой работает с равной эффективностью в обоих случаях, так как при **M << N** работает почти так же, как бинарный, а при **M -> N** работает, как линейный (с поправкой на константу).


### По тестам алгоритма 3 на двух типах данных:

1.  Из-за особенности генерации данных и выбора искомого числа (цель, `target = 16 * N + 1`), алгоритм экспоненциального поиска лесенкой на вторых данных при **M <= 256** выполняет только M операций (проходит линейно по правой стенке матрице вниз), так как все элементы матрицы меньше искомого. Однако, при **M > 512** в матрице появляются элементы больше искомого и алгоритм начинает выполняться медленнее, так как начинает применять экспоненциальный поиск по строкам.

2. При генерации данных и искомого числа первым способом, в матрице уже со второй строчки появляются элементы, превосходящие искомое число, и экспоненциальный поиск по строке используется уже при **M = 2**.

4. Как следствие, алгоритм всегда работает быстрее на вторых данных, однако до **M = 256**, отношение времени выполнения растет (в пике достигает 11.5), так как на первых данных требуется `ConstExp * 2^8 * (log(2^13) - log(2^8)  + 1 ) =  6 * ConstExp * 2^8` операций, в то время как на вторых данных алгоритм завершается за **M = 256** операций, а при **M > 256** уменьшается до 2.5 раз, так как при этих значениях и на вторых данных требуется экспоненциальный поиск по строкам.
