using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gauss_Squad
{
    class Program
    {
        public static double eps = 0.000001;
        public static int sign = 1;

        /// <summary>
        /// Вывод на экран
        /// </summary>
        /// <param name="a">Исходная матрица</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="x">Матрица искомых значений</param>
        /// <param name="u">Параметр невязки</param>
        /// <param name="text">Текст для вывода на экран</param>
        public static void Print(double[,] a, double[] b, double[] x, double[] u, string text)
        {
            Console.WriteLine("\t" + text);
            Console.WriteLine();

            for (int i = 0; i < a.GetLength(0); i++)
            {
                Console.Write("|");
                for (int j = 0; j < a.GetLength(0); j++)
                {
                    Console.Write(a[i, j].ToString("0.###") + "\t");
                    if (j + 1 == a.GetLength(0)) Console.Write("|");
                }

                {
                    if (x != null) Console.Write(x[i].ToString("0.###") + "\t|");
                    if (b != null) Console.Write(b[i].ToString("0.###") + "\t|");
                    if (u != null) Console.Write(u[i].ToString("0.###") + "\t|");
                    Console.WriteLine();
                }
            }
            Console.WriteLine("---------------------------------------");
            Console.WriteLine();
        }

        /// <summary>
        /// Вычисление определителя матрицы
        /// </summary>
        /// <param name="matrix">Исходная матрица</param>
        /// <returns>Возвращает значение детерминанта матрицы</returns>
        public static double Determinant(double[,] matrix)
        {
            if (matrix.Length == 4)
            {
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            }
            if (matrix.Length == 1)
                return matrix[0, 0];
            double sign = 1, result = 0;
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                double[,] minor = GetMinor(matrix, i);
                result += sign * matrix[0, i] * Determinant(minor);
                sign = -sign;
            }
            return result;
        }

        /// <summary>
        /// Вычисление минора матрицы (определителя n-1 порядка)
        /// </summary>
        /// <param name="matrix">Матрица поиска</param>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Возвращает матрицу n-1 порядка</returns>
        private static double[,] GetMinor(double[,] matrix, int n)
        {
            double[,] result = new double[matrix.GetLength(0) - 1, matrix.GetLength(0) - 1];
            for (int i = 1; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < n; j++)
                    result[i - 1, j] = matrix[i, j];
                for (int j = n + 1; j < matrix.GetLength(0); j++)
                    result[i - 1, j - 1] = matrix[i, j];
            }
            return result;
        }

        /// <summary>
        /// Вычисление невязки решения
        /// U = b - x * A
        /// </summary>
        /// <param name="a">Исходная матрица</param>
        /// <param name="b">Матрица результатов</param>
        /// <param name="x">вектор неизвестных, полученных одним из методом</param>
        /// <returns>Возвращает матрицу невязок для каждого уравнения</returns>
        public static double[] Discrepancy(double[,] a, double[] b, double[] x, bool Sovm)
        {
            if (!Sovm)
                return null;
            double[] u = new double[a.GetLength(0)];
            for (int i = 0; i < a.GetLength(0); ++i)
            {
                double actual_b_i = 0.0;   // результат перемножения i-строки 
                // исходной матрицы на вектор x
                for (int j = 0; j < a.GetLength(0); ++j)
                    actual_b_i += a[i, j] * x[j];
                // i-й элемент вектора невязки
                u[i] = b[i] - actual_b_i;
            }
            return u;
        }

        /// <summary>
        /// Вычисление числа обусловленности
        /// </summary>
        /// <param name="first">Первая матрицы</param>
        /// <param name="second">Вторая матрица</param>
        public static void ConditionNumber(double[,] first, double[,] second)
        {
            Console.WriteLine("\tЧисло обусловленности в l-норме:");
            Console.WriteLine("\t{0:0.####}", LNorma(first) * LNorma(second));
            Console.WriteLine();
        }

        /// <summary>
        /// Вычисление l-нормы матрицы
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>
        /// <returns>Возвращает значение Lнормы для матрицы</returns>
        public static double LNorma(double[,] a)
        {
            Double norma = 0;
            for (int i = 0; i < a.GetLength(0); i++)
            {
                double sum = 0;
                double max = 0;
                for (int j = 0; j < a.GetLength(0); j++)
                    sum += Math.Abs(a[i, j]);
                if (sum > max)
                    norma = sum;
            }
            return norma;
        }

        /// <summary>
        /// Вычисление строки единичной матрицы
        /// </summary>
        /// <param name="size">Размер матрицы</param>
        /// <param name="j">Номер единичного столбца</param>
        /// <returns>Возвращает единичную матрицу</returns>
        public static double[] Matrix_E_Column(int j, int size)
        {
            double[] E = new double[size];
            for (int i = 0; i < size; i++)
                if (i == j)
                    E[i] = 1;
                else
                    E[i] = 0;
            return E;
        }

        /// <summary>
        /// Вычисление детерминанта матрицы
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>
        /// <returns>Возвращает значение найденного детерминанта</returns>
        public static double Det(double[,] a)
        {
            double determinant = 1;
            for (int i = 0; i < a.GetLength(0); i++)
                determinant *= a[i, i];
            return determinant;
        }

        /// <summary>
        /// Замена рядов в матрице
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="i">Номер первой переставляемой строки</param>
        /// <param name="n">Номер второй переставляемой строки</param>
        public static void ReplaceRows(double[,] a, double[] b, int i, int n)
        {
            double Temp;
            for (int j = 0; j < a.GetLength(0); j++)
            {
                Temp = a[i, j];
                a[i, j] = a[n, j];
                a[n, j] = Temp;
            }
            double TempB;
            for (int j = 0; j < b.Length; j++)
            {
                TempB = b[i];
                b[i] = b[n];
                b[n] = TempB;
            }
        }

        /// <summary>
        /// Замена столбцов в матрице
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="index">Матрица индексов столбцов</param>
        /// <param name="i">Номер первого переставляемого столбца</param>
        /// <param name="n">Номер второго переставляемого столбца</param>
        public static void ReplaceColumns(double[,] a, double[] b, int[] index, int i, int n)
        {
            double Temp;
            for (int j = 0; j < a.GetLength(0); j++)
            {
                Temp = a[j, i];
                a[j, i] = a[j, n];
                a[j, n] = Temp;
            }

            //фиксируем перестановку столбцов 
            int Elem = index[i];
            index[i] = index[n];
            index[n] = Elem;
        }


        /// <summary>
        /// Поиск главного элемента в матрице по строкам
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>        
        /// <param name="b">Матрица результата</param>
        /// <param name="index">Матрица индексов столбцов</param>
        /// <param name="eps">Заданная точность вычисления</param>
        /// <param name="row">Номер строки</param>
        /// <param name="column">Номер столбца</param>
        /// <returns>Возвращает значение найденного главного элемента</returns>
        private static double FindR(double[,] a, double[] b, int[] index, double eps, int row, int column)
        {
            int max_i = row;
            int max_j = column;
            double max = a[row, column];
            double max_abs = Math.Abs(max);
            for (int cur_i = row; cur_i < a.GetLength(0); cur_i++)
            {
                double cur = a[cur_i, column];
                double cur_abs = Math.Abs(cur);
                if (cur_abs > max_abs)
                {
                    max_i = cur_i;
                    //max_j = cur_j;
                    max = cur;
                    max_abs = cur_abs;
                }
            }

            if (max_i != row)
            {
                ReplaceRows(a, b, row, max_i);
                sign *= -1;
            }
            else
            {
                if (max_abs < eps)
                {
                    FindC(a, b, index, row, column);
                }
            }
            return max;
        }

        /// <summary>
        /// Поиск главного элемента в матрице по столбцам
        /// </summary>
        /// <param name="a">Исходная матрица поиска</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="index">Матрица индексов столбцов</param>
        /// <param name="row">Номер строки</param>
        /// <param name="column">Номер столбца</param>
        private static void FindC(double[,] a, double[] b, int[] index, int row, int column)
        {
            int max_i = row;
            int max_j = column;
            double max = a[row, column];
            double max_abs = Math.Abs(max);
            for (int cur_j = column; cur_j < a.GetLength(0); cur_j++)
            {
                double cur = a[row, cur_j];
                double cur_abs = Math.Abs(cur);
                if (cur_abs > max_abs)
                {
                    max_j = cur_j;
                    max = cur;
                    max_abs = cur_abs;
                }
            }
            //замена в случае необходимости строк
            if (max_j != column)
            {
                ReplaceColumns(a, b, index, row, max_j);
            }
        }

        /// <summary>
        /// Расчет СЛАУ методом Гаусса
        /// </summary>
        /// <param name="a">Исходная матрица</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="index">Матрица индексов столбцов</param>
        /// <param name="Invers">Параметр обратной матрицы</param>
        /// <param name="Sovm">Параметр совместной матрицы</param>
        /// <param name="det">Значение детерминанта</param>
        /// <returns>Возвращает матрицу искомых значений</returns>
        private static double[] CalculateGauss(double[,] a, double[] b, int[] index, bool Invers, bool Sovm, double det)
        {
            // Прямой ход метода Гаусса
            // перемещаемся по каждой строке сверху вниз
            int size = a.GetLength(0);
            //int sign = 1;
            double[] x = new double[size];
            double r = 0;
            for (int i = 0; i < size; ++i)//k
            {
                // 1) выбор главного элемента
                r = FindR(a, b, index, eps, i, i);
                double p;
                // 2) Вычитание текущей строки из всех нижерасположенных строк
                if (a[i, i] == 0)
                    p = 0;
                else
                {
                    for (int k = i + 1; k < size; ++k)
                    {
                        p = a[k, i] / a[i, i];
                        for (int j = i; j < size; ++j)
                            a[k, j] -= a[i, j] * p;
                        b[k] -= b[i] * p;
                    }
                }
            }

            if (!Invers)
            {
                Invers = false;
                int[] count = new int[size];
                for (int i = size - 1; i >= 0; i--)
                {
                    for (int j = size - 1; j >= 0; j--)
                    {
                        if (Math.Abs(a[i, j]) < eps)
                        {
                            count[i]++;
                        }
                    }
                    if (count[i] == size && Math.Abs(b[i]) > eps)
                    {
                        Console.WriteLine("Система несовместна");
                        Sovm = false;
                    }
                    if (count[i] == size && Math.Abs(b[i]) <= eps)
                    {
                        Invers = true;
                        break;
                    }
                }
                //Вывод треугольной матрицы 
                Print(a, b, null, null, "Треугольная матрица:");
            }
            if (!Sovm)
                return null;
            // Обратный ход метода Гаусса
            x[size - 1] = b[size - 1] / a[size - 1, size - 1];
            // перемещаемся по каждой строке снизу вверх
            for (int i = size - 2; i >= 0; --i)
            {
                double dSum = 0.0;
                for (int j = size - 1; j > i; j--)
                    dSum += a[i, j] * x[j];
                if (Math.Abs(b[i] - dSum) < eps)
                    x[i] = 0;
                else
                    x[i] = (b[i] - dSum) / a[i, i];
            }
            //восстановление порядка переменных 
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (i == index[j])
                    {
                        double Temp;
                        Temp = x[j];
                        x[j] = x[i];
                        x[i] = Temp;
                    }
            return x;
        }

        /// <summary>
        /// Нахождение обратной матрицы
        /// </summary>
        /// <param name="a">Исходная матрица</param>
        /// <param name="index">Матрица индексов столбцов</param>
        /// <param name="det">Значение детерминанта</param>
        /// <returns>Возвращает обратную матрицу</returns>
        public static double[,] Inverse(double[,] a, int[] index, double det)
        {
            if (Math.Abs(det) < eps)
            {
                return null;
            }
            int size = a.GetLength(0);
            double[,] E = new double[size, size];//  еденичная матрица
            double[,] inverse = new double[size, size];//  еденичная матрица
            double[,] A_Copy = new double[size, size];
            double[] x_ = new double[size]; // стобец обратной матрицы
            bool Invers = true;
            for (int k = 0; k < size; k++)
            {
                //находим методом Гаусса, вектор b - столбец матрицы E
                A_Copy = (double[,])a.Clone();
                x_ = CalculateGauss(A_Copy, Matrix_E_Column(k, size), index, Invers, true, det);
                for (int i = 0; i < size; i++)
                    inverse[i, k] = x_[i];
            }
            return inverse;
        }

        /// <summary>
        /// Расчет матрицы Гильберта
        /// </summary>
        /// <param name="H">Гильбертова матрица</param>
        /// <param name="det">Значение детерминанта</param>
        /// <param name="n">Размер матрицы</param>
        public static void Hilbert(double[,] H, double det, int n)
        {
            int[] index = new int[n];
            for (int i = 0; i < index.Length; i++)
                index[i] = i;
            det = Determinant(H);
            double[,] E = new double[n, n];//  еденичная матрица
            double[,] A_Copy = new double[n, n];
            double[,] H_Inverse = new double[n, n];
            double[] x_ = new double[n]; // стобец обратной матрицы
            bool Invers = true;
            for (int k = 0; k < n; k++)
            {
                //находим методом Гаусса, вектор b - столбец матрицы E
                A_Copy = (double[,])H.Clone();
                x_ = CalculateGauss(A_Copy, Matrix_E_Column(k, n), index, Invers, true, det);
                for (int i = 0; i < n; i++)
                    H_Inverse[i, k] = x_[i];
            }
            Print(H, null, null, null, "Гильбертова матрица");
            ConditionNumber(H, H_Inverse);

        }

        /// <summary>
        /// Проверка на семметричность
        /// </summary>
        /// <param name="a">Исходный массив поиска</param>
        /// <returns>Возвращает true если симметрична и false если нет</returns>
        private static bool Symmetry(double[,] a)
        {
            bool symm = true;
            for (int i = 0; i < a.GetLength(0); i++)
            {
                for (int j = 0; j < a.GetLength(0); j++)
                    if (a[i, j] != a[j, i])
                    {
                        symm = false;
                        break;
                    }
                if (!symm) break;
            }
            return symm;
        }

        public static bool sg(double x)
        {
            return x > 0;
        }

        /// <summary>
        /// Расчет СЛАУ методом квадратного корня
        /// </summary>
        /// <param name="a">Исходная матрица</param>
        /// <param name="b">Матрица результата</param>
        /// <param name="Sovm">Параметр совместной матрицы</param>
        /// <returns>Возвращает матрицу искомых значений</returns>
        public static double[] CalculateSquareRoot(double[,] a, double[] b, bool Sovm)
        {
            if (!Sovm)
                return null;
            int size = a.GetLength(0);
            double[,] S = new double[size, size];//Верняя треугольная матрица
            double[,] St = new double[size, size];//Нижняя треугольная матрица
            double[,] D = new double[size, size];
            double[] y = new double[size]; //Вектор промежуточного решения
            double[] z = new double[size]; //Вектор промежуточного решения
            double[] x = new double[size];

            S[0, 0] = Math.Sqrt(Math.Abs(a[0, 0]));//Задаем первый элемент матрицы L как корень из первого элемента матрицы А
            if (sg(a[0, 0]))
                D[0, 0] = 1;
            else
                D[0, 0] = -1;
            for (int i = 1; i < size; i++)
            {
                for (int j = 1; j < size; j++)
                {
                    if (i - 1 < j)
                    {
                        double SumJ_El = 0;
                        for (int k = 0; k < i; k++)
                            SumJ_El += S[k, i - 1] * S[k, j] * D[k, k];
                        S[i - 1, j] = (a[i - 1, j] - SumJ_El) / (S[i - 1, i - 1] * D[i - 1, i - 1]);
                    }

                    if (i == j)
                    {
                        double SumI_El = 0;
                        for (int k = 0; k < i; k++)
                        {
                            SumI_El += (S[k, i] * S[k, i]) * D[k, k];
                        }
                        if (sg(a[i, i] - SumI_El))
                            D[i, i] = 1;
                        else
                            D[i, i] = -1;
                        S[i, i] = Math.Sqrt(Math.Abs(a[i, i] - SumI_El));
                        //if (S[i,i] == 0)
                    }
                }

            }
            Print(S, null, null, null, "Верхняя треугольная матрица");

            //Транспонируем матрицу L и получаем матрицу LT
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    St[i, j] = S[j, i];// По сути транспонирование - записываем рядки, как столбцы
            Print(St, null, null, null, "Нижняя треугольная матрица");

            //Рассчитываем вектор у
            double summa = 0;
            for (int i = 0; i < size; i++)
            {
                summa = 0;
                for (int j = 0; j < i; j++)
                {
                    summa += (S[j, i] * z[j]);
                }
                z[i] = (b[i] - summa) / S[i, i];
            }
            for (int i = 0; i < size; i++)
                y[i] = z[i] / D[i, i];


            //Находим значение наших иксов
            for (int i = size - 1; i >= 0; i--)
            {
                summa = 0;
                for (int j = size - 1; j > i; j--)
                {
                    summa += S[i, j] * x[j];
                }
                x[i] = (y[i] - summa) / S[i, i];
            }
            //Print(L, x, null, null, "");
            return x;
        }

        /// <summary>
        /// Инициализация матрицы
        /// </summary>
        /// <param name="allStr">Массив введенных данных</param>
        /// <param name="a_matrix">Исходная матрица</param>
        /// <param name="vector">Матрица результата</param>
        /// <param name="size">Размерность матрицы</param>
        static void Initialization(string[] allStr, double[,] a_matrix, double[] vector, int size)
        {
            String[] str; // одна строка из введённых данных

            for (int i = 1; i <= size; i++) // в остальных координаты для фигуры
            {
                str = allStr[i].Split(new String[] { " " }, StringSplitOptions.RemoveEmptyEntries);
                for (int j = 0; j < size; j++)
                    a_matrix[i - 1, j] = Convert.ToDouble(str[j]);
            }
            for (int i = size + 1; i <= size * 2; i++) // в остальных координаты для фигуры
            {
                vector[i - (size + 1)] = Convert.ToDouble(allStr[i]);
            }
        }

        /// <summary>
        /// Решение СЛАУ методом Гаусса
        /// </summary>
        static void Gauss()
        {
            String path = "input.txt";
            string[] allStr = File.ReadAllLines(path); // прочитываем из файла все строки

            double[,] a_matrix;          // матрица A
            double[] solution;           // вектор неизвестных x
            double[] vector;           // вектор результата
            double[] discrepancy;           // вектор невязки U
            double[,] inverse;         // обратная матрица
            int[] index_gauss;
            bool Sovm;
            int size;                  // размерность задачи
            double det;
            bool Invers = false;


            size = Convert.ToInt32(allStr[0]); // в первой строчке размер
            a_matrix = new double[size, size];
            vector = new double[size];
            Initialization(allStr, a_matrix, vector, size);

            solution = new double[size];
            discrepancy = new double[size];
            inverse = new double[size, size];//конечная обратная матрица
            det = 1;
            Sovm = true;
            double[,] work_a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            double[] work_b = (double[])vector.Clone();  // с копией будем производить вычисления

            index_gauss = new int[size];
            for (int i = 0; i < index_gauss.Length; i++)
                index_gauss[i] = i;


            Print(work_a, work_b, null, null, "Исходная система: ");
            solution = CalculateGauss(work_a, work_b, index_gauss, Invers, Sovm, det);
            det = sign * Det(work_a);
            double[,] a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            double[] b = (double[])vector.Clone();  // с копией будем производить вычисления
            double[] x = (double[])solution.Clone();  // с копией будем производить вычисления
            discrepancy = Discrepancy(a, b, x, Sovm);
            Console.WriteLine("\tДетерминант {0}", det);
            Console.WriteLine();
            inverse = Inverse(a, index_gauss, det);
            Print(inverse, null, null, null, "Обратная матрица:");
            double[,] _a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            Print(_a, vector, solution, discrepancy, "Решение системы: ");
            ConditionNumber(_a, inverse);
        }

        /// <summary>
        /// Решение СЛАУ методом квадратного корня
        /// </summary>
        static void SqRoot()
        {
            String path = "input2.txt";
            string[] allStr = File.ReadAllLines(path); // прочитываем из файла все строки

            double[,] a_matrix;          // матрица A
            double[] solution;           // вектор неизвестных x
            double[] vector;           // вектор результата
            double[] discrepancy;           // вектор невязки U
            double[,] inverse;         // обратная матрица
            int[] index_root;
            bool Sovm;
            int size;                  // размерность задачи
            double det;

            size = Convert.ToInt32(allStr[0]); // в первой строчке размер
            a_matrix = new double[size, size];
            vector = new double[size];
            Initialization(allStr, a_matrix, vector, size);

            solution = new double[size];
            discrepancy = new double[size];
            inverse = new double[size, size]; //конечная обратная матрица
            det = 1;
            Sovm = true;
            double[,] work_a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            double[] work_b = (double[])vector.Clone();  // с копией будем производить вычисления

            index_root = new int[size];
            for (int i = 0; i < index_root.Length; i++)
                index_root[i] = i;


            Print(work_a, work_b, null, null, "Исходная система: ");
            solution = CalculateSquareRoot(work_a, work_b, Sovm);

            double[,] a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            double[] b = (double[])vector.Clone();  // с копией будем производить вычисления
            double[] x = (double[])solution.Clone();  // с копией будем производить вычисления

            double[,] _a = (double[,])a_matrix.Clone(); // с копией будем производить вычисления
            Print(_a, vector, solution, discrepancy, "Решение системы: ");
        }

        /// <summary>
        /// Нахождение матрицы гильберта
        /// </summary>
        /// <param name="n">Размер матрицы</param>

        static void Hil(int n)
        {
            double[,] H = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    H[i, j] = (1f / ((i + 1) + (j + 1) - 1));

            Hilbert(H, 1, n);
        }

        /// <summary>
        /// Точка входа в программу
        /// Вызов методов решения СЛАУ 
        /// Гаусс
        /// Квадратный корень
        /// Нахождение матрицы Гильберта
        /// </summary>
        static void Main(string[] args)
        {
            //
            //---------------------------------------------------------------------------------------------------------------------------------------------------
            //
            Console.WriteLine("\tМетод Гаусса");
            Console.WriteLine();
            Gauss();

            //
            //------------------------------------------------------------------------------------------------------------------------------------------------------
            //
            Console.WriteLine("\tМетод Квадратного корня");
            Console.WriteLine();
            SqRoot();

            //
            //-------------------------------------------------------------------------------------------------------------------------------------------------------
            //
            Console.WriteLine("\tМатрица Гильберта");
            Console.WriteLine();
            Hil(5);

            //
            //---------------------------------------------------------------------------------------------------------------------------------------------------------
            //
            Console.ReadKey();
        }
    }
}