using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
 
namespace Jakobi_Seidel
{
    class Program
    {
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

        // Entry point
        // Точка входа.
        static void Main(string[] args)
        {
            String path = "input.txt";
            string[] allStr = File.ReadAllLines(path); // прочитываем из файла все строки

            int n = Convert.ToInt32(allStr[0]); // в первой строчке размер
            // Matrix coefficients
            // Матрица коеффциентов СЛАУ
            double[,] matrix = new double[n, n];


            // matrix of free coefficients
            // матрица свободных членов
            double[] additional = new double[n];

            Initialization(allStr, matrix, vector: additional, size: n);
            showMatrix(n, matrix, additional, null, "Исходная система");

            //if (conv_condition(matrix, n))
            {
                // Create and init.
                // Объявляем и инициализимруем классы.
                Seidel S = new Seidel(n, matrix, additional, 0.001);
                Jacobi J = new Jacobi(n, matrix, additional, 0.001);

                //S.calculateMatrix();
                //J.calculateMatrix();

                // set method args of ThreadStart delegate 
                // Передаем методы потоку через делегат ThreadStart
                Thread Z = new Thread(new ThreadStart(S.calculateMatrix));
                Thread Y = new Thread(new ThreadStart(J.calculateMatrix));

                // Start threads
                // Запускаем потоки.
                Z.Start();
                Y.Start();

                // wait for endings
                // Ожидаем завершения.
                Z.Join();
                Y.Join();

                // Show results of calculations 
                // Выводим на экран.
                showMatrix(n, matrix, additional, S.ResultMatrix, "Метод Зейделя");
                show(S.Iteration);
                showMatrix(n, matrix, additional, J.ResultMatrix, "Метод Якоби");
                show(J.Iteration);
            }
            Console.ReadKey();
        }


        //------------------------------------------
        // showMatrix - method overloads showing results 
        // showMatrix - перегрузки разнотипных выводов
        //------------------------------------------
        public static void showMatrix(int n, double[,] a, double[] b, double[] x, string text)
        {
            Console.WriteLine("\t" + text);
            Console.WriteLine();

            for (int i = 0; i < n; i++)
            {
                Console.Write("|");
                for (int j = 0; j < n; j++)
                {
                    Console.Write(a[i, j].ToString("0.###") + "\t");
                    if (j + 1 == a.GetLength(0)) Console.Write("|");
                }

                {
                    if (x != null) Console.Write(x[i].ToString("0.########") + "\t|");
                    if (b != null) Console.Write(b[i].ToString("0.###") + "\t|");
                    Console.WriteLine();
                }
            }
            Console.WriteLine("---------------------------------------");
            Console.WriteLine();
        }
        static void showMatrix(double[,] x)
        {
            Console.WriteLine("\n Result:");

            for (int i = 0; i < x.GetLength(0); i++)
            {
                for (int j = 0; j < x.GetLength(1); j++)
                {
                    Console.Write(" {0} ", x[i, j]);
                }
                Console.WriteLine();
            }
        }

        static void show(int it)
        {
            Console.WriteLine("\t|Итераций: " + it.ToString("0.###") + "\t|");
            Console.WriteLine("---------------------------------------");
            Console.WriteLine();
        }
        static void showMatrix(double[] x)
        {
            Console.WriteLine("\n Result:");
            for (int i = 0; i < x.Length; i++)
            {
                Console.WriteLine(" {0} ", x[i]);
            }
        }

        //достаточное условие сходимости
        static bool conv_condition(double[,] A, int n)
        {
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    if (i != j)
                        sum += Math.Abs(A[i,j]);
                sum /= Math.Abs(A[i,i]);
                if (sum > 1)
                {
                    Console.WriteLine("Достаточное условие сходимости не выполнено");
                    return false;
                }
            }
            return true;
        }

        public static double[] Column_Mul(int n,double[,] A, double[] B)
        {
            double[] res = new double[n];
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    res[row] += A[row, col] * B[col];
                }
            }
            return res;
        }

        public static double[] Discrepancy(int n, double[] b, double[] x)
        {
            //if (!Sovm)
            //    return null;
            double[] u = new double[n];
            for (int i = 0; i < n; ++i)
                // i-й элемент вектора невязки
                u[i] = b[i] - x[i];
            return u;
        }
    }

    // Base class of all iteration methods
    // Общий класс для итерационных методов.
    abstract class SimpleIterations
    {
        public abstract void calculateMatrix();
    }

    /// <summary>
    /// Class Jacobi 
    /// Класс отвечает за работу метода Якоби
    /// </summary>
    class Jacobi : SimpleIterations
    {
        // Матрица ответов
        private double[] resultMatrix;
        public double[] ResultMatrix
        {
            get
            {
                if (resultMatrix != null)
                    return resultMatrix;
                else
                {
                    return new double[n];
                }
            }
        }

        // Основная матрица и свободные члены.
        private double[,] matrix;
        private double[] addtional;
        private int n;
        private int count_iteration;

        // точность (кол-во итераций)
        private double accuracy;
        // избегаем ошибок с итерациями.
        public double Accuracy
        {
            get
            {
                return accuracy;
            }
            set
            {
                if (value <= 0.0)
                    accuracy = 0.1;
                else
                    accuracy = value;
            }
        }

        // Конструктор. Получает значения при создании.
        public Jacobi(int n, double[,] Matrix, double[] FreeElements, double Accuracy)
        {
            this.n = n;
            this.matrix = Matrix;
            this.addtional = FreeElements;
            this.Accuracy = Accuracy;
            count_iteration = 0;

        }

        public int Iteration
        {
            get { return count_iteration; }
        }

        // Сам метод рассчета.
        public override void calculateMatrix()
        {
            // общий вид:
            // [x1]   [ b1/a11 ]   / 0 x x \ 
            // [x2] = [ b2/a22 ] - | x 0 x |
            // [x3]   [ b3/a33 ]   \ x x 0 /
            // где x - делится на диагональый элемент первоначальной матрицы.
            // где b - эелементы из свободных членов
            // где а - элементы из матрицы

            int iteration = 0, count = 0;


            // матрица коеффициентов + столбец свободных членов.
            double[,] a = new double[n, n + 1];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    a[i, j] = matrix[i, j];

            for (int i = 0; i < n; i++)
                a[i, n] = addtional[i];

            //---------------
            // Метод Якоби.
            //---------------

            // Введем вектор значений неизвестных на предыдущей итерации,
            // размер которого равен числу строк в матрице, т.е. size,
            // причем согласно методу изначально заполняем его нулями

            double[] previousValues = new double[n];
            for (int i = 0; i < n; i++)
            {
                previousValues[i] = 0.0;
            }
            /*for (int i = 0; i < n; i++)
            {
                previousValues[i] = addtional[i]/matrix[i,i];
            }*/

            // Будем выполнять итерационный процесс до тех пор,
            // пока не будет достигнута необходимая точность
            while (true)
            {
                // Введем вектор значений неизвестных на текущем шаге
                double[] currentValues = new double[n];

                // Посчитаем значения неизвестных на текущей итерации
                // в соответствии с теоретическими формулами
                for (int i = 0; i < n; i++)
                {
                    // Инициализируем i-ую неизвестную значением
                    // свободного члена i-ой строки матрицы
                    currentValues[i] = a[i, n];

                    // Вычитаем сумму по всем отличным от i-ой неизвестным
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                        {
                            currentValues[i] -= a[i, j] * previousValues[j];
                        }
                    }

                    // Делим на коэффициент при i-ой неизвестной
                    currentValues[i] /= a[i, i];
                }

                // Посчитаем текущую погрешность относительно предыдущей итерации
                double differency = 0.0;

                for (int i = 0; i < n; i++)
                    differency += Math.Abs(currentValues[i] - previousValues[i]);

                // Если необходимая точность достигнута, то завершаем процесс
                if (differency < accuracy)
                    break;
                    //count++;

                // Переходим к следующей итерации, так
                // что текущие значения неизвестных
                // становятся значениями на предыдущей итерации
                previousValues = currentValues;
                iteration++;
                //if (count == n)
                //    break;
                //count = 0;
            }
            count_iteration = iteration;
            resultMatrix = previousValues;
            //Console.WriteLine(iteration.ToString("|\t0.###") + "\t|");
        }

    }

    /// <summary>
    /// Класс Отвечает за работу метода Зейделя.
    /// Class Seidel
    /// </summary>
    class Seidel : SimpleIterations
    {
        // Матрица ответов
        private double[] resultMatrix;
        public double[] ResultMatrix
        {
            get
            {
                if (resultMatrix != null)
                    return resultMatrix;
                else
                {
                    return new double[n];
                }
            }
        }
        // private int t;

        // Основная матрица и свободные члены.
        private double[,] matrix;
        private double[] addtional;
        private int n;
        private int count_iteration;

        // точность (кол-во итераций)
        private double accuracy;
        // избегаем ошибок с итерациями.
        public double Accuracy
        {
            get
            {
                return accuracy;
            }
            set
            {
                if (value <= 0.0)
                    accuracy = 0.1;
                else
                    accuracy = value;
            }
        }

        // Конструктор. Получает значения при создании.
        public Seidel(int n, double[,] Matrix, double[] FreeElements, double Accuracy)
        {
            this.n = n;
            this.matrix = Matrix;
            this.addtional = FreeElements;
            this.Accuracy = Accuracy;
            count_iteration = 0;

        }

        public int Iteration
        {
            get { return count_iteration; }
        }

        // Сам метод рассчета.
        public override void calculateMatrix()
        {

            // общий вид:
            // [x1]   [ b1/a11 ]   / 0 x x \ 
            // [x2] = [ b2/a22 ] - | x 0 x |
            // [x3]   [ b3/a33 ]   \ x x 0 /
            // где x - делится на диагональый элемент первоначальной матрицы.
            // где b - эелементы из свободных членов
            // где а - элементы из матрицы

            int iteration = 0, count = 0;

            // матрица коеффициентов + столбец свободных членов.
            double[,] a = new double[n, n + 1];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    a[i, j] = matrix[i, j];

            for (int i = 0; i < n; i++)
                a[i, n] = addtional[i];

            //---------------
            // Метод Зейделя.
            //---------------

            // Введем вектор значений неизвестных на предыдущей итерации,
            // размер которого равен числу строк в матрице, т.е. size,
            // причем согласно методу изначально заполняем его нулями
            double[] previousValues = new double[n];
            for (int i = 0; i < n; i++)
            {
                previousValues[i] = 0.0;
            }

            // Будем выполнять итерационный процесс до тех пор,
            // пока не будет достигнута необходимая точность
            while (true)
            {
                // Введем вектор значений неизвестных на текущем шаге
                double[] currentValues = new double[n];

                // Посчитаем значения неизвестных на текущей итерации
                // в соответствии с теоретическими формулами
                for (int i = 0; i < n; i++)
                {
                    // Инициализируем i-ую неизвестную значением
                    // свободного члена i-ой строки матрицы
                    currentValues[i] = a[i, n];

                    // Вычитаем сумму по всем отличным от i-ой неизвестным
                    for (int j = 0; j < n; j++)
                    {
                        // При j < i можем использовать уже посчитанные
                        // на этой итерации значения неизвестных
                        if (j < i)
                        {
                            currentValues[i] -= a[i, j] * currentValues[j];
                        }
                        // При j > i используем значения с прошлой итерации
                        if (j > i)
                        {
                            currentValues[i] -= a[i, j] * previousValues[j];
                        }
                    }
                    // Делим на коэффициент при i-ой неизвестной
                    currentValues[i] /= a[i, i];
                }

                // Посчитаем текущую погрешность относительно предыдущей итерации
                double differency = 0.0;

                for (int i = 0; i < n; i++)
                    differency += Math.Abs(currentValues[i] - previousValues[i]);

                // Если необходимая точность достигнута, то завершаем процесс
                if (differency < accuracy)
                    break;
                    //count++;

                // Переходим к следующей итерации, так
                // что текущие значения неизвестных
                // становятся значениями на предыдущей итерации

                previousValues = currentValues;
                iteration++;
                //if (count == n)
                //    break;
                //count = 0;
            }
            count_iteration = iteration;
            // результат присваиваем матрице результатов.
            resultMatrix = previousValues;
            //Console.WriteLine(iteration.ToString("|\t0.###") + "\t|");
        }

    }
}