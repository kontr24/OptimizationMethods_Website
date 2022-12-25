using Accord.Math.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace OptimizationMethods_Website.Models
{
    public class MethodsRepository
    {
        public int k { get; set; }

        public List<string> _сonstantNewton = new List<string>();
        public List<string> _adjustmentNewton = new List<string>();
        public List<string> _np = new List<string>();


        //метод наискорейшего спуска
        public /*static */void steepestDescent(double[] x, double alpha, double tolerance)
        {
            int n = x.Length; //Размер входного массива
            double h = 1e-6;  //Коэффициент допуска
            double g0 = g(x); //Первоначальная оценка результата

            //Вычислить начальный градиент
            double[] fi = new double[n];
            fi = GradG(x, h);

            //Вычислить начальную норму
            double DelG = 0;
            for (int i = 0; i < n; ++i)
                DelG += fi[i] * fi[i];
            DelG = Math.Sqrt(DelG);

            double b = alpha / DelG;


            while (DelG > tolerance)
            {
                k++;
                for (int i = 0; i < n; ++i)
                    x[i] -= b * fi[i];
                h /= 2;

                //Вычислите следующий градиент
                fi = GradG(x, h);

                //
                DelG = 0;
                for (int i = 0; i < n; ++i)
                    DelG += fi[i] * fi[i];
                DelG = Math.Sqrt(DelG);

                b = alpha / DelG;

                //Проверка значения заданной функции с текущими значениями
                double g1 = g(x);

                //
                if (g1 > g0) alpha /= 2;
                else g0 = g1;


            }

        }

        //Обеспечивает приблизительный расчет градиента g(x).
        public /*static*/ double[] GradG(double[] x, double h)
        {
            int n = x.Length;
            double[] z = new double[n];
            double[] y = (double[])x.Clone();
            double g0 = g(x);
            for (int i = 0; i < n; ++i)
            {
                y[i] += h;
                z[i] = (g(y) - g0) / h;
            }
            return z;
        }

        //
        public /*static*/ double g(double[] x)
        {
            return 2 * (Math.Pow(x[1] - Math.Pow(x[0], 2), 2)) + (Math.Pow(1 - x[0], 2));
        }
        //метод наискорейшего спуска



        //Метод сопряженных градиентов
        static public double Fx(int N, ref double[] X, ref double[] fParam)
        {
            return 2 * (Math.Pow(X[1] - Math.Pow(X[0], 2), 2)) + (Math.Pow(1 - X[0], 2));
        }

        public delegate double MyFxDelegate(int nNumVars, ref double[] fX, ref double[] fParam);

        public class CConjugateGradient
        {
            MyFxDelegate m_MyFx;
            public double MyFxEx(int nNumVars, ref double[] fX, ref double[] fParam, ref double[] fDeltaX, double fLambda)
            {
                int i;
                double[] fXX = new double[nNumVars];

                for (i = 0; i < nNumVars; i++)
                {
                    fXX[i] = fX[i] + fLambda * fDeltaX[i];
                }

                return m_MyFx(nNumVars, ref fXX, ref fParam);
            }

            private void GetGradients(int nNumVars, ref double[] fX, ref double[] fParam, ref double[] fDeriv, ref double fDerivNorm)
            {
                int i;
                double fXX, H, Fp, Fm;
                fDerivNorm = 0;
                for (i = 0; i < nNumVars; i++)
                {
                    fXX = fX[i];
                    H = 0.01 * (1 + Math.Abs(fXX));
                    fX[i] = fXX + H;
                    Fp = m_MyFx(nNumVars, ref fX, ref fParam);
                    fX[i] = fXX - H;
                    Fm = m_MyFx(nNumVars, ref fX, ref fParam);
                    fX[i] = fXX;
                    fDeriv[i] = (Fp - Fm) / 2 / H;
                    fDerivNorm += Math.Pow(fDeriv[i], 2);
                }
                fDerivNorm = Math.Sqrt(fDerivNorm);
            }

            public bool LinSearch_DirectSearch(int nNumVars, ref double[] fX, ref double[] fParam, ref double fLambda, ref double[] fDeltaX, double InitStep, double MinStep)
            {
                double F1, F2;

                F1 = MyFxEx(nNumVars, ref fX, ref fParam, ref fDeltaX, fLambda);

                do
                {
                    F2 = MyFxEx(nNumVars, ref fX, ref fParam, ref fDeltaX, fLambda + InitStep);
                    if (F2 < F1)
                    {
                        F1 = F2;
                        fLambda += InitStep;
                    }
                    else
                    {
                        F2 = MyFxEx(nNumVars, ref fX, ref fParam, ref fDeltaX, fLambda - InitStep);
                        if (F2 < F1)
                        {
                            F1 = F2;
                            fLambda -= InitStep;
                        }
                        else
                        {
                            //уменьшить размер шага поиска
                            InitStep /= 10;
                        }
                    }
                } while (!(InitStep < MinStep));

                return true;

            }


            public double CalcOptim(int nNumVars, ref double[] fX, ref double[] fParam, double fEpsFx, int nMaxIter, ref int nIter, ref string sErrorMsg, MyFxDelegate MyFx)
            {

                int i;
                double[] fDeriv = new double[nNumVars];
                double[] fDerivOld = new double[nNumVars];
                double F, fDFNormOld, fLambda, fLastF, fDFNorm = 0;

                m_MyFx = MyFx;

                //вычисление функции значения в начальной точке
                fLastF = MyFx(nNumVars, ref fX, ref fParam);

                GetGradients(nNumVars, ref fX, ref fParam, ref fDeriv, ref fDFNorm);

                fLambda = 0.1;
                if (LinSearch_DirectSearch(nNumVars, ref fX, ref fParam, ref fLambda, ref fDeriv, 0.1, 0.000001))
                {
                    for (i = 0; i < nNumVars; i++)
                    {
                        fX[i] += fLambda * fDeriv[i];
                    }
                }
                else
                {
                    sErrorMsg = "Неудачный линейный поиск";
                    return fLastF;
                }

                nIter = 1;
                do
                {
                    nIter++;
                    if (nIter > nMaxIter)
                    {
                        sErrorMsg = "Достигнут максимальный предел итераций";
                        break;
                    }
                    fDFNormOld = fDFNorm;
                    for (i = 0; i < nNumVars; i++)
                    {
                        fDerivOld[i] = fDeriv[i]; // сохранить старый градиент
                    }
                    GetGradients(nNumVars, ref fX, ref fParam, ref fDeriv, ref fDFNorm);
                    for (i = 0; i < nNumVars; i++)
                    {
                        fDeriv[i] = Math.Pow((fDFNorm / fDFNormOld), 2) * fDerivOld[i] - fDeriv[i];
                    }
                    if (fDFNorm <= fEpsFx)
                    {
                        sErrorMsg = "Норма градиента соответствует критериям сходимости";
                        break;
                    }
                    fLambda = 0;
                    if (LinSearch_DirectSearch(nNumVars, ref fX, ref fParam, ref fLambda, ref fDeriv, 0.1, 0.000001))
                    {
                        for (i = 0; i < nNumVars; i++)
                        {
                            fX[i] += fLambda * fDeriv[i];
                        }
                        F = MyFx(nNumVars, ref fX, ref fParam);
                        if (Math.Abs(F - fLastF) < fEpsFx)
                        {
                            sErrorMsg = "Последовательные значения функции удовлетворяют критериям сходимости";
                            break;
                        }
                        else
                        {
                            fLastF = F;
                        }
                    }
                    else
                    {
                        sErrorMsg = "Неудачный линейный поиск";
                        break;
                    }
                } while (true);

                return fLastF;
            }

        }
        //Метод сопряженных градиентов


        //Метод с дроблением шага
        public double F(double x1, double x2)
        {
            return 2 * (Math.Pow(x2 - Math.Pow(x1, 2), 2)) + (Math.Pow(1 - x1, 2)); //заданная функция
        }
        public double prOne(double x1, double x2)//1-я производная по х1
        {
            return 8 * Math.Pow(x1, 3) - 8 * x1 * x2 + 2 * x1 - 2;
        }
        public double prTwo(double x1, double x2)//1-я производная по х2
        {
            return -4 * Math.Pow(x1, 2) + 4 * x2;
        }

        //Метод с дроблением шага




        //Метод Ньютона
        double f(double x1, double x2)
        {
            return 2 * (Math.Pow(x2 - Math.Pow(x1, 2), 2)) + (Math.Pow(1 - x1, 2)); //заданная функция

        }

        double f1(double x1, double x2)
        {

            return -8 * x1 * (-(Math.Pow(x1, 2)) + x2) + 2 * x1 - 2; //1-я производная по х1
        }

        double f2(double x1, double x2)
        {
            return -4 * Math.Pow(x1, 2) + 4 * x2; //1-я производная по х2
        }

        double f11(double x1, double x2)
        {
            return 24 * Math.Pow(x1, 2) - 8 * x2 + 2;
        }

        double f12(double x1, double x2)
        {
            return -8 * x1;
        }

        double f22(double x1, double x2)
        {
            return -8;
        }

        double fi(int ki, double dalf, double x1, double x2, double z)

        {
            return f(x1 - dalf * z, x2);
        }
        //Метод Ньютона


        //Метод Ньютона с постоянным шагом
        public int NewtonConstant(double[,] x, double accuracy)
        {

            double[,] y = new double[3, 3];

            int k = 0, N1, N2 = 0;

            double y1, y2, del, det, p1, p2;

            y1 = f1(x[0, 1], x[0, 2]);

            y2 = f2(x[0, 1], x[0, 2]);

            N1 = 2;

        p1:

            y[1, 1] = f11(x[k, 1], x[k, 2]);

            y[1, 2] = f12(x[k, 1], x[k, 2]);

            y[2, 1] = y[1, 2];

            y[2, 2] = f22(x[k, 1], x[k, 2]);

            N2 = N2 + 3;

            //x ̅[k]-∇f(x ̅[k])H^(-1) (x[k])
            det = y[1, 1] * y[2, 2] - y[1, 2] * y[2, 1];

            p1 = (y1 * y[2, 2] - y2 * y[1, 2]) / det;

            p2 = (y[1, 1] * y2 - y[2, 1] * y1) / det;

            x[k + 1, 1] = x[k, 1] - p1;

            x[k + 1, 2] = x[k, 2] - p2;
            //x ̅[k]-∇f(x ̅[k])H^(-1) (x[k])

            y1 = f1(x[k + 1, 1], x[k + 1, 2]);

            y2 = f2(x[k + 1, 1], x[k + 1, 2]);

            N1 = N1 + 2;

            del = Math.Pow(y1 * y1 + y2 * y2, 0.5); //‖∇f(x ̅[k]‖

            k++;

            if (del < accuracy /*k<100*/) goto p1;

            else
            {
                //function.Add("k=" + k + " N=" + N1 + N2 + " x=(" + x[k, 1] + ", " + x[k, 2] + ")\ny=" + f(x[k, 1], x[k, 2]));
                _сonstantNewton.Add("Количество итераций: " + k);
                _сonstantNewton.Add("x(" + " " + x[k, 1] + " " + "; " + " " + x[k, 2] + " " + ")");
                _сonstantNewton.Add("F(x) = " + f(x[k, 1], x[k, 2]));
                return 0;
            }
        }
        //Метод Ньютона с постоянным шагом



        //Метод Ньютона с регулировкой шага
        public int NewtonAdjustment(double[,] x, double accuracy)
        {

            double[,] y = new double[3, 3];

            int k = 0, N1, N2 = 0, i, l, n = 10;

            double y1, y2, del, det, p1, p2;

            double m, a1, b1, ym1, ym2, min, ma, kan;

            double dalf = 1;
            y1 = f1(x[0, 1], x[0, 2]);

            y2 = f2(x[0, 1], x[0, 2]);

            N1 = 2;


        p1:

            y[1, 1] = f11(x[k, 1], x[k, 2]);

            y[1, 2] = f12(x[k, 1], x[k, 2]);

            y[2, 1] = y[1, 2];

            y[2, 2] = f22(x[k, 1], x[k, 2]);

            N2 = N2 + 3;

            det = y[1, 1] * y[2, 2] - y[1, 2] * y[2, 1];

            p1 = (y1 * y[2, 2] - y2 * y[1, 2]) / det;

            p2 = (y[1, 1] * y2 - y[2, 1] * y1) / det;

            m = 0;

            ym1 = f(x[k, 1], x[k, 2]);

            N1++;

        ab:

            ym2 = fi(k, dalf * (m + 1), x[k, 1], x[k, 2], y1);

            N1++;

            if (ym2 < ym1)

            { m++; ym1 = ym2; goto ab; }

            else

            {
                b1 = (m + 1) * dalf;

                if (m == 0) a1 = 0; else a1 = (m - 1) * dalf;
            }

            ma = a1 + (b1 - a1) / n;

            min = f(x[k, 1] - ma * y1, x[k, 2] - ma * y2); N1++;

            l = n;

            for (i = 2; i <= n; i++)

            {
                ma = a1 + i * (b1 - a1) / n;

                kan = f(x[k, 1] - ma * y1, x[k, 2] - ma * y2); N1++;

                if (min > kan) { min = kan; l = i; }

            }
            //x ̅[k+1]=x ̅[k]-α_k ∇f(x ̅[k])H^(-1) (x[k])
            x[k + 1, 1] = x[k, 1] - (l * (b1 - a1) / n) * p1;

            x[k + 1, 2] = x[k, 2] - (l * (b1 - a1) / n) * p2;

            y1 = f1(x[k, 1], x[k, 2]);

            y2 = f2(x[k, 1], x[k, 2]);

            del = Math.Pow(y1 * y1 + y2 * y2, 0.5); //‖∇f(x ̅[k]‖

            k++;

            if (del < accuracy  /*k < 1*/) goto p1;

            else
            {
                //label1.Text += "k=" + k + " N=" + N1 + N2 + " x=(" + x[k, 1] + ", " + x[k, 2] + ")\ny=" + f(x[k, 1], x[k, 2]);

                _adjustmentNewton.Add("Количество итераций: " + k);
                _adjustmentNewton.Add("Количество шагов: " + l);
                _adjustmentNewton.Add("x(" + " " + x[k, 1] + " " + "; " + " " + x[k, 2] + " " + ")");
                _adjustmentNewton.Add("F(x) = " + f(x[k, 1], x[k, 2]));
                return 0;
            }
        }
        //Метод Ньютона с регулировкой шага



        //ЗНП//3 практическая работа (задание 3.23)
        public void Solve(NonlinearObjectiveFunction f, NonlinearConstraint[] constraints, double accuracy, int number)
        {
            var cobyla = new Cobyla(f, constraints);

            if (number == 0)
            {
                bool success = cobyla.Maximize(); //находит максимальное значение функции
            }
            else
            {
                bool success = cobyla.Minimize();
            }
            //_np.Add(success.ToString());

            double minimum = cobyla.Value; //минимальное значение функции
            double[] solution = cobyla.Solution; //Получает текущее найденное решение, значения параметров, которые оптимизируют функцию
            int counter = 0;

            foreach (double item in solution)
            {
                _np.Add("x" + (++counter) + " = " + item);
            }

            _np.Add("F(Xmin) = " + minimum.ToString());// минимум функции 

        }

        //public NonlinearConstraint[] MakeByHands()
        //{
        //    //Нелинейные ограничения
        //    var constraints = new[]
        //     {
        //        new NonlinearConstraint(3, x => 2* x[0] + x[1] -x[2] <= -6),
        //        new NonlinearConstraint(3, x => 4* x[0] + x[1] - x[2] <= 20),
        //        new NonlinearConstraint(3, x => x[1] >=0),
        //     };
        //    return constraints;
        //    //Нелинейные ограничения
        //}

        
        public NonlinearConstraint[] MakeAutomatic(int signFirst, double firstRestriction, int signTwo, double twoRestriction, int signThree, double threeRestriction)
        {
            List<NonlinearConstraint> ConstraintList = new List<NonlinearConstraint>();
            //Нелинейные ограничения
            if (signFirst == 0)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 2 * x[0] + x[1] - x[2] <= firstRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 2 * x[0] + x[1] - x[2] >= firstRestriction));
            }
            if (signTwo == 0)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 4 * x[0] + x[1] - x[2] <= twoRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 4 * x[0] + x[1] - x[2] >= twoRestriction));
            }
            if (signThree == 1)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[1] >= threeRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[1] <= threeRestriction));
            }
            //Нелинейные ограничения

            foreach (var j in Enumerable.Range(0, 3))
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[j] >= 0));
            }
            return ConstraintList.ToArray();
        }
        //ЗНП//3 практическая работа (задание 3.23)



        //ЗНП//3 практическая работа (задание 3.24)
        public NonlinearConstraint[] MakeAutomaticTwo(int signFirst, double firstRestriction, int signTwo, double twoRestriction, int signThree, double threeRestriction)
        {
            List<NonlinearConstraint> ConstraintList = new List<NonlinearConstraint>();

            //Нелинейные ограничения
            if (signFirst == 0)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 2 * x[0] + x[1] - x[2] + 6 <= firstRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 2 * x[0] + x[1] - x[2] + 6 >= firstRestriction));
            }
            if (signTwo == 0)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 4 * x[0] - x[1] + 3 * x[2] <= twoRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => 4 * x[0] - x[1] + 3 * x[2] >= twoRestriction));
            }
            if (signThree == 1)
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[1] >= threeRestriction));
            }
            else
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[1] <= threeRestriction));
            }
            ////Нелинейные ограничения
            foreach (var j in Enumerable.Range(0, 3))
            {
                ConstraintList.Add(new NonlinearConstraint(3, x => x[j] >= 0));
            }
            return ConstraintList.ToArray();
        }
        //ЗНП//3 практическая работа (задание 3.24)










        //Метод Ньютена (ненужный код)
        public int Nyuton(double[] x, double step, double tolerance, int countIteration)
        {
            const int n = 2;            //размерность вектора    
            //x = new double[n];                //вектор координат исходного минимума
            int ntrial = countIteration;        //максимально число шагов для метода Ньютона
            double eps = tolerance;    //точность 
            int iter;                    //номер итерационного шага в методе Ньютона
            int i, j, k;                //переменные циклов
            double tmp;                    //временная переменная
            double[] fvec = new double[n];                //вектор для хранения частных производных 
            double[] p = new double[n];                //вектор для хранения вектора смещений 
            double[,] fjac = new double[n, n];            //матрица Гессе 

            //Шаг 1
            //Можно задать любое начальное приближение
            //x[0] = 3;
            //x[1] = 4;
            //Итерационный процесс метода Ньютона

            for (iter = 1; iter <= ntrial; iter++)
            {
                _сonstantNewton.Add("Итерация " + (iter - 1) + ": " + "x(" + " " + x[0] + " " + "; " + " " + x[1] + " " + ")" + "; F(x) = " + F(x[0], x[1]));

                //Шаг 2  
                fvec[0] = 8 * Math.Pow(x[0], 3) - 8 * x[0] * x[1] + 2 * x[0] - 2;    //df/dx1 = 8 * Math.Pow(x[0], 3) - 8 * x[0] * x[1] + 2 * x[0] - 2
                fvec[1] = -4 * Math.Pow(x[0], 2) + 4 * x[1];    //df/dx2 = -4 * Math.Pow(x[0], 2) + 4 * x[1]

                tmp = 0;
                for (i = 0; i < n; i++)
                {
                    tmp += Math.Abs(fvec[i]);
                }
                if (tmp <= eps)
                    break;

                //Шаг 3
                fjac[0, 0] = 1;        //dF1(x1, x2)/dx1 = 1        
                fjac[0, 1] = 1;         //dF1(x1, x2)/dx1 = 1
                fjac[1, 0] = 1;       //dF1(x1, x2)/dx1 = 1    
                fjac[1, 1] = 8;       //dF1(x1, x2)/dx2 = 8

                //NewtonInfo += "Матрица Гессе в векторе x:";
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        //NewtonInfo += fjac[i, j] + "\t";
                    }
                    //NewtonInfo += "\n";
                }
                //Шаг 4
                //Правая часть линейной системы 
                for (i = 0; i < n; i++)
                    p[i] = -fvec[i];
                //Прямой ход метода
                for (i = 0; i < n - 1; i++)
                {
                    tmp = fjac[i, i];

                    if (tmp == 0)
                    {
                        _сonstantNewton.Add("Ошибка: нулевой элемент на главной диагонали в матрице Гессе");
                        return 0;
                    }

                    // рабочий элемент корректен - обрабатываем все нижлежащие строки
                    for (j = i + 1; j < n; j++)
                    {
                        p[j] -= p[i] / tmp;
                        for (k = n - 1; k >= i; k--)
                        {
                            fjac[j, k] -= fjac[i, k] * fjac[j, i] / tmp;
                        }
                    }
                }
                //Нахождение вектора смещений (обратный ход метода)
                for (i = n - 1; i >= 0; i--)
                {
                    tmp = 0;
                    for (j = i + 1; j < n; j++)
                    {
                        tmp += fjac[i, j] * p[j];
                    }
                    p[i] = (p[i] - tmp) / fjac[i, i];
                }
                //Шаг 5
                for (i = 0; i < n; i++)
                {
                    x[i] += p[i];
                }

                //Шаг 6
            }
            _сonstantNewton[0] = "";
            //if (iter < ntrial)
            //{
            //    //NewtonInfo += "\n\nМинимальное значение\n";
            //    //NewtonInfo += "(" + x[0] + ", " + x[1] + ")\n";
            //    return 0;
            //}
            _сonstantNewton.Add("Максимальное количество итераций превышено");
            return 0;
        }
        //Метод Ньютена (ненужный код)

    }
}
