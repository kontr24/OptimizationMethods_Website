using Accord.Math.Optimization;
using Microsoft.AspNetCore.Mvc;
using Microsoft.Extensions.Logging;
using OptimizationMethods_Website.Models;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using static OptimizationMethods_Website.Models.MethodsRepository;

namespace OptimizationMethods_Website.Controllers
{
    public class HomeController : Controller
    {
        private readonly ILogger<HomeController> _logger;

        public HomeController(ILogger<HomeController> logger)
        {
            _logger = logger;
        }

        public IActionResult Index()
        {
            return View();
        }

        [HttpPost]
        public IActionResult Index(double pointOne, double pointTwo, double accuracy)
        {
            ViewBag.pointOne = pointOne;
            ViewBag.pointTwo = pointTwo;
            ViewBag.accuracy = accuracy;
            //Метод наискорейшего спуска
            double tolerance = ViewBag.accuracy; //точность
            double alpha = 0.1;// шаг
            double[] x = new double[2];
            //начальная точка
            x[0] = pointOne;
            x[1] = pointTwo;
            //начальная точка
            MethodsRepository methodsRepository = new MethodsRepository();
            methodsRepository.steepestDescent(x, alpha, tolerance);


            //ViewBag.result = "Результат";
            ViewBag.minimum = "x (" + " " + Math.Round(x[0], 4) + " " + ";" + " " + Math.Round(x[1], 4) + " " + ")";
            ViewBag.function = "F(x)min: " + Math.Round(methodsRepository.g(x), 4);

            ViewBag.iterations = "Количество итераций: " + methodsRepository.k;
            methodsRepository.k = 0;
            //Метод наискорейшего спуска


            //Метод сопряженных градиентов
            int nNumVars = 2;
            double[] fX = new double[] { pointOne, pointTwo };//начальная точка
            double[] fParam = new double[] { 0, 0 };
            int nIter = 0;
            int nMaxIter = 100;
            double fEpsFx = accuracy; //точность
            int i;
            double fBestF;
            string sErrorMsg = "";
            CConjugateGradient oOpt;
            MyFxDelegate MyFx = new MyFxDelegate(Fx);
            oOpt = new CConjugateGradient();

            fBestF = oOpt.CalcOptim(nNumVars, ref fX, ref fParam, fEpsFx, nMaxIter, ref nIter, ref sErrorMsg, MyFx);

            ViewBag.errorMsg = sErrorMsg;
            //ViewBag.resultG = "Результат";
            for (i = 0; i < nNumVars; i++)
            {
                ViewBag.minimumG = "x (" + " " + Math.Round(fX[0], 4) + " " + ";" + " " + Math.Round(fX[1], 4) + " " + ")";
            }
            ViewBag.functionG = "F(x)min: " + Math.Round(fBestF, 4);
            ViewBag.iterationsG = "Количество итераций: " + nIter;
            //Метод сопряженных градиентов

            //Метод с дроблением шага
            double x1_prev, x2_prev;
            double x1 = pointOne, x2 = pointTwo, eps = accuracy, alphaC = 0.5;//alphaС - шаг; x1,x2 - начальные точки; eps - точность
            double dx1;
            double dx2;
            //double Y1, Y2;
            int max_iter = 100;

            int k = 0, l = 0;
            List<string> function = new List<string>();
            do
            {
                x1_prev = x1;
                x2_prev = x2;

                dx1 = methodsRepository.prOne(x1, x2);//вычисляем 1-ю производную по х1
                dx2 = methodsRepository.prTwo(x1, x2);//вычисляем 1-ю производную по х2    

                x1 = x1 - alphaC * dx1;//новое х1 
                x2 = x2 - alphaC * dx2;//новое х2

                //Y1 = methodsRepository.F(x1_prev, x2_prev);
                //Y2 = methodsRepository.F(x1, x2);

                k++;
                alphaC = alphaC / 10;

                function.Add("Итерация " + k + ": x (" + " " + x1 + " " + ";" + " " + x2 + " " + ")" + "; F(x) = " + methodsRepository.F(x1, x2));

                if (k > max_iter) break;

            } while ((Math.Abs(x1 - x1_prev) >= eps) && (Math.Abs(x2 - x2_prev) > eps));

            function.Add("Решение: Количество итераций = " + k + "; x (" + " " + x1 + " " + ";" + " " + x2 + " " + ")" + "; F(x) = " + methodsRepository.F(x1, x2));
            ViewBag.functionC = function;
            //Метод с дроблением шага
            return View();
        }


        public IActionResult Privacy()
        {
            return View();
        }
        [HttpPost]
        public IActionResult Privacy(double pointOne, double pointTwo, double accuracy/*, int countIteration*/)
        {
            ViewBag.pointOne = pointOne;
            ViewBag.pointTwo = pointTwo;
            ViewBag.accuracy = accuracy;
            //ViewBag.countIteration = countIteration;

            MethodsRepository methods = new MethodsRepository();
            double tolerance = accuracy; //точность
            double[,] x = new double[10000, 3];
            //начальная точка
            x[0, 1] = pointOne;
            x[0, 2] = pointTwo;
            //начальная точка

            methods.NewtonConstant(x, tolerance);//Метод Ньютона с постоянным шагом
            ViewBag.сonstantNewton = methods._сonstantNewton;


            methods.NewtonAdjustment(x, tolerance);//Метод Ньютона с регулировкой шага
            ViewBag.adjustmentNewton = methods._adjustmentNewton;
            //ViewBag.step = 

            return View();
        }

        //3 практическая работа (задание 3.23)
        public IActionResult NonlinearProgramming()
        {
            ViewBag.accuracy = 0.001;
            ViewBag.signFirst = 0;
            ViewBag.firstRestriction = -6;

            ViewBag.signTwo = 0;
            ViewBag.twoRestriction = 20;

            ViewBag.signThree = 1;
            ViewBag.threeRestriction = 0;

            return View();
        }
        [HttpPost]
        public IActionResult NonlinearProgramming(double accuracy, int signFirst, double firstRestriction, int signTwo, double twoRestriction, int signThree, double threeRestriction, int number = 0)
        {
            MethodsRepository methods = new MethodsRepository();

            ViewBag.accuracy = accuracy;
            ViewBag.signLess = signFirst;
            ViewBag.firstRestriction = firstRestriction;

            ViewBag.signTwo = signTwo;
            ViewBag.twoRestriction = twoRestriction;

            ViewBag.signThree = signThree;
            ViewBag.threeRestriction = threeRestriction;


            var f = new NonlinearObjectiveFunction(3, x => 2 * Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 4 * x[1] - 3 * x[2]); //функция

            //var constraintsOne = methods.MakeByHands();
            //methods.Solve(f, constraintsOne);

            var constraintsTwo = methods.MakeAutomatic(signFirst, firstRestriction, signTwo, twoRestriction, signThree, threeRestriction);
            methods.Solve(f, constraintsTwo, accuracy, number);

            ViewBag.np = methods._np;


            return View();
        }
        //3 практическая работа (задание 3.23)

        //3 практическая работа (задание 3.24)
        public IActionResult NonlinearProgrammingTwo()
        {
            ViewBag.accuracy = 0.001;
            ViewBag.signFirst = 0;
            ViewBag.firstRestriction = 0;

            ViewBag.signTwo = 0;
            ViewBag.twoRestriction = 20;

            ViewBag.signThree = 1;
            ViewBag.threeRestriction = 0;

            return View();
        }
        [HttpPost]
        public IActionResult NonlinearProgrammingTwo(double accuracy, int signFirst, double firstRestriction, int signTwo, double twoRestriction, int signThree, double threeRestriction, int number = 1)
        {
            MethodsRepository methods = new MethodsRepository();

            ViewBag.accuracy = accuracy;
            ViewBag.signLess = signFirst;
            ViewBag.firstRestriction = firstRestriction;

            ViewBag.signTwo = signTwo;
            ViewBag.twoRestriction = twoRestriction;

            ViewBag.signThree = signThree;
            ViewBag.threeRestriction = threeRestriction;


            var f = new NonlinearObjectiveFunction(3, x => 2 * Math.Pow(x[0], 2) + Math.Pow(x[1], 2) + Math.Pow(x[2], 2) - 4 * x[1] - 3 * x[2]); //функция
            //var constraintsOne = methods.MakeByHands();
            //methods.Solve(f, constraintsOne);

            var constraintsTwo = methods.MakeAutomaticTwo(signFirst, firstRestriction, signTwo, twoRestriction, signThree, threeRestriction);
            methods.Solve(f, constraintsTwo, accuracy, number);

            ViewBag.np = methods._np;


            return View();
        }
        //3 практическая работа (задание 3.24)

        //[ResponseCache(Duration = 0, Location = ResponseCacheLocation.None, NoStore = true)]
        //public IActionResult Error()
        //{
        //    return View(new ErrorViewModel { RequestId = Activity.Current?.Id ?? HttpContext.TraceIdentifier });
        //}
    }
}
