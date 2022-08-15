using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
namespace NumSharp
{
    /// <summary>
    /// 
    /// </summary>
   public class SpecialFunctions
    {
        private static readonly Vector factorialLookup =new double[]
       { 
            1.0,
            1.0,
            2.0,
            6.0,
            24.0,
            120.0,
            720.0,
            5040.0,
            40320.0,
            362880.0,
            3628800.0,
            39916800.0,
            479001600.0,
            6227020800.0,
            87178291200.0,
            1307674368000.0,
            20922789888000.0,
            355687428096000.0,
            6402373705728000.0,
            121645100408832000.0,
            2432902008176640000.0,
            51090942171709440000.0,
            1124000727777607680000.0,
            25852016738884976640000.0
        };
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double Factorial(int n)
        {
            if (n < 0)
            {
                throw new Exception("Input value must be > 0");
            }
            else if (n < factorialLookup.Length)
            {
                return factorialLookup[n];
            }
            else
            {
                return Floor(Exp(FactorialLn(n)) + 0.5);
            }
        }
        private static double FactorialLn(int n)
        {
            if (n < 0)
            {
                throw new Exception("Input value must be > 0");
            }
            else
            {
                return GammaLn(n + 1.0);
            }
        }
        private static double GammaLn(double x)
        {
            if (x <= 0) throw
            new Exception($"Input value must be > 0");
            Vector coef = new double[14]
            {
                57.1562356658629235,
                -59.5979603554754912,
                14.1360979747417471,
                -0.491913816097620199,
                0.339946499848118887E-4,
                0.465236289270485756E-4,
                -0.983744753048795646E-4,
                0.158088703224912494E-3,
                -0.210264441724104883E-3,
                0.217439618115212643E-3,
                -0.164318106536763890E-3,
                0.844182239838527433E-4,
                -0.261908384015814087E-4,
                0.368991826595316234E-5
            };
            double denominator = x;
            double series = 0.999999999999997092;
            double temp = x + 5.24218750000000000;
            temp = (x + 0.5) *Log(temp) - temp;
            for (int j = 0; j < 14; j++)
                series += coef[j] / ++denominator;
            return (temp + Log(2.5066282746310005 * series / x));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double Gamma(double x)
        {
            if (Sign(x) == -1)//Apply the repflection rule
            {
                if (x == (int)x) throw new Exception("Input value must not be a negative interger");
                double xx = 1 + Abs(x);
                double gg = Exp(GammaLn(xx));
                double res = PI / (gg * Sin(PI * x));
                return res;

            }
            return Exp(GammaLn(x));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double Combination(int n,int k)
        {
            if ((n < 0 || k < 0) || (n < k))throw new Exception("Input value must be > 0");
            if (n > 170) return Floor( (0.5 + Exp(FactorialLn(n) - FactorialLn(k) - FactorialLn(n - k))));
            return (Factorial(n) / (Factorial(k) * Factorial(n - k)));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static double Permutation(double n, double k)
        {
            return (Gamma(n) * Gamma(k)) / Gamma(n + k);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double Beta(int n, int k)
        {
            if ((n < 0 || k < 0) || (n < k))throw new Exception("Input value must be > 0");
            return Math.Floor(0.5 + (Math.Exp(FactorialLn(n) -FactorialLn(k) - FactorialLn(n - k))));
        }
    }
}
