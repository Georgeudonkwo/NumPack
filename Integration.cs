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
    public enum OrthogonalPolynomial
    {
        /// <summary>
        /// 
        /// </summary>
        Legendre,
        /// <summary>
        /// 
        /// </summary>
        Laguerre,
        /// <summary>
        /// 
        /// </summary>
        Hermite,
        /// <summary>
        /// 
        /// </summary>
        Chebyshev
    }
    /// <summary>
    /// 
    /// </summary>
   public class Integration
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="F"></param>
        /// <param name="Limits"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double RombergIntegration(Func<double, double> F,( double Lower, double upper)Limits, int n)
        {
            double[][] R = new double[n][];
            for (int i = 0; i < R.GetLength(0); i++)
            {
                double h = (Limits.upper - Limits.Lower) / Pow(2, i);
                if (i == 0)
                {
                    double[] c = new double[1];
                    c[0] = h / 2 * (F(Limits.Lower) - F(Limits.upper));
                    R[i] = c; continue;
                }
                else
                {
                    double sum = 0.0;
                    for (int j = 0; j < Pow(2, i - 1); j++)
                    {
                        sum += F(Limits.Lower + (2 * j - 1) * h);
                    }
                    double[] c1 = new double[1];
                    c1[0] = 1.0 / 2.0 * R[i - 1][0] + h * sum;
                    R[i] = c1;
                }
                double[] inneraaray = new double[i + 1];
                inneraaray[0] = R[i][0];
                for (int j = 1; j <= i; j++)
                {
                    R[i] = inneraaray;
                    double value = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (Pow(4, j) - 1.0);
                    inneraaray[j] = value;
                }
                R[i] = inneraaray;
            }
            double result = R.Last().Last();
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="F"></param>
        /// <param name="Limits"></param>
        /// <param name="level"></param>
        /// <param name="maxlevel"></param>
        /// <returns></returns>
        public static double AdaptiveSimpson(Func<double, double> F, (double LB, double UB)Limits, int level, int maxlevel)
        {
            //double a, b, c, d, e, h, simpRes;
            double c, d, e, h, simpRes;
            level++;
            h = Limits.UB - Limits.LB;
            c = (Limits.LB + Limits.UB) / 2.0;
            double simp1 = h * (F(Limits.LB) + 4 * F(c) + F(Limits.UB)) / 6.0;
            d = (Limits.LB + c) / 2.0;
            e = (c + Limits.UB) / 2.0;
            double simp2 = h * (F(Limits.LB) + 4 * F(d) + 2 * F(c) + 4 * F(e) + F(Limits.UB)) / 12.0;
            if (level >= maxlevel)
            {
                simpRes = simp2;
                return simpRes;
            }
            else
            {
                if (Abs(simp2 - simp1) < 15 * 0.001)
                {
                    simpRes = simp2 + (simp2 - simp1) / 15.0;
                }
                else
                {
                    double Lsimp = AdaptiveSimpson(F, (Limits.LB, c), level, maxlevel);
                    double Rsimp = AdaptiveSimpson(F, (c, Limits.UB), level, maxlevel);
                    simpRes = Rsimp + Lsimp;
                }
            }
            return simpRes;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="Limits"></param>
        /// <param name="QuadratureType"></param>
        /// <returns></returns>
        public static  double GaussQuadrature(Func<double,double>Fun,(double Lower,double Upper) Limits, OrthogonalPolynomial QuadratureType)
        {
            (Vector R, Vector W) WN = default(ValueTuple<Vector,Vector>);
            double Integrand = 0.0;
            double IntegrandC = 0.0;
            int n = 2;
            bool exit = true;
            switch (QuadratureType)
            {
                case OrthogonalPolynomial.Legendre:
                    WN = Polynomial.LegendreNodesWeights(n);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x =  0.5* (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        double funcVal = Fun(x);
                        Integrand += WN.W[i] * funcVal;
                    }
                    Integrand = 0.5 * (Limits.Upper - Limits.Lower) * Integrand;
                    do
                    {
                        n++;
                        IntegrandC = 0.0;
                        WN = Polynomial.LegendreNodesWeights(n);
                        for (int i = 0; i < WN.R.Length; i++)
                        {
                            double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                            double funcVal = Fun(x);
                            IntegrandC += WN.W[i] * funcVal;
                        }
                         IntegrandC = 0.5 * (Limits.Upper - Limits.Lower) * IntegrandC;
                        exit = Abs(Integrand - IntegrandC) < 1e-6;
                        Integrand = IntegrandC;
                    } while (!exit);
                    break;
                case OrthogonalPolynomial.Laguerre:
                    WN = Polynomial.LaguerreNodesWeights(n);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = WN.R[i];
                        double funcVal = Fun(x);
                        Integrand += WN.W[i] * funcVal;
                    }
                    do
                    {
                        n++;
                        IntegrandC = 0.0;
                        WN = Polynomial.LaguerreNodesWeights(n);
                        for (int i = 0; i < WN.R.Length; i++)
                        {
                            double x = WN.R[i];
                            double funcVal = Fun(x);
                            IntegrandC += WN.W[i] * funcVal;
                        }
                        exit = Abs(Integrand - IntegrandC) < 1e-6;
                        Integrand = IntegrandC;

                    } while (!exit);
                    break;
                case OrthogonalPolynomial.Hermite:
                    WN = Polynomial.HermiteNodesWeights(n);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = WN.R[i];
                        double funcVal = Fun(x);
                        Integrand += WN.W[i] * funcVal;
                    }
                    do
                    {
                        n++;
                        IntegrandC = 0.0;
                        WN = Polynomial.HermiteNodesWeights(n);
                        for (int i = 0; i < WN.R.Length; i++)
                        {
                            double x = WN.R[i];
                            double funcVal = Fun(x);
                            IntegrandC += WN.W[i] * funcVal;
                        }
                        exit = Abs(Integrand - IntegrandC) < 1e-6;
                        Integrand = IntegrandC;

                    } while (!exit);
                    break;
                case OrthogonalPolynomial.Chebyshev:
                    WN = Polynomial.ChebyshevRoots(n);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        //double x = WN.R[i];
                        double funcVal = Fun(x);
                        Integrand += funcVal;
                    }
                    Integrand = WN.W[0] * Integrand;
                    do
                    {
                        n++;
                        IntegrandC = 0.0;
                        WN = Polynomial.ChebyshevRoots(n);
                        for (int i = 0; i < WN.R.Length; i++)
                        {
                            double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                            //double x = WN.R[i];
                            double funcVal = Fun(x);
                            IntegrandC += funcVal;
                        }
                        IntegrandC = WN.W[0] * IntegrandC;
                        exit = Abs(Integrand - IntegrandC) < 1e-6;
                        Integrand = IntegrandC;
                    } while (!exit);
                    break;
                default:
                    break;
            }
            return Integrand;
        }
         static double GaussQuadrature(Func<double,double, double> Fun,double y, (double Lower, double Upper) Limits, OrthogonalPolynomial QuadratureType)
        {
            (Vector R, Vector W) WN = default(ValueTuple<Vector, Vector>);
            double Integrand = 0.0;
            switch (QuadratureType)
            {
                case OrthogonalPolynomial.Legendre:
                    WN = Polynomial.LegendreNodesWeights(5);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        double funcVal = Fun(x,y);
                        Integrand += WN.W[i] * funcVal;
                    }
                    Integrand = 0.5 * (Limits.Upper - Limits.Lower) * Integrand;
                    break;
                case OrthogonalPolynomial.Laguerre:
                    WN = Polynomial.LaguerreNodesWeights(3);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        //double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        double x = WN.R[i];
                        double funcVal = Fun(x,y);
                        Integrand += WN.W[i] * funcVal;
                    }
                    //Integrand = 0.5 * (Limits.Upper - Limits.Lower) * Integrand;
                    break;
                case OrthogonalPolynomial.Hermite:
                    WN = Polynomial.HermiteNodesWeights(5);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        //double x = WN.R[i];
                        double funcVal = Fun(x,y);
                        Integrand += WN.W[i] * funcVal;
                    }
                    Integrand = 0.5 * (Limits.Upper - Limits.Lower) * Integrand;
                    break;
                case OrthogonalPolynomial.Chebyshev:
                    WN = Polynomial.ChebyshevRoots(8);
                    for (int i = 0; i < WN.R.Length; i++)
                    {
                        double x = 0.5 * (Limits.Upper + Limits.Lower) + 0.5 * WN.R[i] * (Limits.Upper - Limits.Lower);
                        //double x = WN.R[i];
                        double funcVal = Fun(x,y);
                        Integrand += funcVal;
                    }
                    Integrand = WN.W[0] * Integrand;
                    Integrand = 0.5 * (Limits.Upper - Limits.Lower) * Integrand;
                    break;
                default:
                    break;
            }
            return Integrand;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="L1"></param>
        /// <param name="L2"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="QuadratureType"></param>
        /// <returns></returns>
        public static double MultipleIntegral(Func<double, double, double> Fun, (double Lower, double Upper) L1, (double Lower, double Upper) L2, double dx, double dy, OrthogonalPolynomial QuadratureType)
        {
            Func<double,(double,double), double> F = new Func<double, (double,double), double>((x,h) =>
               {
                   return GaussQuadrature(Fun,dx, h, QuadratureType);
               });
            return F(dy,L2);

        }
    }
}
