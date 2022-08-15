using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using static System.Math;
namespace NumSharp
{
    /// <summary>
    /// 
    /// </summary>
  public class Polynomial
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Coeffs"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        public static double polyEval(Vector Coeffs, double val)
        {
            Vector CoeffsRev = Coeffs;//.Reverse();
            double lastterm = CoeffsRev[CoeffsRev.Length - 1];
            double secondToLastterm = CoeffsRev[CoeffsRev.Length - 2];
            double innerProd = secondToLastterm + (val * lastterm);
            if (Coeffs.Length == 2) return innerProd;
            double sum = val * innerProd;

            for (int i = CoeffsRev.Length - 3; i >= 1; i--)
            {
                sum = val * (CoeffsRev[i] + sum);
            }
            sum = CoeffsRev[0] + sum;
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Coeffs"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        public static (double value,double FirstDerivative,double SecondDerivative) polyEval2(Vector Coeffs, double val)
        {
            Vector CoeffsRev = Coeffs.Copy();
            double p = CoeffsRev[CoeffsRev.Length - 1];
            double dp= 0.0;
            double ddp = 0.0;
            for (int i = CoeffsRev.Length - 2; i >= 0; i--)
            {
                ddp = dp * val + 2.0 * dp;
                dp = dp * val + p;
                p = CoeffsRev[i] + p*val;
            }
           double der = dp;
           double der2 = ddp;
            return (p,der,der2);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="coeffs"></param>
        /// <param name="R"></param>
        /// <returns></returns>
        public static Vector PolyDeflate(Vector coeffs,double R)
        {
            Vector b = new Vector(coeffs.Length - 1);
            b[b.Length - 1] = coeffs[coeffs.Length - 1];
            for (int i = b.Length-2; i >=0; i--)
            {
                b[i] = coeffs[i + 1] + R * b[i + 1];
            }
            return b;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Coeffs"></param>
        /// <returns></returns>
        public static Vector PolyRoots(Vector Coeffs)
        {
            int n = Coeffs.Length;
            int nn = n-1;
            double tol = 1e-12;
          
            Vector roots = new Vector(n - 1);
           
            for (int i = 0; i < nn; i++)
            {
                int count = 0;
                Vector c = Coeffs.Absolute() / Abs(Coeffs[Coeffs.Length - 1]);
                double Cma = c.Max()+1;
                double x = Cma;
                double diff = 0.0;
                bool converged = true;
                do
                {
                    count++;
                    var Peval = polyEval2(Coeffs, x);
                    if (Abs(Peval.value) < tol)
                    {
                        roots[i] = x;
                        break;
                    }
                    double g = Peval.FirstDerivative / Peval.value;
                    double h = Pow(g, 2) - Peval.SecondDerivative / Peval.value;
                    double s = (n - 1) * (n * h - Pow(g, 2));
                    double f = Sqrt(Abs(s));
                    ////var fff = Complex.Sqrt(new Complex(s,0));
                    ////double f = fff.Real;
                    double dx = 0.0;
                    if (Abs(g + f) > Abs(g - f))
                    {
                        dx = n / (g + f);
                    }
                    else
                    {
                        dx = n / (g - f);
                    }
                    x = x - dx;
                    diff = Abs(dx);
                    converged = Abs(dx) > tol;
                    if (count > 200) break;

                } while (converged);
                roots[i] = x;
                Coeffs = PolyDeflate(Coeffs, x);
                n = Coeffs.Length;
            }
            return roots;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Coeffs"></param>
        /// <param name="mat"></param>
        /// <returns></returns>
        public static Vector PolyRoots2(Vector Coeffs,out Matrix mat)
        {
            Vector c = (-1)*( Coeffs / Coeffs[0]);
            c.Remove(0);
            Matrix M = Matrix.Zeros(c.Length, c.Length);
            M[0, "row"] = c;
            Matrix eye = Matrix.IdentityMatrix(c.Length - 1);
            Matrix.InsectSubMatrix(M, eye, 1, 0);
            Vector roots = EigenSystem.EigenPair(M, QRFactorizationMethod.HouseHolder).eigenValus.GetDiagonalElements();
            //double cccond = M.Cond();
            double cdcd = M.Condest();
            mat = M;
            return roots;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double LegendrePolynomial(double x,int order)
        {
            
            if (order == 0) {return 1.0; }
            else if (order == 1) { return x;}
            else
            {
                int n = order - 1;
                double c0 = (double)((2 * n) + 1.0) / (double)(n + 1.0);
                double c1 = (double)n / ((double)(n + 1.0));
                return c0 * x * LegendrePolynomial(x, n) - c1 * LegendrePolynomial(x, n - 1);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static (Vector Root,Vector Weight) LegendreNodesWeights(int n)
        {
            double c, d, p1, p2, p3, dp;
            double[] x = new double[n];
            double[] w = new double[n];
            for (int i = 0; i < (n + 1) / 2; i++)
            {
                c = Math.Cos(Math.PI * (4 * i + 3) / (4 * n + 2));
                do
                {
                    p2 = 0;
                    p3 = 1;
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = ((2 * j + 1) * c * p2 - j * p1) / (j + 1);
                    }
                    dp = n * (c * p3 - p2) / (c * c - 1);
                    d = c;
                    c -= p3 / dp;
                }
                while (Math.Abs(c - d) > 1e-12);
                x[i] = c;
                x[n - 1 - i] = -c;
                w[i] = 2 * (1 - x[i] * x[i]) / (n + 1) / (n + 1) / LegendrePolynomial(x[i], n + 1) /
                LegendrePolynomial(x[i], n + 1);
                w[n - 1 - i] = w[i];
            }
            return (x, w);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double LaguerrePolynomial(double x, int order)
        {

            if (order == 0) { return 1.0; }
            else if (order == 1) { return 1.0-x; }
            else
            {
                int k = order - 1;
                double c0 = (double)(((2 * k) + 1.0) - x);
                double c1= (double)(k + 1.0);
                return (c0 * LaguerrePolynomial(x, k) -k* LaguerrePolynomial(x, k - 1))/c1;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static (Vector Root, Vector Weight) LaguerreNodesWeights(int n)
        {
            double c = 0.0;
            double d, p1, p2, p3, dp;
            double[] x = new double[n];
            double[] w = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    c = 3 / (1 + 2.4 * n);
                }
                else
                {
                    if (i == 1)
                    {
                        c += 15 / (1 + 2.5 * n);
                    }
                    else
                    {
                        c += (1 + 2.55 * (i - 1)) / (1.9 * (i - 1)) * (c - x[i - 2]);
                    }
                }
                do
                {
                    p2 = 0;
                    p3 = 1;
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = ((-c + 2 * j + 1) * p2 - j * p1) / (j + 1);
                    }
                    dp = (n * p3 - n * p2) / c;
                    d = c;
                    c = c - p3 / dp;
                }
                while (Math.Abs(c - d) > 1e-12);
                x[i] = c;
                w[i] = x[i] / (n + 1) / (n + 1) / LaguerrePolynomial(x[i], n + 1) / LaguerrePolynomial(x[i], n + 1);
            }
            return (x, w);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="node"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double ChebyshevPolynonial(double node, int order)
        {
            if (order == 0) { return 1.0; }
            else if (order == 1) { return node; }
            else
            {
                return 2*node * ChebyshevPolynonial(node, order - 1) - ChebyshevPolynonial(node, order - 2);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="order"></param>
        /// <returns></returns>
        public static (Vector Root,Vector Weight) ChebyshevRoots(int order)
        {
            Vector roots = new Vector();
            Vector w = new Vector();
            for (int i = 0; i < order; i++)
            {
                w.Add(PI / order);
                double theta = ((double)(2 * i + 1.0) / (2 * order)) * PI;
                double r = Cos(theta);
                roots.Add(r);
            }
            return (roots,w);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double HermitePolynomial(double x, int order)
        {

            if (order == 0) { return 1.0; }
            else if (order == 1) { return 2*x; }
            else
            {
                int k = order - 1;
                return 2 * x * HermitePolynomial(x, k) - 2*k * HermitePolynomial(x, k - 1);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static (Vector Roots, Vector Weight) HermiteNodesWeights(int n)
        {
            double c = 0.0;
            double d, p1, p2, p3, dp;
            double[] x = new double[n];
            double[] w = new double[n];
            for (int i = 0; i < (n + 1) / 2; i++)
            {
                if (i == 0)
                {
                c = Sqrt(2 * n + 1) - 1.85575 * Pow(2 * n + 1, -((double)(1)/ (double)(6)));
                }
                else
                {
                    if (i == 1)
                    {
                        c = c - 1.14 * Pow(n, 0.426) / c;
                    }
                    else
                    {
                        if (i == 2)
                        {
                            c = 1.86 * c - 0.86 * x[0];
                        }
                        else
                        {
                            if (i == 3)
                            {
                                c = 1.91 * c - 0.91 * x[1];
                            }
                            else
                            {
                                c = 2 * c - x[i - 2];
                            }
                        }
                    }
                }
                do
                {
                    p2 = 0;
                    p3 = Pow(PI, -0.25);
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = p2 * c * Sqrt((double)(2) / ((double)(j + 1))) -
                        p1 * Sqrt((double)(j) / ((double)(j + 1)));
                    }
                    dp = Sqrt(2 * n) * p2;
                    d = c;
                    c -= p3 / dp;
                }
                while (Abs(c - d) > 1e-12);
                x[i] = c;
                w[i] = Pow(2, n + 1) *SpecialFunctions.Gamma(n + 1) * Math.Sqrt(Math.PI) / HermitePolynomial(x[i], n + 1) / HermitePolynomial(x[i], n + 1);
                x[n - 1 - i] = -x[i];
                w[n - 1 - i] = w[i];
            }
            return (x, w);
        }
    }
}
