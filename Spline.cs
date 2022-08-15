using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumSharp
{
    /// <summary>
    /// 
    /// </summary>
    public class Spline
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public Vector CubicSpline1(Vector x, Vector y)
        {
            Vector h = new Vector();
            Vector La = new Vector();
            Vector mue = new Vector();
            Matrix Mat = Matrix.Diagonal(x.Length, 2);
            Vector d = new Vector();
            for (int i = 1; i < x.Length; i++)
            {
                h.Add((x[i] - x[i - 1]));
            }
            for (int i = 0; i < h.Length - 1; i++)
            {
                La.Add((h[i + 1] / (h[i] + h[i + 1])));
                mue.Add(1.0 - h[i]);
            }
            for (int i = 1; i < y.Length - 1; i++)
            {
                double T1 = (y[i + 1] - y[i]) / h[i];
                double T2 = (y[i] - y[i - 1]) / h[i - 1];
                double term = (T1 - T2) / (h[i - 1] + h[i]);
                d.Add(6 * term);
            }
            int n = y.Length - 1;
            La.Append(1);
            mue.Prepend(1);
            d.Append((3 * (y[1] - y[0]) / h[1]));
            d.Append(3 * (y[n] - y[n - 1]) / h[h.Length - 1]);
            //form the tridiagonal matrix
            for (int i = 0; i < Mat.Nrow; i++)
            {
                if (i == 0) { Mat[i, i + 1] = mue.First();continue; }
                else if (i == Mat.Nrow - 1) { Mat[i, i - 1] = mue.Last();continue; }
                for (int j = 0; j < Mat.Ncol; j++)
                {
                    if (i == j)
                    {
                        Mat[i, j - 1] = mue[i - 1];
                        Mat[i, j + 1] = La[i - 1];
                        break;
                    }
                }
            }
            Vector theta = Solvers.LinearSolvers.SolveTriad(Mat, d);
            return theta;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public Vector CubicSpline(Vector x, Vector y)
        {
            Vector h = new Vector();
            Matrix Mat = Matrix.Zeros(x.Length-2,x.Length);
            Vector d = new Vector();
            Vector a = new Vector(y);
            Vector adiff = new Vector();
            for (int i = 1; i < x.Length; i++)
            {
                h.Add((x[i] - x[i - 1]));
            }
            for (int i = 1; i < a.Length; i++)
            {
                adiff.Add((a[i] - a[i - 1]));
            }
            Vector bb = new Vector();
            for (int i = 1; i < 3; i++)
            {
                bb.Add(3 * (a[i + 1] * h[i - 1] - a[i] * (x[i + 1] - x[i - 1]) +
                a[i - 1] * h[i]) / (h[i] * h[i - 1]));
            }
            int n = 4;
            double[] XA = new double[n + 1];
            double[] XL = new double[n + 1];
            double[] XU = new double[n + 1];
            double[] XZ = new double[n + 1];
            XL[0] = 1; XU[0] = 0; XZ[0] = 0;
            for ( int i = 1; i < 3; i++)
            {
                XL[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * XU[i - 1];
                XU[i] = h[i] / XL[i];
                XZ[i] = (XA[i] - h[i - 1] * XZ[i - 1]) / XL[i];
            }
            for (int i = 1; i < a.Length - 1; i++)
            {
                double T1 = (a[i + 1] - a[i]) / h[i];
                double T2 = (a[i] - a[i - 1]) / h[i - 1];
                double term = (3.0 * T1) - (3.0 * T2);
                d.Add(term);
            }
            //int n = y.Length - 1;
            d.Prepend(0.0);
            d.Append(0.0);
            //form the tridiagonal matrix
            Vector R1 = Vector.Zeros(d.Length).GetRowVector(); R1[0] = 1.0;
            Vector Rn = Vector.Zeros(d.Length).GetRowVector(); Rn[Rn.Length - 1] = 1.0;
            Mat = Mat.PrePend(R1);
            Mat = Mat.Apend(Rn);
            for (int i = 1; i < Mat.Nrow-1; i++)
            {
                for (int j = 0; j < Mat.Ncol; j++)
                {
                    if (i == j)
                    {
                        Mat[i, j] = 2 * (h[i-1] + h[i]);
                        Mat[i, j - 1] = h[i-1];
                        Mat[i, j + 1] = h[i];
                        break;
                    }
                   
                }
            }
            
            Vector theta = Solvers.LinearSolvers.SolveTriad(Mat, d);
           // Vector theta2 = Solvers.LinearSolvers.Solve(Mat, d);
            return theta;
        }
    }
}
