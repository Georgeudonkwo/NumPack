using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using NumSharp.Solvers;
namespace NumSharp
{
    /// <summary>
    /// 
    /// </summary>
   public class Interpolation
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static Matrix NewtonMatrix(Vector x)
        {
            Matrix N = Matrix.IdentityMatrix(x.Length);
            Vector ones = Vector.Ones(x.Length);
            N[0, "col"] = ones;
            
            for (int i = 1; i < N.Nrow; i++)
            {
               
                for (int j = 0;j< N.Nrow; j++)
                {
                    double prod = 1.0;
                    if (j >= i)
                    {
                        
                        for (int k = 0; k < i; k++)
                        {
                            prod *= (x[j] - x[k]);
                        }
                        N[j, i] = prod;
                    }
                   
                }
               
            }
            return N;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double NewtonInterpolation3(Vector xval, Vector yval, double p,int order=0)
        {
            if (order == 0) { order = xval.Length - 1; };
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Matrix Nmat = NewtonMatrix(x);
            Vector c = LinearSolvers.ForwardSubstitution(Nmat, y);
            double sum = 0.0;
            for (int k = 0; k < c.Length; k++)
            {
                double prod = 1.0;
                for (int j = 0; j <= k - 1; j++)
                {
                    prod *= (p - x[j]);
                }
                sum += c[k] * prod;
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix NewtonDividedDiffernce(Vector x, Vector y, int order = 0)
        {
            if (x.Length != y.Length) { throw new Exception("The dependent and independent variables must be the same size"); }
            if (order >= x.Length) { throw new Exception("The polynomial degree must ne less than the number of data points"); }
            Matrix newtonDiff = new Matrix();
            Vector divideDiff = new Vector(y);
            if (order == 0) { order = x.Length - 1; };

            bool exit = true;
            int i = 0;
            do
            {
                Vector nextDivideDiff = new Vector();
                newtonDiff.Add(divideDiff);
                for (int k = 0; k < divideDiff.Length ; k++)
                {
                    if (k < divideDiff.Length - 1-i)
                    {
                        nextDivideDiff.Add((divideDiff[k + 1] - divideDiff[k]) / (x[k + 1 + i] - x[k]));
                    }
                    else
                    {
                        nextDivideDiff.Add(0.0);
                    }
                   
                }

                int count = nextDivideDiff.Length;
                divideDiff = nextDivideDiff;
                //exit = count == order;
                exit = i >= order;
                i++;

            } while (!exit);
            return newtonDiff;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="yFunc"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix NewtonDividedDiffernce(Vector x, Func<double, double> yFunc, int order = 0)
        {
            Vector y =(Vector) x.Select(v => yFunc(v)) ;
            if (x.Length != y.Length) { throw new Exception("The dependent and independent variables must be the same size"); }
            if (order >= x.Length) { throw new Exception("The polynomial degree must ne less than the number of data points"); }
            Matrix newtonDiff = new Matrix();
            Vector divideDiff = new Vector(y);
            if (order == 0) { order = x.Length - 1; };

            bool exit = true;
            int i = 0;
            do
            {
                Vector nextDivideDiff = new Vector();
                newtonDiff.Add(divideDiff);
                for (int k = 0; k < divideDiff.Length; k++)
                {
                    if (k < divideDiff.Length - 1 - i)
                    {
                        nextDivideDiff.Add((divideDiff[k + 1] - divideDiff[k]) / (x[k + 1 + i] - x[k]));
                    }
                    else
                    {
                        nextDivideDiff.Add(0.0);
                    }

                }

                int count = nextDivideDiff.Length;
                divideDiff = nextDivideDiff;
                //exit = count == order;
                exit = i >= order;
                i++;

            } while (!exit);
            return newtonDiff;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double NewtonInterpolation(Vector  xval, Vector yval, double p, int order = 0)
        {
            if (order == 0) { order = xval.Length - 1; };
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Matrix dividedDiff = NewtonDividedDiffernce(x, y, order);
            double res = dividedDiff[0][0];
            for (int i = 1; i < dividedDiff.Ncol; i++)
            {
                double prod = 1.0;
                for (int k = 0; k < i; k++)
                {
                    prod *= (p - x[k]);
                }
                res += dividedDiff[0, i] * prod;
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public  static (Vector IndependentVec,Vector DependentVec)  InterpolationInterval(Vector xval, Vector yval, double p, int order)
        {
            Vector x = new Vector();
            Vector y = new Vector();
            if (p< xval[0] || p > xval[yval.Length - 1]) { order = xval.Length - 1; }
            if (order >= xval.Length ) { throw new Exception(""); }
                if (order == xval.Length - 1) { x = xval.Copy(); y = yval.Copy();return (x, y); }
            for (int i = 0; i < xval.Length; i++)
                {
                    if (xval[i] >= p )
                    {
                        if(xval[i]==p && i == 0)//loop forward
                        {
                            for (int k = 0; k <= order; k++)
                            {
                                x.Add(xval[i+k]);
                                y.Add(yval[i+k]);
                            }
                            break;
                        }
                        else if (xval[i] == p && i == xval.Length - 1)//loop backward
                        {
                            for (int k = xval.Length - 1; k >= xval.Length- order-1; k--)
                            {
                                x.Add(xval[k]);
                                y.Add(yval[k]);
                            }
                            x = x.Reverse();
                            y = y.Reverse();
                            break;
                        }
                        else if (xval[i] >= p && i != 0)//loop backward and forward
                        {
                            if (xval.Length - i >= order)//loop forward only
                            {
                                for (int k = 0; k < order; k++)
                                {
                                    x.Add(xval[i + k]);
                                    y.Add(yval[i + k]);
                                }
                                x.Prepend(xval[i - 1]);
                                y.Prepend(yval[i - 1]);
                                break;
                            }
                            else
                            {
                                int LenF= (xval.Length)-i;
                                int lenB =Abs( order - LenF)+1;
                                for (int j = 0; j < lenB; j++)
                                {
                                    x.Add(xval[i-1 - j]);
                                    y.Add(yval[i-1 - j]);
                                }
                                    x = x.Reverse();
                                    y = y.Reverse();
                            for (int k = 0; k < LenF; k++)
                                {
                                    x.Add(xval[i + k]);
                                    y.Add(yval[i + k]);
                                }
                                break;
                            }
                        }
                    }
            }

            return (x,y);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yFunc"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double NewtonInterpolation(Vector xval, Func<double,double>yFunc, double p, int order = 0)
        {
            Vector yval =(Vector) xval.Select(v => yFunc(v));
            if (order == 0) { order = xval.Length - 1; };
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Matrix dividedDiff = NewtonDividedDiffernce(x, y, order);
            double res = dividedDiff[0][0];
            for (int i = 1; i < dividedDiff.Ncol; i++)
            {
                double prod = 1.0;
                for (int k = 0; k < i; k++)
                {
                    prod *= (p - x[k]);
                }
                res += dividedDiff[0, i] * prod;
            }
            return res;
        }
        static Vector NewtonPolynomialCoefficient(Vector x, Vector y)
        {
            int n = x.Length;
            Vector coeffs = new Vector();
            for (int i = 0; i < n; i++)
            {
                if (i == 0) { coeffs.Add(y[i]); continue; }
                double sum = 0.0;
                for (int j = 0; j <= i - 1; j++)
                {
                    double prod = 1.0;
                    for (int k = 0; k <= j - 1; k++)
                    {
                        prod *= (x[i] - x[k]);
                    }
                    sum += coeffs[j] * prod;

                }
                double prod2 = 1.0;
                for (int k = 0; k <= i - 1; k++)
                {
                    prod2 *= (x[i] - x[k]);
                }
                coeffs.Add((y[i] - sum) / prod2);
            }
            return coeffs;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double NewtonInterpolation2(Vector xval, Vector yval, double p,int order=0)
        {
           
            if (order == 0) { order = xval.Length - 1; };
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Vector a = NewtonPolynomialCoefficient(x, y);
            double sum = 0.0;
            for (int k = 0; k < a.Length; k++)
            {
                double prod = 1.0;
                for (int j = 0; j <= k - 1; j++)
                {
                    prod *= (p - x[j]);
                }
                sum += a[k] * prod;
            }
            return sum;
        }
        private static Vector BaryCentricWeight(Vector xVals)
        {
            Vector res = new Vector( xVals.Length);
            for (int i = 0; i < xVals.Length; i++)
            {
                double w = 1.0;
                for (int j = 0; j < xVals.Length; j++)
                {
                    if (i != j)
                    {
                        w *= 1.0 / (xVals[i] - xVals[j]);
                    }
                }
                res[i] = w;
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="P"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static (double polyVal ,Vector Coefficients) BaryCentricInterp(Vector xval, Vector yval, double P,int order=0)
        {
            if (order == 0) { order = xval.Length - 1; }
            var TrimVectors = InterpolationInterval(xval, yval, P, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Vector w = BaryCentricWeight(x);
           Vector  coefficients = new Vector(w.Length);
            double sum = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                sum += w[i] / (P - x[i]);
            }
            double Ysum = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                Ysum += (w[i] * y[i]) / (P - x[i]);
            }
            double res = Ysum / sum;
            for (int i = 0; i < x.Length; i++)
            {
                coefficients[i] = (w[i] / (P - x[i])) / sum;
            }
            return (res,coefficients);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double LangrageInterpolation(Vector xval, Vector yval, double p,int order=0)
        {
            if (order == 0) { order = xval.Length - 1; }
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Vector CardinalPoly = CardinalPolynomial(x, p);
            double result = (CardinalPoly * y).Sum(); 
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="fval"></param>
        /// <param name="xpoint"></param>
        /// <param name="ypoint"></param>
        /// <returns></returns>
        public static double BivariateInterpolation(Vector x, Vector y, Vector fval, double xpoint, double ypoint)
        {
            
            Vector xCardinalPoly = CardinalPolynomial(x, xpoint);
            Vector ycardinalPoly = CardinalPolynomial(y, ypoint);
            double sum = 0.0;
            for (int i = 0; i < xCardinalPoly.Length; i++)
            {
                for (int j = 0; j < ycardinalPoly.Length; j++)
                {
                    sum += fval[i] * xCardinalPoly[i] * ycardinalPoly[j];
                }
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="DependentVariable"></param>
        /// <param name="xpoint"></param>
        /// <param name="ypoint"></param>
        /// <returns></returns>
        public static double BivariateInterpolation(Vector x, Vector y, Matrix DependentVariable, double xpoint, double ypoint)
        {
            //x=p,y=t
            Matrix A = new Matrix(x.Length + 1, x.Length + 1);
            A[0, 1, 0, "col"] = x;
            A[0, 1, 0, "row"] = ~y;
            Matrix.InsectSubMatrix(A, DependentVariable, 1, 1);
            int xindex = -1;
            int yindex = -1;
            bool breakOut = false;
            bool breakOut2 = false;
            int lastCol = A.Ncol - 1;
            int lastRow = A.Nrow - 1;
            for (int i = 1; i < A.Nrow; i++)
            {
                if (A[i, 0] >= xpoint)
                {
                    if (i == 1)
                    {
                        xindex = i;
                        breakOut = true;
                    }
                    else if (i == lastRow)
                    {
                        xindex = i-1;
                        breakOut = true;
                    }
                    else
                    {
                        xindex = i;
                        breakOut = true;
                    }
                }
                if (breakOut) { break; }
            }
            for (int i = 1; i < A.Nrow; i++)
            {

                
                if (A[0, i] >= ypoint)
                {
                    if (i == 1)
                    {
                        yindex = i;
                        breakOut2 = true;
                    }
                    else if (i == lastCol)
                    {
                        yindex = i-1;
                        breakOut2 = true;
                    }
                    else
                    {
                        yindex = i;
                        breakOut2 = true;
                    }
                }
                    if (breakOut2) { break; }
            }
            Matrix B = Matrix.SubMatrix(A, xindex, yindex, xindex+1, yindex+1);
            xindex = xindex - 1;yindex = yindex - 1;
            Vector C0 = Vector.Ones(B.Ncol+2);
            Matrix BiMat = new Matrix();
            BiMat.Add(C0);
            Vector C1 = new Vector() {y[yindex],y[yindex],x[xindex], x[xindex] };
            BiMat.Add(C1);
            Vector C2 = new Vector() { x[xindex], x[xindex+1], x[xindex ], x[xindex + 1] };
            BiMat.Add(C2);
            Vector C3 = new Vector() { y[yindex]* x[xindex], y[yindex]* x[xindex + 1], y[yindex+1]*x[xindex], y[yindex+1]* x[xindex + 1] };
            BiMat.Add(C3);
            Vector Rhs = new double[] {B[0,0],B[1,0],B[0,1],B[1,1] };
            Vector coefs = LinearSolvers.Solve(BiMat, Rhs);
            double result = coefs[0] + coefs[1] * ypoint + coefs[2]*xpoint + coefs[3] * ypoint * xpoint;
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="DependentVariable"></param>
        /// <param name="xpoint"></param>
        /// <param name="ypoint"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double BivariateInterpolation2(Vector x, Vector y, Matrix DependentVariable, double xpoint, double ypoint,int order=0)
        {
            double res = 0.0;
            Vector Vlist = new Vector();
            DependentVariable.Orientation = "row";
            foreach (var item in DependentVariable)
            {
                double inter1 = NewtonInterpolation(y, ~item, ypoint,order);
                Vlist.Add(inter1);
            }
            res = NewtonInterpolation(x, Vlist, xpoint,order);
            return res;
        }
        private static Vector CardinalPolynomial(Vector x, double p)
        {
            int n = x.Length;
            double l = 0.0;
            Vector CardinalPoly = new Vector();
            for (int i = 0; i < n; i++)
            {
                l = 1.0;
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                    {
                        l *= (p - x[j]) / (x[i] - x[j]);
                    }
                }
                CardinalPoly.Add(l);
            }

            return CardinalPoly;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="order"></param>
        /// <returns></returns>
        public static Vector ChebyshevNodes(int order)
        {
           Vector nodes = new Vector(order + 1);
            for (int i = 0; i < nodes.Length; i++)
            {
                double theta = ((2.0 * order) + 1.0 - (2.0 * i)) / (2.0 * (order + 1));
                nodes[i] = Cos(theta * PI);
            }
            return nodes;
        }
       
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double AitkenInterpolation(Vector xval, Vector yval, double p,int order=0)
        {
            if (order == 0) { order = xval.Length - 1; }
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Vector xx = x.Copy();Vector yy = y.Copy();
            int L = 1;
            do
            {
                Vector vec = new Vector();
                Matrix AitMat = new Matrix(2, 2);
                for (int i = 0; i < yy.Length-1; i++)
                {
                    AitMat[0, 0] = yy[0];
                    AitMat[0, 1] = xx[0] - p;
                    AitMat[1, 0] = yy[i + 1];
                    AitMat[1, 1] = xx[i + 1] - p;
                    double fval = (1.0 / (xx[i + 1] - xx[0])) * AitMat.Det();
                    vec.Add(fval);
                }
                yy = vec.Copy();
                xx = x.LastNth(L++);
            } while (yy.Length>1);
            return yy[0];
        }
        /// <summary>
        /// buugy do not use
        /// </summary>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <param name="p"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double NevilleInterpolation(Vector xval, Vector yval, double p,int order=0)
        {
            if (order == 0) { order = xval.Length - 1; }
            var TrimVectors = InterpolationInterval(xval, yval, p, order);
            Vector x = TrimVectors.IndependentVec;
            Vector y = TrimVectors.DependentVec;
            Matrix NevilleMat = new Matrix(y.Length,y.Length);
            NevilleMat.Add(y);
            for (int i = 1; i < NevilleMat.Nrow; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    NevilleMat[i, j] = (((p - x[i - j]) * NevilleMat[i, j - 1]) - ((p - x[i]) * NevilleMat[i - 1, j - 1])) / (x[i] - x[i - j]);
                }
            }
            double res = NevilleMat[NevilleMat.Nrow-1, NevilleMat.Nrow-1];
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="yval"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static double InverseInterpolation(Vector x,Vector y,double yval,OptimizationAndSolverSettings setting=null!)
        {

            int insertIndex = -1;
            for (int i = 0; i < y.Length; i++)
            {
                if (y[i] >= yval)
                {
                    insertIndex = i;
                    if (insertIndex == y.Length - 1) { insertIndex = insertIndex - 1; }
                    break;
                }
            }
            Func<double, double> fun = new Func<double, double>(c =>
                {
                    double res =yval- NewtonInterpolation2(x, y, c);
                    return res;
                });
            double initialG = x[insertIndex+1];
            double result = NonLinearSolver.NewtonFzero(fun, initialG,setting);
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="dydx"></param>
        /// <returns></returns>
        public static Matrix HermiteCoefficients(Vector x, Vector y, Vector dydx)
        {
            int n = x.Length;
            Matrix coeffs = new Matrix(2 * n, 2 * n);
            Vector z = new Vector(2 * n);
            for (int i = 0; i < n; i++)
            {
                z[2 * i] = x[i];
                z[2 * i + 1] = x[i];
                coeffs[2 * i, 0] = y[i];
                coeffs[2 * i + 1, 0] = y[i];
                coeffs[2 * i + 1, 1] = dydx[i];
                if (i != 0)
                {
                    coeffs[2 * i, 1] = (coeffs[2 * i, 0] - coeffs[2 * i - 1, 0]) / (z[2 * i] - z[2 * i - 1]);
                }
            }
            for (int i = 2; i < 2 * n; i++)
            {
                for (int j = 2; j <= i; j++)
                {
                    coeffs[i, j] = (coeffs[i, j - 1] - coeffs[i - 1, j - 1]) / (z[i] - z[i - j]);
                }
            }
            return coeffs;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="dydx"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        public static double HermiteInterpolation(Vector x, Vector y, Vector dydx, double xval)
        {
            double res = 0.0;
            Matrix coeff = HermiteCoefficients(x, y, dydx);
            Vector CD = coeff.GetDiagonalElements();
           Vector  xd = xval - x;
            res = coeff[0, 0];
            for (int i = 1; i < CD.Length; i++)
            {
                double prod = 1.0;
                int j = 0;
                for (int k = 1; k <= i; k++)
                {
                    prod *= xd[j];
                    if (k % 2 == 0) { j++; }
                }
                res += CD[i] * prod;
            }
            return res;
        }

    }
}
