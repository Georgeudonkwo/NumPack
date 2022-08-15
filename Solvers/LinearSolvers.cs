using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
namespace NumSharp.Solvers
{
    /// <summary>
    /// 
    /// </summary>
 public   class LinearSolvers
    {

        internal static Vector BackSubstition(Matrix A, Vector b)
        {
            Vector x = new Vector(b.ColVector.Length);
            for (int i = A.Nrow-1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int j = i + 1; j < A.Ncol; j++)
                {
                    sum += A[i, j] * x[j];
                }               
                double divisor = A[i, i];
                x[i] = (b[i] - sum) / divisor;
            }
            return x;
        }
        internal static Vector ForwardSubstitution(Matrix A, Vector b)
        {
            Vector c = new Vector(b.ColVector.Length);
            for (int i = 0; i < A.Nrow; i++)
            {
                c[i] = b[i] / A[i, i];
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sum += A[i, j] * c[j];
                    c[i] = (b[i] - sum) / A[i, i];
                }
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector SolveLU(Matrix A, Vector b)
        {
            try
            {
                Matrix Aclone = A.Copy();
                var lup = MatrixFactorization.LUP(Aclone);
                Vector uDia = Matrix.GetDiagonalElements(lup.U);
                if  (uDia.Contains(0.0)) { throw new Exception(); }
                Vector x = lup.P * b;
                Vector c = ForwardSubstitution(lup.L, x);
                Vector res = BackSubstition(lup.U, c);
                bool isnanorifinity = res.Any(d => double.IsNaN(d) || double.IsInfinity(d));
                if (isnanorifinity) { throw new Exception(); }
                return res;
            }
            catch (Exception)
            {

                return SolveQR(A, b,QRFactorizationMethod.HouseHolder);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        
        public static Vector SolveLDLt(Matrix A, Vector b)
        {
            try
            {
                Matrix Aclone = A.Copy();
                var ldl = MatrixFactorization.LDLt(Aclone);
                Vector z = ForwardSubstitution(ldl.L, b);
                for (int i = 0; i < ldl.D.Nrow; i++)
                {
                    for (int j = 0; j < ldl.D.Ncol; j++)
                    {
                        if (i == j)
                        {
                            ldl.D[i, j] = Pow(ldl.D[i, j], -1);
                        }

                    }
                }
                Vector y = ldl.D * z;
                Vector x = BackSubstition(ldl.Lt, y);
                return x;
            }
            catch (Exception)
            {

                return SolveLU(A, b);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="QRScheme"></param>
        /// <returns></returns>
        public static Vector SolveQR(Matrix A, Vector b, QRFactorizationMethod QRScheme=QRFactorizationMethod.HouseHolder)
        {
            try
            {
                Vector res = new Vector(b.Length);
                Matrix Aclone = A.Copy();
                var qr = MatrixFactorization.QR(Aclone, QRScheme);
                int row = Aclone.Nrow;
                int col = Aclone.Ncol;
                int diff = row - col;
                Vector x = new Vector(b.Length);

                if (diff != 0 && QRScheme != QRFactorizationMethod.MGS)
                {
                    Matrix QMod = new Matrix(qr.Q.Nrow, qr.Q.Ncol - diff);
                    Matrix RMod = new Matrix(qr.R.Nrow - diff, qr.R.Ncol);
                    for (int i = 0; i < qr.Q.Nrow; i++)
                    {
                        for (int j = 0; j < qr.Q.Ncol - diff; j++)
                        {
                            QMod[i, j] = qr.Q[i, j];
                        }
                    }
                    for (int i = 0; i < qr.R.Nrow - diff; i++)
                    {
                        for (int j = 0; j < qr.R.Ncol; j++)
                        {
                            RMod[i, j] = qr.R[i, j];
                        }
                    }
                    if (qr.R.Nrow == 1 && qr.R.Ncol == 1)
                    {
                        if (qr.R[0, 0] == 0.0) { qr.R[0, 0] = 1.0; }
                    }
                    x = Matrix.Transpose(QMod) * b;
                    res = BackSubstition(RMod, x);
                }
                else
                {
                    if (qr.R.Nrow == 1 && qr.R.Ncol == 1)
                    {
                        if (qr.R[0, 0] == 0.0) { qr.R[0, 0] = 1.0; }
                    }

                    x = Matrix.Transpose(qr.Q) * b;
                    var maam = qr.R * x;
                    res = BackSubstition(qr.R, x);
                }
                if (res.Any(r => double.IsInfinity(r) || double.IsNaN(r))) { throw new Exception(); }
                return res;
            }
            catch (Exception)
            {
                return SolveSVD(A, b);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector SolveChol(Matrix A,Vector b)
        {
            try
            {
                Matrix Aclone = A.Copy();
                var cho = MatrixFactorization.Cholesky(Aclone);
                Vector c = ForwardSubstitution(cho.RT, b);
                Vector res = BackSubstition(cho.R, c);
                return res;
            }
            catch (Exception)
            {
                return SolveLDLt(A, b);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector SolveTriad(Matrix A, Vector b)
        {
            var LuMat =MatrixFactorization.TriDiagonalLU(A);
           Vector L = LuMat.LowerDiag;
            Vector u = LuMat.UperDiag;
            Vector superd = LuMat.SuperDiag;
            Vector y = new Vector(b.Length);
            Vector x = new Vector(b.Length);
            y[0] = b[0];
            for (int i = 1; i < y.ColVector.Length; i++)
            {
                y[i] = b[i] - (L[i - 1] * y[i - 1]);
            }
            x[x.Length - 1] = y[y.ColVector.Length - 1] / u[u.ColVector.Length - 1];
            for (int i = x.ColVector.Length - 2; i >= 0; i--)
            {
                x[i] = (y[i] - (superd[i] * x[i + 1])) / u[i];
            }
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector Solve(Matrix A,Vector b)
        {
            if (A.Nrow != A.Ncol)
            {
                return SolveSVD(A, b);
            }
            else
            {
                if (Matrix.IsSymmetric(A))
                {
                    return SolveChol(A, b);
                }
                else
                {
                    return SolveLU(A, b);
                }
            }
           
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector SolveSVD(Matrix A,Vector b)
        {
            var svd = MatrixFactorization.SVD(A);
            Matrix U = svd.LeftSingularVec;
            Matrix S = svd.SingularValues;
            Matrix V = svd.RightSingularVec;
            Vector c = Matrix.Transpose(U) * b;

            Matrix sD = S;
            double maxVal = Matrix.GetDiagonalElements(S).Max(m => Abs(m));
            double bound= 1e-6 * maxVal;
            for (int i = 0; i < sD.Nrow; i++)
            {
                for (int j = 0; j < sD.Ncol; j++)
                {
                    if (i == j)
                    {
                        if (S[i, j] > bound)
                        {
                            sD[i, j] = 1.0 / S[i, j];
                        }
                        else
                        {
                            sD[i, j] = 0.0;
                        }
                    }
                }
            }
            if (A.Nrow > A.Ncol)
            {
                Matrix utrim = new Matrix(A.Ncol, A.Ncol);
                Matrix sDtrim = new Matrix(A.Ncol, A.Ncol);
                Vector Ctrim = new Vector(sDtrim.Nrow);
                for (int i = 0; i < utrim.Nrow; i++)
                {
                    for (int j = 0; j < utrim.Ncol; j++)
                    {
                        utrim[i, j] = U[i, j];
                        sDtrim[i, j] = sD[i, j];
                        Ctrim[i] = c[i];
                    }
                }
                Vector xtrim = V * sDtrim * Ctrim;
                return xtrim;
            }
          
            Vector x = V *sD * c;
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="QRScheme"></param>
        /// <returns></returns>
        public static (Vector OptimalVector,double residual) LeastSquare(Matrix A, Vector b, QRFactorizationMethod QRScheme = QRFactorizationMethod.HouseHolder)
        {
            double residual = 0.0;
            Vector res = new Vector(b.Length);
            try
            {
                Matrix Aclone = A.Copy();
                var qr = MatrixFactorization.QR(Aclone, QRScheme);
                int row = Aclone.Nrow;
                int col = Aclone.Ncol;
                int diff = row - col;
                Vector x = new Vector(b.Length);

                if (diff != 0 && QRScheme != QRFactorizationMethod.MGS)
                {
                    Matrix QMod = new Matrix(qr.Q.Nrow, qr.Q.Ncol - diff);
                    Matrix RMod = new Matrix(qr.R.Nrow - diff, qr.R.Ncol);
                    for (int i = 0; i < qr.Q.Nrow; i++)
                    {
                        for (int j = 0; j < qr.Q.Ncol - diff; j++)
                        {
                            QMod[i, j] = qr.Q[i, j];
                        }
                    }
                    for (int i = 0; i < qr.R.Nrow - diff; i++)
                    {
                        for (int j = 0; j < qr.R.Ncol; j++)
                        {
                            RMod[i, j] = qr.R[i, j];
                        }
                    }
                    if (qr.R.Nrow == 1 && qr.R.Ncol == 1)
                    {
                        if (qr.R[0, 0] == 0.0) { qr.R[0, 0] = 1.0; }
                    }
                    x = Matrix.Transpose(qr.Q) * b;
                    res = BackSubstition(RMod, x);
                    res = res.TakeWhile((d, i) => i < col).ToArray();
                    Vector vec = x[diff + 1, "col"];
                    residual = vec.Norm();
                }
                else
                {
                    if (qr.R.Nrow == 1 && qr.R.Ncol == 1)
                    {
                        if (qr.R[0, 0] == 0.0) { qr.R[0, 0] = 1.0; }
                    }

                    x = Matrix.Transpose(qr.Q) * b;
                    var maam = qr.R * x;
                    res = BackSubstition(qr.R, x);
                }

                return (res, residual); ;
            }
            catch (Exception)
            {
                res = SolveSVD(A, b);
                return (res, 0.0);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="InitialGuess"></param>
        /// <param name="specIndex"></param>
        /// <param name="qrScheme"></param>
        /// <param name="SpecificationVar"></param>
        /// <param name="normalized"></param>
        /// <returns></returns>
        public static Vector NullSpace(Matrix A, Vector InitialGuess, int specIndex, QRFactorizationMethod qrScheme, double SpecificationVar = 1.0, bool normalized = false)
        {
            Matrix B = new Matrix(A.Nrow - 1, A.Ncol - 1);
            int k = 0;
            for (int i = 0; i < A.Nrow; i++)
            {
                if (i != specIndex)
                {
                    int o = 0;
                    for (int j = 0; j < A.Ncol; j++)
                    {
                        if (j != specIndex)
                        {
                            B[k, o] = A[i, j];
                            o++;
                        }
                    }
                    k++;
                }
            }
           Vector vec = new Vector();
            for (int i = 0; i < InitialGuess.Length; i++)
            {
                if (i != specIndex)
                {
                    vec.Add(InitialGuess[i] * SpecificationVar);
                }
            }
            Vector TrimVec = vec.Copy();
           Vector x = Solve(B, TrimVec);
            Vector res = new Vector(x.Length + 1);
            //res[specIndex] =Sign(V[specIndex])* SpecificationVar;
            res[specIndex] = SpecificationVar;
            for (int i = 0; i < res.Length; i++)
            {
                if (specIndex == 0)
                {
                    if (i == specIndex) continue;
                    res[i] = x[i - (specIndex + 1)];
                }
                else
                {
                    if (i == specIndex) continue;
                    if (i < specIndex)
                    {
                        res[i] = x[i];
                    }
                    else
                    {
                        res[i] = x[i - 1];
                    }
                }
            }
            if (normalized)
            {
                res = res/res.Norm();
            }
            return res;
        }

        #region iterative solvers
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="x0"></param>
        /// <param name="b"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        public static Vector Gauss_Seidel(Matrix A,Vector b,  Vector x0,double tol = 1e-8)
        {
            Matrix mat = A.Copy();
            Vector x = new Vector(x0);
            bool conv = true;
            int iter = 0;
            do
            {
                iter++;
                for (int i = 0; i < mat.Nrow; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < mat.Ncol; j++)
                    {
                        if (i != j) sum += mat[i, j] * x[j];
                    }
                    x[i] = (b[i] - sum)/mat[i,i];
                }
                
                conv = (x - x0).Norm() > tol ;
                x0 = new Vector(x);
                if (iter > 1000) throw new Exception("This is a non convergent system");
            } while (conv);
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="x0"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Vector Jacobi(Matrix A, Vector b, Vector x0, double tol = 1e-8)
        {
            Matrix mat = A.Copy();
            Vector x = new Vector(x0);
            bool conv = true;
            int iter = 0;
            Matrix eye = Matrix.IdentityMatrix(A.Nrow);
            Matrix iterationMatrix = (eye - mat);
            //var eig=EigenSystem.EigenPair(iterationMatrix,QRFactorizationMethod.HouseHolder);
            //double spectralradius=eig.eigenValus.GetDiagonalElements().Max(x=>Math.Abs(x));
            //if (spectralradius > 1.0) throw new Exception("This is a non convergent system");
            do
            {
                iter++;
                x = b + iterationMatrix * x;
                conv = (x - x0).Norm() > tol;
                x0 = new Vector(x);
                if (iter > 1000) throw new Exception("This is a non convergent system");
            } while (conv);
            return x;
        }
        #endregion
    }
}
