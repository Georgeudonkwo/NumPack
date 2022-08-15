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
  public  class EigenSystem
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (double lamda, Vector vec) PowerMethod(Matrix A)
        {
            int dominantDiaColumn(Matrix B)
            {
                double max = Abs(B[0, 0]);
                int maxColIndex = 0;
                for (int i = 0; i < B.Nrow; i++)
                {
                    for (int j = 0; j < B.Ncol; j++)
                    {
                        if (i == j)
                        {
                            if (Abs(B[i, j]) > max)
                            {
                                max = Abs(B[i, j]);
                                maxColIndex = j;
                            }
                            break;
                        }
                    }
                }
                return maxColIndex;
            }
            int k = dominantDiaColumn(A);
            Vector v = A[k, "col"];//  ExtractColumnWise(A, 0, k);
            Vector y = v / v.Norm();// ScalarVectorMult(1.0 / Norm(v), v);
            double Lamnda = 0.0;
            double eps = 1e-6;
            bool exit = false;
            int iter = 0;
            do
            {
                iter++;
                Vector ynext = A * y;//  MatrixVectorMult(A, y);
                double Lamndanext = y.Dot(ynext);// DotProduct(y, ynext);
                y = ynext.Copy();
                y = y / y.Norm();
                exit = Abs(Lamnda - Lamndanext) >= eps;
                Lamnda = Lamndanext;
                if (iter > 500) break;
            } while (exit);
            return (Lamnda, y);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="shift"></param>
        /// <returns></returns>
        public static (double MinLamda, Vector EigVector) InversePowerMethod(Matrix A, double shift = 0.0)
        {
            Matrix eye = Matrix. IdentityMatrix(A.Nrow);
            Matrix B = A.Copy();
            B = B - shift * eye;
            Vector vecini = B[0, "col"];// ExtractColumnWise(A, 0, 0);
            var lup = MatrixFactorization.LUP(B);
            for (int i = 0; i < vecini.ColVector.Length; i++)
            {
                vecini[i] = vecini[i] / vecini.Norm();
            }
            bool exit = false;
            int iter = 0;
            Vector v = vecini.Copy();
            double Nnorm = 0;
            double rayleighQuotient = 0.0;
            Vector q = vecini / vecini.Norm();// vecini.Select(u => u / Norm(vecini)).ToArray();
            Vector q2 = new Vector(q.ColVector.Length);
            q2 = q.Copy();
            do
            {
                iter++;

               Vector b =lup.P* q;
                Vector y = LinearSolvers.ForwardSubstitution(lup.L, b);
               Vector z = LinearSolvers.BackSubstition(lup.U, y);
                double zNorm = z.Norm();
                q = z / z.Norm();//  z.Select(item => item / zNorm).ToArray();
                z =A* q;
                rayleighQuotient = q.Dot( z);
                b = q2.Copy();
                y = LinearSolvers.ForwardSubstitution(Matrix. Transpose(lup.U), b);
                Vector w = LinearSolvers.BackSubstition(Matrix. Transpose(lup.L), y);

                Vector Av = B* z;
                rayleighQuotient = z.Dot(Av) / v.Norm();
                Vector LanV = rayleighQuotient* z;
                Vector diff = Av - LanV;// Av.Select((a, i) => a - LanV[i]).ToArray();

                Nnorm = diff.Norm() / v.Norm();
                exit = Nnorm > 1e-8;
                if (iter > 500) break;
            } while (exit);
            return (rayleighQuotient, v);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (double MinLamda, Vector EigVector) InversePowerMethod(Matrix A)
        {
            Matrix eye = Matrix. IdentityMatrix(A.Nrow);
            Matrix B = A.Copy();
            int dominantDiaColumnIndex(Matrix C)
            {
                double max = Abs(C[0, 0]);
                int maxColIndex = 0;
                for (int i = 0; i < C.Nrow; i++)
                {
                    for (int j = 0; j < C.Ncol; j++)
                    {
                        if (i == j)
                        {
                            if (Abs(C[i, j]) > max)
                            {
                                max = Abs(C[i, j]);
                                maxColIndex = j;
                            }
                            break;
                        }
                    }
                }
                return maxColIndex;
            }
            int k = dominantDiaColumnIndex(A);
            Vector vecini = B[k, "col"];
               var Lup=MatrixFactorization.LUP(B);
            vecini = vecini / vecini.Norm();
            bool exit = false;
            int iter = 0;
           Vector v = vecini.Copy();
            double Nnorm = 0;
            double rayleighQuotient = 0.0;
            do
            {
                iter++;
                v =Lup.P* v;
                Vector c = LinearSolvers.ForwardSubstitution(Lup.L, v);
                Vector vnew = LinearSolvers.BackSubstition(Lup.U, c);
                vnew = vnew / vnew.Norm();
                v = vnew.Copy();

                Vector Av = A * v;
                rayleighQuotient = v.Dot(Av) / v.Norm();
                Vector LanV = rayleighQuotient * v;
                Vector diff = Av - LanV;

                Nnorm = diff.Norm() / v.Norm();
                exit = Nnorm > 1e-8;
                if (iter > 500) break;
            } while (exit);
            return (rayleighQuotient, v);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="qrScheme"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        public static (Matrix eigenValus, Matrix eigenVectors) EigenPair(Matrix A, QRFactorizationMethod qrScheme,double tol=1e-8)
        {
            Matrix B = Matrix. Hessenberg(A);
            // double[,] Cmat = Duplicate(A);
           Vector calculate2by2MatEigVal(Matrix C)
            {
                double a = 1.0;
                double b = -1 * (C[0, 0] + C[1, 1]);
                double c = (C[0, 0] * C[1, 1]) - (C[0, 1] * C[1, 0]);
                double lamda1 = (-b + Sqrt(Pow(b, 2) - 4 * a * c)) / (2 * a);
                double lamda2 = (-b - Sqrt(Pow(b, 2) - 4 * a * c)) / (2 * a);
                Vector res = new Vector(2);
                res[0] = lamda1;
                res[1] = lamda2;
                return res;
            }
            double WilkinsonShift()
            {

                int row = B.Nrow;
                int col = B.Ncol;
                Matrix mat2by2 = new Matrix(2, 2);
                mat2by2[0, 0] = B[row - 2, col - 2];
                mat2by2[0, 1] = B[row - 2, col - 1];
                mat2by2[1, 0] = B[row - 1, col - 2];
                mat2by2[1, 1] = B[row - 1, col - 1];
                Vector vals = calculate2by2MatEigVal(mat2by2);
                double lamda = (vals[0] - mat2by2[1, 1]) <= (vals[1] - mat2by2[1, 1]) ? vals[0] : vals[1];
                if (double.IsNaN(lamda) || double.IsInfinity(lamda)) { lamda = mat2by2[1, 1]; }
                return lamda;
            }
            Vector Lamdas = new Vector();
            while (B.Ncol >= 2)
            {
                Matrix eye = Matrix. IdentityMatrix(B.Nrow);
                double s = WilkinsonShift();
                bool exit = false;
                int iter = 0;

                do
                {
                    iter++;
                    B = B- s* eye;
                   var qr=MatrixFactorization.QR(B, qrScheme);
                    B = qr.R* qr.Q+ s* eye;
                    int r = B.Nrow - 1;
                    int c = B.Ncol - 2;
                    exit = Abs(B[r, c]) < tol;
                    double lam = 0.0;
                    if (exit)
                    {
                        lam = B[B.Nrow - 1, B.Ncol - 1];
                        Lamdas.Add(lam);
                    }
                    if (iter > 500)
                    {
                        Lamdas.Add(lam);
                        break;
                    }

                } while (!exit);
                B = deflate(B);
                Matrix deflate(Matrix C)
                {
                    Matrix DefMat = new Matrix(C.Nrow - 1, C.Ncol - 1);
                    for (int i = 0; i < DefMat.Nrow; i++)
                    {
                        for (int j = 0; j < DefMat.Ncol; j++)
                        {
                            DefMat[i, j] = C[i, j];
                        }
                    }
                    return DefMat;
                }
            }
            Lamdas.Add(B[0, 0]);
            Lamdas = Lamdas.OrderByDescending(x => x).ToList();
            Matrix eigVals = Matrix. IdentityMatrix(A.Nrow);
            Matrix eigVec = new Matrix(A.Nrow, A.Ncol);
            calculateEigenValuesAndEigenVectors();
            void calculateEigenValuesAndEigenVectors()
            {
                int o = 0;
                for (int i = 0; i < eigVals.Nrow; i++)
                {
                    for (int j = 0; j < eigVals.Ncol; j++)
                    {
                        if (i == j) { eigVals[i, j] = Lamdas[o]; }
                    }
                    o++;
                }

                for (int i = 0; i < eigVals.Nrow; i++)
                {
                    Matrix eye = Matrix. IdentityMatrix(eigVec.Nrow);

                    for (int j = 0; j < eigVals.Ncol; j++)
                    {

                        if (i == j)
                        {
                            double land = eigVals[i, j];
                            Matrix va = A- land* eye;
                            Vector v2 = va[i, "col"];
                            Vector vecx = EigenVectorsCalc2(A, land);
                            eigVec[i, "col"] = vecx;
                        }
                    }
                }
            }
              Vector EigenVectorsCalc2(Matrix D, double lamda)
            {
                if (lamda == 0.0) { lamda = 1e-16; }
                double sigma = lamda;// - 0.05;
                Matrix eye = Matrix. IdentityMatrix(D.Nrow);
                Matrix E =D- sigma* eye;
                Vector vecini = D[0, "col"];
               var Lup=MatrixFactorization.LUP(E);
                double pertube = 0.0001;int count = 0;
                while (Lup.U.GetDiagonalElements().Contains(0.0)&&count<=10)
                {
                    count++;
                    sigma = lamda - pertube;
                    E = D - sigma * eye;
                    Lup = MatrixFactorization.LUP(E);
                    pertube += 0.05;
                }
                vecini = vecini / vecini.Norm();
                bool exit = false;
                int iter = 0;
               Vector v = vecini.Copy();
                double Nnorm = 0;
                double rayleighQuotient = 0.0;
                do
                {
                    iter++;
                    v = Lup.P* v;
                    Vector c = LinearSolvers.ForwardSubstitution(Lup.L, v);
                    Vector vnew = LinearSolvers.BackSubstition(Lup.U, c);
                    //Vector vnew = LinearSolvers.Solve(E,v);
                    vnew = vnew / vnew.Norm();
                    v = vnew.Copy();

                    Vector Av = D* v;
                    rayleighQuotient = v.Dot(Av) / v.Norm();
                  Vector LanV = rayleighQuotient* v;
                    Vector diff = Av - LanV;//  Av.Select((a, i) => a - LanV[i]).ToArray();
                    Nnorm = diff.Norm() / v.Norm();
                    exit = Nnorm > 1e-8;// Abs(( lam-lamda))> 1e-8;
                    if (iter > 500) break;
                } while (exit);
                return v;
            }
            return (eigVals, eigVec);
        }
    }

}
