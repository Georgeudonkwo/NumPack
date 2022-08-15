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
    public enum QRFactorizationMethod {
        /// <summary>
        /// 
        /// </summary>
        MGS,
        /// <summary>
        /// 
        /// </summary>
        HouseHolder,
        /// <summary>
        /// 
        /// </summary>
        Givens }
    /// <summary>
    /// 
    /// </summary>
    public class MatrixFactorization
    {

        #region LU Factorization
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="PermMatrix"></param>
        /// <param name="colIndex"></param>
        static void pivoting(Matrix A, Matrix PermMatrix, int colIndex)
        {

            int maxValueRowIndex = colIndex;
            double maxVal = Abs(A[colIndex, colIndex]);
            for (int j = colIndex; j < A.Nrow - 1; j++)
            {
                if (!(maxVal >= Abs(A[j + 1, colIndex])))
                {
                    maxVal = Abs(A[j + 1, colIndex]);
                    maxValueRowIndex = j + 1;
                }
            }
            if (colIndex == maxValueRowIndex) return;

            for (int j = 0; j < A.Ncol; j++)
            {
                double temp = A[maxValueRowIndex, j];
                A[maxValueRowIndex, j] = A[colIndex, j];
                A[colIndex, j] = temp;

                double tempPerm = PermMatrix[maxValueRowIndex, j];
                PermMatrix[maxValueRowIndex, j] = PermMatrix[colIndex, j];
                PermMatrix[colIndex, j] = tempPerm;
            }
            NumberOfRowInterChangeDuringLU++;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="colIndex"></param>
        static void pivoting(Matrix A, Vector b, int colIndex)
        {
            int maxValueRowIndex = colIndex;
            maxValueRowIndex = Matrix.MaxValIndex(A[colIndex, "col"]);
            if (colIndex == maxValueRowIndex) return;
            A.SwapRowVectors(colIndex, maxValueRowIndex);
            double item1 = b[colIndex];
            double item2 = b[maxValueRowIndex];
            b.Swap(ref item1, ref item2);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="P"></param>
        private static void Elimination(Matrix A, Matrix P)
        {
            for (int i = 0; i < A.Ncol; i++)
            {
                if (i < A.Ncol - 1)
                {
                    pivoting(A, P, i);
                }

                for (int j = i + 1; j < A.Nrow; j++)
                {

                    double pivot = A[i, i];
                    double lamda = A[j, i] / pivot;
                    if (double.IsInfinity(lamda) || double.IsNaN(lamda)) { lamda = 0.0; }
                    for (int k = i; k < A.Ncol; k++)
                    {
                        A[j, k] = A[j, k] - lamda * A[i, k];
                    }
                    A[j, i] = lamda;//replaces the zeros with lamda
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (Matrix L, Matrix U, Matrix P,int NumberOfRowExchange) LUP(Matrix A)
        {
            NumberOfRowInterChangeDuringLU = 0;
            Matrix Acopy = A.Copy();
            Matrix PermutationMatrix = Matrix.IdentityMatrix(Acopy.Nrow, Acopy.Ncol);//initializes the permutation matrix
            Elimination(Acopy, PermutationMatrix);
            var Lu = MakeLU(Acopy);//the permutation matrix is made during the partial pivoting process
            return (Lu.L, Lu.U, PermutationMatrix,NumberOfRowInterChangeDuringLU);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (Vector LowerDiag, Vector UperDiag, Vector SuperDiag) TriDiagonalLU(Matrix A)
        {
            Matrix Acopy = A.Copy();
            Vector upper = new Vector(Acopy.Nrow);
            Vector lower = new Vector(Acopy.Nrow - 1);
            Vector superD = new Vector(Acopy.Nrow - 1);
           upper[0] = Acopy[0, 0];
            for (int i = 0; i < upper.ColVector.Length - 1; i++)
            {
                for (int j = i; j < Acopy.Ncol - 1; j++)
                {
                    double ai = Acopy[i + 1, j];
                    double ci = Acopy[i, j + 1];
                    superD[i] = ci;
                    double bnext = Acopy[i + 1, j + 1];
                    lower[i] = ai / upper[i];
                    upper[i + 1] = bnext - (lower[i] * ci);
                    break;
                }
            }
            return (lower, upper, superD);
        }
        private static (Matrix L, Matrix U) MakeLU(Matrix A)
        {
            Matrix LTM = new Matrix(A.Nrow, A.Ncol);
            Matrix UTM = new Matrix(A.Nrow, A.Ncol);
            for (int i = 0; i < A.Nrow; i++)
            {
                for (int j = 0; j < A.Ncol; j++)
                {
                    if (i == j)
                    {
                        LTM[i, j] = 1.0;
                        UTM[i, j] = A[i, j];
                    }
                    else if (i > j)
                    {
                        LTM[i, j] = A[i, j];
                    }
                    else
                    {
                        UTM[i, j] = A[i, j];
                    }
                }
            }
            return (LTM, UTM);
        }
        #endregion
        #region Cholesky Factorization
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static (Matrix R, Matrix RT)Cholesky(Matrix A)
        {
            Matrix B = A.Copy();
            if (!Matrix.IsSymmetric(B)) { throw new Exception("Matrix must be symmetric"); }

            Matrix L = new Matrix(B.Nrow, B.Ncol);
            for (int i = 0; i < L.Nrow; i++)
            {
                double sum = 0.0;
                for (int k = 0; k <= i - 1; k++)
                {
                    sum += Pow(L[i, k], 2);
                }
                if (B[i, i] < sum) { throw new Exception("Matrix is not Positive Definite"); }
                L[i, i] = Sqrt(B[i, i] - sum);
                for (int j = i + 1; j < B.Nrow; j++)
                {
                    double sum2 = 0.0;
                    for (int k = 0; k <= i - 1; k++)
                    {
                        sum2 += L[i, k] * L[j, k];
                    }
                    if (L[i, i] == 0.0) { throw new Exception("Matrix is not Positive Definite"); }
                    L[j, i] = (B[i, j] - sum2) / L[i, i];
                }
            }
            //Matrix R = Matrix.Transpose(L);// Transpose(L);
            return (~L, L);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (Matrix L, Matrix D, Matrix Lt)LDM(Matrix A)
        {
            var luu = LUP(A);
            Matrix ll = luu.L;
            Matrix dd = Matrix.Diagonal(luu.U);
            for (int i = 0; i < dd.Nrow; i++)
            {
                for (int j = 0; j < dd.Ncol; j++)
                {
                    if (i == j)
                    {
                        dd[i, j] = Pow(dd[i, j], -1);
                    }
                  
                }
            }
            Matrix llT =  dd * luu.U;
            Matrix fff = luu.L * luu.U;
            return (ll, dd, llT);
          
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (Matrix L, Matrix D, Matrix Lt) LDLt(Matrix A)
        {
            Matrix D = Matrix.Zeros(A.Nrow, A.Ncol);
            Matrix L = Matrix.IdentityMatrix(A.Nrow, A.Ncol);
            Matrix Lt = Matrix.Zeros(A.Nrow, A.Ncol);

            for (int j = 0; j < D.Nrow; j++)
            {
                double sum1 = 0.0;
                for (int k = 0; k <= j - 1; k++)
                {
                    sum1 += Pow(L[j, k], 2) * D[k, k];
                }
                D[j, j] = A[j, j] - sum1;
                double sum2 = 0.0;
                for (int i = j+1; i < D.Nrow; i++)
                { 
                  
                    for (int k = 0; k <= j - 1; k++)
                    {
                        sum2 += L[i, k] * D[k, k] * L[j, k];
                    }
                    L[i, j] = (A[i, j] - sum2) / D[j, j];
                }
               
            }
            Lt = Matrix.Transpose(L);
            return (L, D, Lt);
        }
        /// <summary>
        /// Implements the GMW Algorithm
        /// </summary>
        /// <param name="A"></param>
        /// <param name="ApplyPivot"></param>
        /// <returns></returns>
        public static (Matrix L, Matrix D, Matrix E, Matrix P) ModifiedCholeskyDecomp(Matrix A, bool ApplyPivot = false)
        {
            
                double eps = new OptimizationAndSolverSettings().Eps;
                double MaxDiagEle = A.GetDiagonalElements().Max(m => Abs(m));
                double maxoffdiag = A.MaxOffDiagonalElement();
                Vector v1 = new double[] { maxoffdiag + MaxDiagEle, 1 };
            double delta = eps * v1.Max(x => Abs(x));
                Vector v2 = new double[] { maxoffdiag, MaxDiagEle / Sqrt(A.Ncol - 1), eps };
                v2 = v2.Pow(0.5);
            double beta = v2.Max<double>();
            Matrix B = A.Copy();
            Matrix D = Matrix.Diagonal(B);
            Matrix E = Matrix.IdentityMatrix(B.Nrow);
            Matrix P = Matrix.IdentityMatrix(B.Nrow);
            Matrix L = Matrix.IdentityMatrix(B.Nrow, B.Ncol);
            Matrix Lt = Matrix.IdentityMatrix(B.Nrow, B.Ncol);
            Matrix C = D.Copy();
            int n = B.Nrow;
            for (int j = 0; j < n; j++)
            {
                int q = j;
                if (ApplyPivot)
                {

                    for (int i = j; i < n; i++)
                    {
                        double maxval = Abs(C[j, j]);
                        if (maxval < Abs(C[i, i]))
                        {
                            maxval = Abs(C[i, i]);
                            q = i;
                        }
                    }
                    if (q != j)
                    {
                        B.swap(j, q, SwapState.BothRowAndColumnWise);
                        P.swap(j, q, SwapState.RowWise);
                    }
                }
            }
            for (int j = 0; j < n; j++)
            {
                #region
                ////int q = j;
                ////if (ApplyPivot)
                ////{

                ////    for (int i = j; i < n; i++)
                ////    {
                ////        double maxval = Abs(C[j, j]);
                ////        if (maxval < Abs(C[i, i]))
                ////        {
                ////            maxval = Abs(C[i, i]);
                ////            q = i;
                ////        }
                ////    }
                ////    if (q != j)
                ////    {
                ////        B.swap(j, q, SwapState.BothRowAndColumnWise);
                ////        P.swap(j, q, SwapState.RowWise);
                ////    }
                ////}
                #endregion    
                for (int s = 0; s < j; s++)
                {
                    L[j, s] = C[j, s] / D[s, s];
                }
                for (int i = j; i < n; i++)
                {
                    double sum1 = 0.0;
                    for (int s = 0; s <= j-1; s++)
                    {
                        sum1 += L[j, s] * C[i, s];
                    }
                    C[i, j] = B[i, j] - sum1;
                }
                double thetaj = 0.0;
                if (j < n-1)
                {
                    Vector Vj = C[j, j + 1, 0, "col"];
                    thetaj = Vj.Max((k) => Abs(k));
                }
                D[j, j] = new Vector(new double[] { Abs(C[j, j]), Pow((thetaj / beta), 2), delta }).Max();
                E[j, j] = D[j, j] - C[j, j];
                if (j < n)
                {
                    for (int i = j; i < n; i++)
                    {
                        C[i, i] = C[i, i] - Pow(C[i, j], 2) / D[j, j];
                    }
                }
            }
            return (L, D, E,P);
        }
        #endregion

        #region QR Factorization
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        private static (Matrix Q, Matrix R) HouseHolderQR(Matrix A)
        {
            Matrix Qmat = new Matrix(A.Nrow, A.Ncol);
            Matrix Rmat = new Matrix(A.Nrow, A.Ncol);
            var HouseHolderFactorization = Matrix.HouseHolder(A);
            var Hlist = HouseHolderFactorization.ListofHouseHolderMatrix;
            if (Hlist.Count == 1)
            {
                Qmat = Hlist[0];
            }
            else
            {
                Qmat = Hlist[0];
                for (int i = 1; i < Hlist.Count; i++)
                {
                    Qmat = Qmat * Hlist[i];
                }
            }

            Rmat = HouseHolderFactorization.UpperTriangle;
            for (int i = 0; i < Rmat.Nrow; i++)
            {
                for (int j = 0; j < Rmat.Ncol; j++)
                {
                    if (j < i) { Rmat[i, j] = 0.0; }
                }
            }
            return (Qmat, Rmat);
        }
        private static (Matrix Qmatric, Matrix Rmatrix) GivensQR(Matrix A)
        {
            #region Givens Implementation
            Matrix B = A.Copy();
            Matrix copy = A.Copy();
            Matrix Qmat = Matrix.IdentityMatrix(A.Nrow);
            Matrix GivenRotMat = Matrix.IdentityMatrix(B.Nrow);
            for (int i = 0; i < A.Ncol - 1; i++)
            {
                Matrix fff = Matrix.SubMatrix(B, i);
                int p = fff.Nrow;
                for (int k = 1; k < p; k++)
                {
                    B = Matrix.SubMatrix(B, i);
                    GivenRotMat = Matrix.GivensMatrix(B, 0, k);
                    GivenRotMat = GivenModdified(GivenRotMat, i, k + i);
                    B = GivenRotMat * copy;
                    copy = B.Copy();
                    Qmat = Qmat * Matrix.Transpose(GivenRotMat);
                }

            }
            Matrix GivenModdified(Matrix GG, int rowindex, int colindex)
            {
                Matrix C = Matrix.IdentityMatrix(A.Nrow);
                C[rowindex, rowindex] = GG[0, 0];
                C[rowindex, colindex] = GG[0, 1];
                C[colindex, rowindex] = GG[1, 0];
                C[colindex, colindex] = GG[1, 1];
                return C;
            }
            #endregion

            return (Qmat, copy);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="Qrscheme"></param>
        /// <returns></returns>
        public static (Matrix Q, Matrix R) QR(Matrix A, QRFactorizationMethod Qrscheme)
        {
            Matrix Qmat = new Matrix(A.Nrow, A.Ncol);
            Matrix Rmat = new Matrix(A.Nrow, A.Ncol);
            switch (Qrscheme)
            {
                case QRFactorizationMethod.MGS:
                    Qmat = Matrix.ModifiedGramSchmidtOthogonalization(A);
                    Rmat = Matrix.Transpose(Qmat)* A;
                    break;
                case QRFactorizationMethod.HouseHolder:
                 var HH=   HouseHolderQR(A);
                    Qmat = HH.Q;
                    Rmat = HH.R;
                    break;
                case QRFactorizationMethod.Givens:
                    var qr = GivensQR(A);
                    Qmat = qr.Qmatric;
                    Rmat = qr.Rmatrix;
                    break;
                default:
                    break;
            }
            for (int i = 0; i < Rmat.Nrow; i++)
            {
                for (int j = 0; j < Rmat.Ncol; j++)
                {
                    if (j < i) { Rmat[i, j] = 0.0; }
                }
            }
            return (Qmat, Rmat);
        }
        #endregion

        #region Singular Value Decomposition
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="qrSch"></param>
        /// <returns></returns>
        public static (Matrix LeftSingularVec, Matrix SingularValues, Matrix RightSingularVec) SVD(Matrix A, QRFactorizationMethod qrSch = QRFactorizationMethod.HouseHolder)
        {
            int row = A.Nrow;
            int col = A.Ncol;
            bool containsZeroSingularValue = false;
            Matrix B = null!;
            Matrix V = null!;
            Matrix U = null!;
            Matrix S = null!;
            if (row > col)
            {
                B =A* Matrix.Transpose(A);
                U = Matrix.IdentityMatrix(B.Nrow, B.Ncol);
                V = new Matrix(A.Ncol,A.Ncol);
                S = new Matrix(row, col);
            }
            else
            {
                B = Matrix.Transpose(A) * A;
                V = new Matrix(B.Nrow, B.Ncol);
                U = new Matrix(A.Nrow, A.Nrow);
                S = new Matrix(col, row);
            }
            var eigSys = EigenSystem.EigenPair(B, qrSch);
             Vector lamdas = Matrix.GetDiagonalElements(eigSys.eigenValus);
            for (int i = 0; i < lamdas.ColVector.Length; i++)
            {
                if (Abs( lamdas[i]) <= 1e-10) { lamdas[i] = Abs(lamdas[i]); }
            }
            Vector LamdasSquaredRoot = lamdas ^ 0.5;
            if (row > col)
            {
                for (int i = 0; i < S.Ncol; i++)
                {
                    for (int j = 0; j < S.Nrow; j++)
                    {
                        if (i == j)
                        {
                            S[i, j] = LamdasSquaredRoot[i];
                            break;
                        }
                    }
                }
                U = eigSys.eigenVectors;
                U = Matrix.ModifiedGramSchmidtOthogonalization(U);
                for (int i = 0; i < row; i++)
                {
                    for (int j = 0; j < row; j++)
                    {
                        if (j >= col) { U[i, j] = 0.0; }
                    }
                }
               
                for (int i = 0; i < S.Ncol; i++)
                {
                    if (lamdas[i] != 0.0)
                    {
                        Vector vi = Matrix.Transpose(A) * (U[i, "col"] / lamdas[i]);
                        V[i, "col"] = vi;
                    }
                    else
                    {
                        V[i, "col"] = Vector.Ones(V.Nrow);
                        containsZeroSingularValue = true;
                    }
                }
                V = Matrix.ModifiedGramSchmidtOthogonalization(V);
                if (containsZeroSingularValue)
                {
                    V = Matrix.ModifiedGramSchmidtOthogonalization(V);
                }
            }
            else
            {
                for (int i = 0; i < S.Nrow; i++)
                {
                    for (int j = 0; j < S.Ncol; j++)
                    {
                        if (i == j)
                        {
                            S[i, j] = LamdasSquaredRoot[i];
                            break;
                        }
                    }
                }
                V = eigSys.eigenVectors;
                V = Matrix.ModifiedGramSchmidtOthogonalization(V);
                for (int i = 0; i < S.Ncol; i++)
                {
                    if (S[i, i] != 0.0)
                    {
                        Vector ui = A *( V[i, "col"] / S[i,i]);
                        U[i, "col"] = ui;
                    }
                    else
                    {
                        U[i, "col"] = Vector.Ones(U.Nrow);
                        containsZeroSingularValue = true;
                    }
                }
                U = Matrix.ModifiedGramSchmidtOthogonalization(U);
                if (containsZeroSingularValue)
                {
                    U = Matrix.ModifiedGramSchmidtOthogonalization(U);
                }
            }
           
            return (U, S, V);
        }
        #endregion
        #region Spectral Decomposition
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="qrscheme"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static (Matrix Q, Matrix D)SpectralDecomposition(Matrix A,QRFactorizationMethod qrscheme)
        {
            if (!Matrix.IsSymmetric(A)) { throw new Exception("Matrix is not diagonalizable"); }
            var eigsys = EigenSystem.EigenPair(A, qrscheme);
            Matrix lamda = eigsys.eigenValus;
            Matrix EigVectors = eigsys.eigenVectors;
            Matrix q = EigVectors.Normalized();
            return (q, lamda);
        }
        #endregion
        /// <summary>
        /// 
        /// </summary>
        public static int NumberOfRowInterChangeDuringLU = 0;
    }
}




