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
    public class AddingMultipleOfIdentityMatrix : ISpdModification
    {
        
        double Beta = 0.001;
        double tau = 0.0;
        int count = 0;
        double multiplier = 2.0;
        Matrix nspdMod = null!;
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Nspd"></param>
        /// <returns></returns>
        public Matrix SPD(Matrix Nspd)
        {
            Matrix eye = Matrix.IdentityMatrix(Nspd.Nrow);
                Beta =Sqrt(Matrix.FrobeniusNorm(Nspd));
                double MinDiad = Nspd.GetDiagonalElements().Min<double>();
                //if (MinDiad <= 0) { tau = -MinDiad + Beta; }
            if (MinDiad <= 0) { tau = Beta / 2; }// -MinDiad + Beta; }
            nspdMod = Nspd + tau * Matrix.IdentityMatrix(Nspd.Nrow);
            var choleskyDecomp = Cholesky(nspdMod);
            while(!choleskyDecomp.Success)
            {
                if (count % 2 == 0) { multiplier += 2; }
                //tau = multiplier * tau >= Beta ? multiplier * tau : Beta;
                tau = multiplier * tau >= Beta ? multiplier * tau : Beta/multiplier;
                nspdMod = Nspd + tau * Matrix.IdentityMatrix(Nspd.Nrow);
                choleskyDecomp = Cholesky(nspdMod);
                if (count > 10) break;
                count++;
            } 
                
            return nspdMod;
        }
        private static (Matrix L, Matrix LT,bool Success) Cholesky(Matrix A)
        {
            Matrix B = A.Copy();
            bool succeed = true;
            if (!Matrix.IsSymmetric(B)) { throw new Exception("Matrix must be symmetric"); }
            Matrix L = new Matrix(B.Nrow, B.Ncol);
            try
            {
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
            }
            catch (Exception)
            {
                succeed = false;
                return (null!, null!, succeed);
            }
            Matrix Lt = Matrix.Transpose(L);// Transpose(L);
            return (L,Lt,succeed);
        }
    }
}
