using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumSharp
{
    /// <summary>
    /// Uses the Gill Murray Wright algorithm
    /// </summary>
    public class ModifiedCholeskyGMWalg : ISpdModification
    {
        bool applyPivoting;
        /// <summary>
        /// 
        /// </summary>
        /// <param name="applyPivoting"></param>
        public ModifiedCholeskyGMWalg(bool applyPivoting=false)
        {
            this.applyPivoting = applyPivoting;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Nspd"></param>
        /// <returns></returns>
        public Matrix SPD(Matrix Nspd)
        {
            var res = MatrixFactorization.ModifiedCholeskyDecomp(Nspd, applyPivoting);
            Matrix L = res.L;
            Matrix D = res.D;
            Matrix SqrD = Matrix.TwPow(D, 0.5);
            Matrix M = L * SqrD;
            Matrix P = res.P;
            M = P * M;
            Matrix MMt = M * ~M;
            return MMt;
        }
    }
}
