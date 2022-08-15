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
   public interface ISpdModification
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Nspd"></param>
        /// <returns></returns>
        Matrix SPD(Matrix Nspd);
    }
}
