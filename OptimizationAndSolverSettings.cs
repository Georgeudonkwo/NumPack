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
    public enum UpdateStrategy {
        /// <summary>
        /// 
        /// </summary>
        LineSearch,
        /// <summary>
        /// 
        /// </summary>
        TrustRegion}
    /// <summary>
    /// 
    /// </summary>
  public  class OptimizationAndSolverSettings
    {
        /// <summary>
        /// 
        /// </summary>
        public OptimizationAndSolverSettings()
        {
            tol = Math.Pow(Eps, 1.0 / 3.0);
        }
        private  (double C1, double C2,bool IncludeCurvatureCondition) wolfeConstants ;
        /// <summary>
        /// 
        /// </summary>
        public double Alpha { get;internal set; } = 1.0;
        private double tol;
        /// <summary>
        /// 
        /// </summary>
        public int MaxIteration { get; set; } = 100;
        /// <summary>
        /// 
        /// </summary>
        public double Delta { get; set; } = 1e-8;
        /// <summary>
        /// 
        /// </summary>
        public double rho { get; set; } = 0.1;
        /// <summary>
        /// 
        /// </summary>
        public double tau { get; set; } = 0.001;
        /// <summary>
        /// 
        /// </summary>
        public Func<Vector, Matrix> CustomJac { get; set; } = null!;
        /// <summary>
        /// 
        /// </summary>
        public double Eps
        {
            get
            {
                double machineEps = 1.0;
                bool exit = true;
                do
                {
                    machineEps = machineEps / 2.0;
                    exit = (1.0 - machineEps) != 1.0;
                } while (exit);
                return 2 * machineEps;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        public double TrustRegionRadius { get; set; } = -1;
        /// <summary>
        /// 
        /// </summary>
        public UpdateStrategy UpdateMode { get; set; } = UpdateStrategy.LineSearch;
        /// <summary>
        /// 
        /// </summary>
        public ISpdModification MatrixModificationStrategy { get; set; } = new ModifiedCholeskyGMWalg();
        /// <summary>
        /// 
        /// </summary>
        public (double C1, double C2, bool IncludeCurvatureCondition) WolfeConstants
        {
            get { return wolfeConstants; }
            set
            {
                if (value.C1 > 0.25) { value.C1 = 0.1; }
                wolfeConstants = value;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        public double Tol { get => tol; set => tol = value; }
        internal double MaxStep { get; set; }
    }
}
