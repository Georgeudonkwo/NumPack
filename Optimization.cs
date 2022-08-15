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
   public class Optimization
    {
        #region Unconstrained

        private static Vector NewtonDirectionSearch(Vector g, Matrix H)
        {
            Vector gradient = g.Copy();
            Vector DirSearch = LinearSolvers.Solve(H, -gradient);
            return DirSearch;
        }
        private static Vector SteepestDecent(Vector g)
        {
            Vector DirSearch = g.Copy();
            return -DirSearch/ DirSearch.Norm();
        }
        private static Vector CGDirection(List<Vector> g,Vector S)
        {
            Vector DirSearch = g[1].Copy();
            double B = g[1].Dot((g[1] - g[0])) / (g[0].Dot(g[0]));
            if (Sign(B) == -1) { B = 0.0; }
            Vector res = -DirSearch + B * S;
            return res;
        }
        private static (Vector CurruntSolution, Vector DirSearch) PRCG(Func<Vector, double> fun, Vector x0, OptimizationAndSolverSettings setting = null!)
        {
            Vector x = x0.Copy();
            Vector g = BMath.Grad(fun, x, setting.Delta);
            Vector Dsir = -g;
            List<Vector> glist = new List<Vector>();
            glist.Add(g);
            int iter = 0;
            while (g.Norm() > setting.Tol)
            {
                iter++;
                double alp = BMath.InterpolationLineSearch(fun, x, Dsir, setting);
                x = x + alp * Dsir;
                Vector gnext = BMath.Grad(fun, x, setting.Delta);
                glist.Add(gnext);
                Dsir = CGDirection(glist, Dsir);
                glist.Remove(glist.First());
                if (iter > 100) break;
            }
            return (x, Dsir);
        }
        private static Matrix BFGS(Matrix H, Vector y, Vector s)
        {
            Matrix T2 = Matrix.OuterProduct(y, ~y) / Vector.Dot(~y, s);
            Matrix T3numerator = Matrix.OuterProduct(H * s, ~s) * H;
            double T3denom = Vector.Dot(~s, H * s);
            Matrix T3 = T3numerator / T3denom;
            Matrix bfgsUpdate = H + T2 - T3;
            return bfgsUpdate;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector Newton(Func<Vector, double> fun, Vector x,OptimizationAndSolverSettings setting=null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector xc = x.Copy();
            Vector xnext = x.Copy();
            setting.MaxStep = 1000 * xc.Norm();
            int counter = 0;
            Vector g = BMath. Grad(fun, xc);
            Matrix Hess =BMath.Hessian(fun, xc);
            bool exitCriterian = false;
            double fval = 0.0;
            double norm = 0.0;
            double alpha = setting.Alpha;
            double Trr = setting.TrustRegionRadius;
            if (Trr <= 0)
            {
                Vector ScaledG = g / g.Norm();
                Trr = ScaledG.Norm();
            }
            do
            {
                counter++;
                Hess = setting.MatrixModificationStrategy.SPD(Hess);
                Vector SearchDir = NewtonDirectionSearch(g, Hess);
                bool isDecentDir = Vector.Dot(~g, SearchDir) < 0.0;
                if (!isDecentDir)
                {
                    SearchDir = SteepestDecent(g);
                }
                switch (setting.UpdateMode)
                {
                    case UpdateStrategy.LineSearch:
                        if (!(BMath.WolfeCondition(fun, xc, SearchDir, alpha, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition)))
                        {
                            alpha = BMath.InterpolationLineSearch(fun, xc, SearchDir, setting);
                        }
                        xnext = xc + alpha * SearchDir;
                        break;
                    case UpdateStrategy.TrustRegion:
                        xnext = BMath.TrustRegionDriver(fun, SearchDir, g, xc, Hess, ref Trr, setting);
                        break;
                    default:
                        break;
                }
                Vector xdiff = xnext - xc;
                xc = xnext.Copy();
                g =BMath.Grad(fun, xc);
                Hess =BMath. Hessian(fun, xc);
                norm = g.Norm();
                exitCriterian = counter < setting.MaxIteration && norm >= 1e-6;
            } while (exitCriterian);
            fval = fun(xc);
            return xc;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector QuasiNewton(Func<Vector, double> fun, Vector x, OptimizationAndSolverSettings setting = null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector xc = x.Copy();
            Vector xnext = x.Copy();
            int counter = 0;
            Vector g = BMath.Grad(fun, xc);
            Matrix Hess = BMath.Hessian(fun, xc);
            bool exitCriterian = false;
            double fval = 0.0;
            double norm = 0.0;
            double alpha = setting.Alpha;
            double Trr = setting.TrustRegionRadius;
            if (Trr <= 0)
            {
                Vector ScaledG = g / g.Norm();
                Trr = ScaledG.Norm();
            }
            List<Vector> gList = new List<Vector>();
            
            do
            {
                counter++;
                Hess = setting.MatrixModificationStrategy.SPD(Hess);
                Vector SearchDir = NewtonDirectionSearch(g, Hess);
                bool isDecentDir = Vector.Dot(~g, SearchDir) < 0.0;
               
                if (!isDecentDir)
                {
                    gList.Add(g);
                    if (gList.Count > 2)
                    {
                        gList.Remove(gList.First());
                    }
                    if (!isDecentDir && counter > 1)
                    {
                        SearchDir = CGDirection(gList, SearchDir);
                    }
                    else
                    {
                        SearchDir = SteepestDecent(g);
                    }
                }
                switch (setting.UpdateMode)
                {
                    case UpdateStrategy.LineSearch:
                        if (!(BMath.WolfeCondition(fun, xc, SearchDir, alpha, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition)))
                        {
                            alpha = BMath.InterpolationLineSearch(fun, xc, SearchDir, setting);
                        }
                        xnext = xc + alpha * SearchDir;
                        break;
                    case UpdateStrategy.TrustRegion:
                        xnext = BMath.TrustRegionDriver(fun, SearchDir, g, xc, Hess, ref Trr, setting);
                        break;
                    default:
                        break;
                }
                Vector xdiff = xnext - xc;
                xc = xnext.Copy();
                Vector gnext = BMath.Grad(fun, xc);
                Vector gdiff = gnext - g;
                g = gnext.Copy();
                if (counter % 5 == 0)
                {
                    Hess = BMath.Hessian(fun, xc);
                }
                else
                {
                    Hess =BFGS(Hess,gdiff,alpha* SearchDir);
                }
               
                norm = g.Norm();
                exitCriterian = counter < setting.MaxIteration && norm >= setting.Tol;
            } while (exitCriterian);
            fval = fun(xc);
            return xc;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector CGNewton(Func<Vector, double> fun, Vector x, OptimizationAndSolverSettings setting = null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector xx = x.Copy();
            int counter = 0;
            Vector g = BMath.Grad(fun, xx);
            Matrix Hess = BMath.Hessian(fun, xx);
            bool exitCriterian = false;
            double fval = 0.0;
            double norm = 0.0;
            double alpha = 1.0;
            List<Vector> gList = new List<Vector>();

            do
            {
                counter++;
                Vector SearchDir = NewtonDirectionSearch(g, Hess);
                bool isDecentDir = Vector.Dot(~g, SearchDir) < 0.0;
                
                if (!isDecentDir )
                {
                    Vector iniGuess = Vector.Zeros(SearchDir.ColVector.GetLength(0));
                    SearchDir = PRCG(fun, iniGuess, setting).DirSearch;
                }
               
                if (!(BMath.WolfeCondition(fun, xx, SearchDir, alpha, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition)))
                {
                    alpha = BMath.InterpolationLineSearch(fun, xx, SearchDir, setting);
                }
                xx = xx + alpha * SearchDir;
                Vector gnext = BMath.Grad(fun, xx);
                Vector y = gnext - g;
                g = gnext.Copy();
                if (counter % 5 == 0)
                {
                    Hess = BMath.Hessian(fun, xx);
                }
                else
                {
                    Hess = BFGS(Hess, y, alpha * SearchDir);
                }

                norm = g.Norm();
                exitCriterian = counter < setting.MaxIteration && norm >= 1e-6;
            } while (exitCriterian);
            fval = fun(xx);
            return xx;
        }

        #endregion
        #region constrained


        #endregion
    }
}
