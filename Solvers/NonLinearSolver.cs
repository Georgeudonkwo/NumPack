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
 public  class NonLinearSolver
    {
        #region Multi variable 
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="InitialGuess"></param>
        /// <param name="tol"></param>
        /// <param name="MaxIter"></param>
        /// <param name="CusJac"></param>
        /// <returns></returns>
        public static Vector Newton(Func<Vector, Vector> Fun, Vector InitialGuess,double tol= 1e-6,int MaxIter=100, Func<Vector, Matrix> CusJac = null!)
        {
            Vector x = InitialGuess.Copy();
            Matrix jac = null!;
            bool exitCriteria = false;
            double functionNorm = 0.0;
            Vector fx = Fun(x);
            double initialFunctionNorm = fx.Norm();
            int iter = 0;
            do
            {
                iter++;
                if (CusJac != null) { jac = CusJac(x); }
                else { jac = BMath.Jacobian(Fun, x); }
                Vector searchDir = -LinearSolvers.Solve(jac, fx);
                x = x + searchDir;
                 fx = Fun(x);
                functionNorm = fx.Norm();
                bool FxConverged = functionNorm > tol * (initialFunctionNorm + 1.0);
                bool stepConverged = (searchDir / x).Norm() > tol;
                exitCriteria = stepConverged || FxConverged;
                bool stepNoteffective = Abs(functionNorm - initialFunctionNorm) < 1e-16;
                if (stepNoteffective) { break; }
                else { initialFunctionNorm = functionNorm; }
                if (iter > MaxIter) break;
            } while (exitCriteria);
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="InitialGuess"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector ModifiedNewton(Func<Vector, Vector> Fun, Vector InitialGuess, OptimizationAndSolverSettings setting = null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector xc = InitialGuess.Copy();
            Vector xnext = InitialGuess.Copy();
            Matrix jac = null!;
            bool exitCriteria = false;
            double functionNorm = 0.0;
            Vector fx = Fun(xc);
            double initialFunctionNorm = fx.Norm();
            double alpha = setting.Alpha;
            double xdiffNorm = 0.0;
            int iter = 0;
            Func<Vector, double> objFun = new Func<Vector, double>(c => 0.5 * (Fun(c) ^ 2).Sum());
            Vector g = BMath.Grad(objFun, xc, setting.Delta);
            double Trr = setting.TrustRegionRadius;
            if (Trr <= 0)
            {
                Vector ScaledG = g / g.Norm();
                Trr = ScaledG.Norm();
            }
            Vector Sdir = null!;
            do
            {
                iter++;
                if (setting.CustomJac != null) { jac = setting.CustomJac(xc); }
                else { jac = BMath.Jacobian(Fun, xc, setting.Delta); }
                double condEst = jac.Condest();
                #region
                //if (condEst > 1e8)
                //{
                //    Matrix jacMod = setting.MatrixModificationStrategy.SPD(jac);
                //    jacMod = ~jacMod * jacMod;
                //    Vector Jr = ~jacMod * fx;
                //    Sdir = -LinearSolvers.Solve(jacMod, Jr);
                //    #region
                //    //double norm1 = MATRIX.MatrixNorm_1((~jac * jac));
                //    //double multiplier = Sqrt(setting.Eps) * norm1;
                //    //Matrix eye = multiplier * MATRIX.IdentityMatrix(jac.Mat.GetLength(0));
                //    //Matrix jacMod = ~jac * jac + eye;
                //    //Vector Jr = ~jacMod * fx;
                //    #endregion
                //}
                //else
                //{
                //Sdir = -LinearSolvers.Solve(jac, fx);
                //}
                //double cosTheta = Vector.Dot(~-Sdir, g) / (Sdir.Norm() * g.Norm());
                //if (cosTheta <= setting.rho && condEst > 1e8)
                #endregion
                Sdir = -LinearSolvers.Solve(jac, fx);
                bool isDecentDir = Vector.Dot(~g, Sdir) < 0.0;
                if (isDecentDir && condEst > 1e8)
                {
                    Matrix jacMod = setting.MatrixModificationStrategy.SPD(jac);
                    jacMod = ~jacMod * jacMod;
                    Vector Jr = ~jacMod * fx;
                    Sdir = -LinearSolvers.Solve(jacMod, Jr);
                    #region
                    //bool exit = true;
                    //int count = 0;
                    //double tau = 1.0;
                    //do
                    //{
                    //    count++;
                    //    double norm1 = MATRIX.MatrixNorm_1((~jac * jac));
                    //    tau = Sqrt(setting.Eps) * norm1;
                    //    Matrix eyes = tau * MATRIX.IdentityMatrix(jac.Mat.GetLength(0));
                    //    Matrix jacModd = (~jac * jac) + eyes;
                    //    Vector jrr = ~jacModd * fx;
                    //    Sdir = -LinearSolvers.Solve(jacModd, jrr);
                    //    cosTheta = Vector.Dot((-Sdir).GetRowVector(), g) / (Sdir.Norm() * g.Norm());
                    //    exit = cosTheta < setting.rho;
                    //    tau = tau + norm1;
                    //    if (count > 10) break;
                    //} while (exit);
                    #endregion
                }
                switch (setting.UpdateMode)
                {
                    case UpdateStrategy.LineSearch:
                        if (!(BMath.WolfeCondition(objFun, xc, Sdir, alpha, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition)))
                        {
                            alpha = BMath.InterpolationLineSearch(objFun, xc, Sdir, setting);
                        }
                        xnext = xc + alpha * Sdir;
                        break;
                    case UpdateStrategy.TrustRegion:
                        xnext = BMath.TrustRegionDriver(objFun, Sdir, g, xc, jac, ref Trr,setting);
                        break;
                    default:
                        break;
                }
                Vector xdiff = xnext - xc;
                xdiffNorm = xdiff.Norm();
                xc = xnext.Copy();
                fx = Fun(xc);
                functionNorm = fx.Norm();
                g = BMath.Grad(objFun, xc, setting.Delta);
                bool FxConverged = functionNorm > setting.Tol * (initialFunctionNorm + 1.0);
                bool stepConverged = (Sdir / xc).Norm() > setting.Tol;
                exitCriteria = stepConverged || FxConverged;
                bool stepNoteffective = Abs(functionNorm - initialFunctionNorm) < setting.Eps;
                if (stepNoteffective) { break; }
                else { initialFunctionNorm = functionNorm; }
                if (iter > setting.MaxIteration) break;

            } while (exitCriteria);
            return xc;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="InitialGuess"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector QuasiNewton(Func<Vector, Vector> Fun, Vector InitialGuess, OptimizationAndSolverSettings setting = null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector xc = InitialGuess.Copy();
            Vector xnext = InitialGuess.Copy();
            Matrix jac = null!;
            bool exitCriteria = false;
            double functionNorm = 0.0;
            Vector fx = Fun(xc);
            double initialFunctionNorm = fx.Norm();
            double alpha = setting.Alpha;
            double xdiffNorm = 0.0;
            int iter = 0;
           
            if (setting.CustomJac != null) { jac = setting.CustomJac(xc); }
            else { jac = BMath.Jacobian(Fun, xc, setting.Delta); }
            Func<Vector, double> objFun = new Func<Vector, double>(c => 0.5 * (Fun(c) ^ 2).Sum());
            Vector g = BMath.Grad(objFun, xc, setting.Delta);
            double Trr = setting.TrustRegionRadius;
            if (Trr <= 0)
            {
                Vector ScaledG = g / g.Norm();
                Trr = ScaledG.Norm();
            }
            Vector Sdir = null!;
            do
            {
                iter++;
                double condEst = jac.Condest();
                #region
                //if (condEst > 1e8)
                //{
                //    #region
                //    //double norm1 = MATRIX.MatrixNorm_1((~jac * jac));
                //    //double multiplier = Sqrt(setting.Eps) * norm1;
                //    //MATRIX eye = multiplier * MATRIX.IdentityMatrix(jac.Mat.GetLength(0));
                //    //MATRIX jacMod = ~jac * jac + eye;
                //    //Vector Jr = ~jacMod * fx;
                //    //Sdir = -LinearSolvers.Solve(jacMod, Jr);
                //    #endregion
                //    MATRIX jacMod = setting.MatrixModificationStrategy.SPD(jac);
                //    jacMod = ~jacMod * jacMod;
                //    Vector Jr = ~jacMod * fx;
                //    Sdir = -LinearSolvers.Solve(jacMod, Jr);
                //}
                //else
                //{
                //    Sdir = -LinearSolvers.Solve(jac, fx);
                //}
                //double cosTheta = Vector.Dot(~-g, Sdir) / (Sdir.Norm() * g.Norm());
                // if (cosTheta <= setting.rho && condEst > 1e8)
                #endregion
                Sdir = -LinearSolvers.Solve(jac, fx);
                bool isDecentDir = Vector.Dot(~g, Sdir) < 0.0;
                if (isDecentDir&& condEst > 1e8)
                {
                    Matrix jacMod = setting.MatrixModificationStrategy.SPD(jac);
                    jacMod = ~jacMod * jacMod;
                    Vector Jr = ~jacMod * fx;
                    Sdir = -LinearSolvers.Solve(jacMod, Jr);
                    #region
                    //bool exit = true;
                    //int count = 0;
                    //double tau = 1.0;
                    //do
                    //{
                    //    count++;
                    //    double norm1 = MATRIX.MatrixNorm_1((~jac * jac));
                    //    tau = Sqrt(setting.Eps) * norm1;
                    //    MATRIX eyes = tau * MATRIX.IdentityMatrix(jac.Mat.GetLength(0));
                    //    MATRIX jacModd = (~jac * jac) + eyes;
                    //    Vector jrr = ~jacModd * fx;
                    //    Sdir = -LinearSolvers.Solve(jacModd, jrr);
                    //    cosTheta = Vector.Dot((-Sdir).GetRowVector(), g) / (Sdir.Norm() * g.Norm());
                    //    exit = cosTheta < setting.rho;
                    //    tau = tau + norm1;
                    //    if (count > 10) break;
                    //} while (exit);
                    #endregion
                }
                switch (setting.UpdateMode)
                {
                    case UpdateStrategy.LineSearch:
                        if (!(BMath.WolfeCondition(objFun, xc, Sdir, alpha, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition)))
                        {
                            alpha = BMath.InterpolationLineSearch(objFun, xc, Sdir, setting);
                        }
                        xnext = xc + alpha * Sdir;
                        break;
                    case UpdateStrategy.TrustRegion:
                        xnext = BMath.TrustRegionDriver(objFun, Sdir, g, xc, jac, ref Trr, setting);
                        break;
                    default:
                        break;
                }
                Vector xdiff = xnext - xc;
                xdiffNorm = xdiff.Norm();
                xc = xnext.Copy();
                Vector fnextX = Fun(xc);
                Vector y = fnextX - fx;
                fx = fnextX.Copy();
                 g = BMath.Grad(objFun, xc, setting.Delta);
                if (iter % 5 == 0)
                {
                    if (setting.CustomJac != null) { jac = setting.CustomJac(xc); }
                    else { jac = BMath.Jacobian(Fun, xc, setting.Delta); }
                }
                else
                {
                    jac = BroydenUpdat(jac, y, Sdir);
                }
                functionNorm = fx.Norm();
                bool FxConverged = functionNorm > setting.Tol * (initialFunctionNorm + 1.0);
                bool stepConverged = (Sdir / xc).Norm() > setting.Tol;
                exitCriteria = stepConverged || FxConverged;
                bool stepNoteffective = Abs(functionNorm - initialFunctionNorm) < 1e-16;
                if (stepNoteffective) { break; }
                else { initialFunctionNorm = functionNorm; }
                if (iter > setting.MaxIteration) break;

            } while (exitCriteria);
            return xc;
        }
        private static Matrix BroydenUpdat(Matrix A,Vector y,Vector s)
        {
            Matrix AUpdate = A + Matrix.OuterProduct((y - A * s), ~s) / s.Dot(s);
            return AUpdate;
        }
        #endregion
        #region Univariable
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="initialGues"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static double NewtonFzero(Func<double,double>fun,double initialGues,OptimizationAndSolverSettings setting=null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            double x = initialGues;
            double f = fun(x);
            int iter = 0;
            while (Abs(f)>setting.Tol && iter<=setting.MaxIteration)
            {
                iter++;
                double dfdx = BMath.Diff(fun, x,setting.Delta);
                double dx = (-1) *( f / dfdx);
                x = x + dx;
                f = fun(x);
            }
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Fun"></param>
        /// <param name="Lb"></param>
        /// <param name="Ub"></param>
        /// <param name="Tol"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double BrentFzero(Func<double, double> Fun, double Lb = 0.0, double Ub = 1.0, double Tol = 1e-8)
        {
            try
            {
                double A, tmp, B, C, d, fa, fb;
                double fc, s, fs, tmp2;
                bool mflag;
                int i;
                if (Lb == default(double)) { Lb = 0.0; }
                if (Ub == default(double)) { Lb = 100.0; }
                double lowerLimit = Lb;
                double upperLimit = Ub;
                double errorTol = Tol;

                A = lowerLimit;
                B = upperLimit;
                C = 0.0;
                d = double.MaxValue;

                fa = Fun(A);
                fb = Fun(B);

                fc = 0;
                s = 0;
                fs = 0;


                if ((fa * fb) >= 0)
                {
                    if (fa < fb)
                    {
                        if (B <= Lb) { B = Lb; }
                        if (B >= Ub) { B = Ub; }
                        return B;
                    }
                    else
                    {
                        if (A <= Lb) { A = Lb; }
                        if (A >= Ub) { A = Ub; }
                        return A;
                    }
                }

                // if |f(a)| < |f(b)| then swap (a,b) end if
                if (Math.Abs(fa) < Math.Abs(fb))
                {

                    tmp = A;
                    A = B;
                    B = tmp;
                    tmp = fa;
                    fa = fb;
                    fb = tmp;
                }

                C = A;
                fc = fa;
                mflag = true;

                i = 0;


                while (!(fb == 0) && (Math.Abs(A - B) > errorTol))
                {

                    if ((fa != fc) && (fb != fc))
                    {

                        // Inverse quadratic interpolation

                        s = A * fb * fc / (fa - fb) / (fa - fc) + B * fa * fc / (fb - fa) / (fb - fc) + C * fa * fb / (fc - fa) / (fc - fb);
                    }

                    else

                    // Secant Rule
                    {
                        s = B - fb * (B - A) / (fb - fa);

                    }

                    tmp2 = (3 * A + B) / 4.0;


                    if (!(((s > tmp2) && (s < B)) || ((s < tmp2) && (s > B))) || ((mflag && (Math.Abs(s - B) >= (Math.Abs(B - C) / 2))) || (!mflag && (Math.Abs(s - B) >= (Math.Abs(C - d) / 2)))))
                    {

                        s = (A + B) / 2;

                        mflag = true;
                    }

                    else
                    {

                        if (mflag && (Math.Abs(B - C) < errorTol) || (!mflag && (Math.Abs(C - d) < errorTol)))
                        {

                            s = (A + B) / 2.0;

                            mflag = true;
                        }

                        else
                        {
                            mflag = false;

                        }

                    }

                    fs = Fun(s);

                    d = C;
                    C = B;
                    fc = fb;

                    if (fa * fs < 0)
                    {

                        B = s;
                        fb = fs;
                    }

                    else
                    {
                        A = s;
                        fa = fs;

                    }

                    // if |f(a)| < |f(b)| then swap (a,b) end if

                    if (Math.Abs(fa) < Math.Abs(fb))
                    {


                        tmp = A;
                        A = B;
                        B = tmp;
                        tmp = fa;
                        fa = fb;
                        fb = tmp;

                    }

                    i = i + 1;
                    if (i > 2000)
                    {
                        break;
                    }

                }
                if (B <= Lb) { B = Lb; }
                if (B >= Ub) { B = Ub; }
                return B;
            }
            catch (Exception e)
            {
                string message = e.Message + " " + "at" + " " + e.TargetSite!.Name + " " + e.Source;
                throw new Exception(message, new Exception(e.Message));
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Bisectfzero(Func<double, double> fun, double a, double b)
        {

            Func<double, double> Fun = new Func<double, double>(x => { return fun(x); });
            double c, fc = 0, fa, Tol = 1e-6, Sol; bool IsConverged = false; int iter = 0;
            while (!IsConverged && iter < 100)
            {
                c = (a + b) / 2; // new midpoint
                fc = Fun(c);
                IsConverged = (Math.Abs(fc) < Tol || (b - a) / 2 < Tol);
                if (IsConverged) //then # solution found
                {
                    Sol = c;
                    return Sol;
                }
                iter++; // increment step counter
                fa = Fun(a);
                if (fc * fa > 0)
                {
                    a = c; // new interval
                }
                else
                {
                    b = c; // new interval
                }
            }
            Sol = a;
            if (Sol <= a) { Sol = a; }
            if (Sol >= b) { Sol = b; }
            return Sol;
        }
        #endregion
    }
}
