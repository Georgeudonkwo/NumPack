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
    public enum TrustregionMethod {
        /// <summary>
        /// 
        /// </summary>
        Hook,
        /// <summary>
        /// 
        /// </summary>
        DogLeg
    }
    /// <summary>
    /// Basic utility mathematical routines
    /// </summary>
    public static class BMath
    {
        /// <summary>
        /// first order central difference routine for differentiating univariate function <paramref name="fun"/>
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="del"></param>
        /// <returns></returns>
        public static double Diff(Func<double, double> fun, double x, double del = 1e-6)
        {
            double delk = del * (1 + Abs(x));
            double next = fun((x + delk));
            double previuos = fun((x - delk));
            double der = (next - previuos) / (2 * delk);
            return der;
        }
        /// <summary>
        /// second order finite difference routine
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="del"></param>
        /// <returns></returns>
        public static double Diff2(Func<double, double> fun, double x, double del = 1e-6)
        {
           
            double delk = del * (1 + Abs(x));
            double f = fun(x);
            double fnext = fun((x + delk));
            double fpreviuos = fun((x - delk));
            double der = (fpreviuos - 2 * f + fnext) / Pow(delk, 2);
            return der;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="del"></param>
        /// <returns></returns>
        public static Matrix Jacobian(Func<Vector, Vector> fun, Vector x, double del = 1e-6)
        {
            Vector xx = x.Copy();
            Vector f1 = fun(xx);
            Matrix jac = new Matrix(f1.Length, x.ColVector.GetLength(0));
            for (int i = 0; i < x.ColVector.GetLength(0); i++)
            {
                double delj = del * (1.0 + Math.Abs(xx[i]));
                xx[i] = xx[i] - delj;
                Vector Dec = fun(xx);
                xx[i] = x[i];
                xx[i] = xx[i] + delj;
                Vector Inc = fun(xx);
                xx[i] = x[i];
                jac[i, "col"] = (Inc - Dec) / (2 * delj);
            }
            return jac;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="del"></param>
        /// <returns></returns>
        public static Vector Grad(Func<Vector, double> fun, Vector x, double del = 1e-6)
        {
            Vector xx = x.Copy();
            Vector g = new Vector(x.ColVector.Length);
            for (int i = 0; i < xx.ColVector.GetLength(0); i++)
            {
                double deli = del * (1.0 + Abs(xx[i]));
                xx[i] = xx[i] - deli;
                double Dec = fun(xx);
                xx[i] = x[i];
                xx[i] = xx[i] + deli;
                double Inc = fun(xx);
                xx[i] = x[i];
                g[i] = (Inc - Dec) / (2 * deli);
            }
            return g;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="x"></param>
        /// <param name="del"></param>
        /// <returns></returns>
        public static Matrix Hessian(Func<Vector, double> fun, Vector x, double del = 1e-5)
        {
            Matrix H = new Matrix(x.ColVector.GetLength(0), x.ColVector.GetLength(0));
            Vector xx = x.Copy();
            for (int i = 0; i < xx.Count(); i++)
            {
                double delj = del * (1.0 + Math.Abs(xx[i]));
                xx[i] = xx[i] - delj;
                Vector Dec = Grad(fun, xx);
                xx[i] = x[i];
                xx[i] = xx[i] + delj;
                Vector Inc = Grad(fun, xx);
                xx[i] = x[i];
                for (int j = 0; j < x.Count(); j++)
                {
                    H[i, j] = (Inc[j] - Dec[j]) / (2 * delj);
                }
            }
            return H;
        }
        internal static bool WolfeCondition(Func<Vector, double> fun, Vector x, Vector Sdir, double alp, double C1, double C2, bool rigorous)
        {
            Vector alpDir = alp * Sdir;
            Vector nextX = x + alpDir;
            double f_nextX = fun(nextX);
            double f_x = fun(x);
            Vector grad_x = Grad(fun, x);
            Vector grad_nextx = Grad(fun, nextX);
            double ScaledDirectionDerivative_x_01 = C1 * alp * Vector.Dot(grad_x.GetRowVector(), Sdir);
            double ScaledDirectionDerivative_x_02 = C2 * Vector.Dot(grad_x.GetRowVector(), Sdir);
            double ScaledDirectionDerivative_alp = Vector.Dot(grad_nextx.GetRowVector(), Sdir);
            bool sufficientDecrease = f_nextX <= (f_x + ScaledDirectionDerivative_x_01);
            bool CurvatureCondition = Abs(ScaledDirectionDerivative_alp) <= Abs(ScaledDirectionDerivative_x_02);
            bool wlofeCond = sufficientDecrease;
            if (rigorous) { wlofeCond = sufficientDecrease && CurvatureCondition; }
            return wlofeCond;
        }
        internal static double BackTrackingLineSearch(Func<Vector, double> fun, Vector x, Vector Sdir, double alp, double C1, double C2, bool rigorous)
        {
            int count = 0;
            while ((!WolfeCondition(fun, x, Sdir, alp, C1, C2, rigorous)))
            {
                count++;
                alp = 1.0 / Pow(2, count);
                if (count > 10) break;
            }
            double alpha = alp;
            return alpha;
        }
        private static double OneDimensionalMinDerivative(Func<Vector, double> fun, Vector X, Vector DirS, double alp)
        {
            double del = 1e-6;
            double delk = del * (1 + Abs(alp));
            double alpDre = alp - delk;
            double alpInc = alp + delk;
            Vector xInc = X + alpInc * DirS;
            Vector xDre = X + alpDre * DirS;
            double fForward = fun(xInc);
            double fBackward = fun(xDre);
            double derivative = (fForward - fBackward) / (2 * delk);
            return derivative;

        }
        private static double OneDimensionalMinFunc(Func<Vector, double> fun, Vector X, Vector DirS, double alp)
        {
            Vector nextX = X + alp * DirS;
            double FunVal = fun(nextX);
            return FunVal;
        }
        internal static double InterpolationLineSearch(Func<Vector, double> fun, Vector x, Vector Sdir, OptimizationAndSolverSettings sets)
        {
            double newtLen = Sdir.Norm();
            if (newtLen > sets.MaxStep )
            {
                Sdir = Sdir * (sets.MaxStep / newtLen);
                newtLen = sets.MaxStep;
            }
            Vector scaledSdir = Sdir / x.Norm_Infinity();
            double relLenght = scaledSdir.Norm_Infinity();
            double MinLamda = sets.Tol / relLenght;
            int count = 0;
            double alpFirst = sets.Alpha;
            double f0 = OneDimensionalMinFunc(fun, x, Sdir, 0.0);
            double f0Derivative = Grad(fun, x, sets.Delta).Dot(Sdir);
            // double foder = OneDimensionalMinDerivative(fun, x, Sdir, 0.0);
            double f_optimumAlpha = OneDimensionalMinFunc(fun, x, Sdir, sets.Alpha);

            bool sufficientDecrease = f_optimumAlpha <= (f0 + sets.WolfeConstants.C1 * sets.Alpha * f0Derivative);
            if (sufficientDecrease)
            {
                if (double.IsNaN(alpFirst) || double.IsInfinity(alpFirst)) { return MinLamda; }
                return alpFirst;
            }
            double alpnext = (-1) * ((f0Derivative * Pow(sets.Alpha, 2)) / (2 * (f_optimumAlpha - f0 - f0Derivative * sets.Alpha)));
            if (alpnext > 0.5) { alpnext = 0.5; }
            if (alpnext < 0.1) { alpnext = 0.1; }
            double funValAtAlpNext = OneDimensionalMinFunc(fun, x, Sdir, alpnext);

            sufficientDecrease = funValAtAlpNext <= (f0 + sets.WolfeConstants.C1 * alpnext * f0Derivative);
            if (sufficientDecrease)
            {
                //if (alpnext < 1e-5) { alpnext = 1e-5; };
                if (alpnext < 0.1* alpnext) { alpnext = 0.1 * alpnext; };
                if (alpnext > 0.5 * alpnext) { alpnext = 0.5 * alpnext; };
                if (double.IsNaN(alpnext) || double.IsInfinity(alpnext)) { return MinLamda; }
                return alpnext;
            }
            List<double> funvals = new List<double>() { f_optimumAlpha, funValAtAlpNext };
            List<double> alpList = new List<double>() { alpFirst, alpnext };
            double alp = 0.0;
            while ((!sufficientDecrease))
            {
                count++;

                double lamdaPrev = alpList[1];
                double lamda2Prev = alpList[0];
                double item001 = 1.0 / Pow(lamdaPrev, 2);
                double item011 = -1.0 / Pow(lamda2Prev, 2);
                double item101 = -lamda2Prev / Pow(lamdaPrev, 2);
                double item111 = lamdaPrev / Pow(lamda2Prev, 2);
                double scale1 = 1.0 / (lamdaPrev - lamda2Prev);
                double item00 = Pow(lamda2Prev, 2);
                double item01 = -Pow(lamdaPrev, 2);
                double item10 = -Pow(lamda2Prev, 3);
                double item11 = Pow(lamdaPrev, 3);
                double scale = 1.0 / (Pow(lamda2Prev, 2) * Pow(lamdaPrev, 2) * (lamdaPrev - lamda2Prev));
                Matrix mat = new Matrix(new double[2, 2] { { item00, item01 }, { item10, item11 } });

                double fPrev = OneDimensionalMinFunc(fun, x, Sdir, lamdaPrev);
                double f2Prev = OneDimensionalMinFunc(fun, x, Sdir, lamda2Prev);

                Vector vec = new Vector(2);
                vec[0] = fPrev - f0 - f0Derivative * lamdaPrev;
                vec[1] = f2Prev - f0 - f0Derivative * lamda2Prev;
                double factor = 1.0 / (Pow(alpList[0], 2) * Pow(alpList[1], 2) * (alpList[1] - alpList[0]));
                Vector ab = factor * (mat * vec);
                if (ab[0] == 0.0)
                {
                    double gtSdir = Grad(fun, x, sets.Delta).Dot(Sdir);
                    alp = gtSdir / (2 * ab[1]);
                }
                else
                {
                    alp = (-ab[1] + Sqrt((Pow(ab[1], 2) - 3 * ab[0] * f0Derivative))) / (3 * ab[0]);
                }
                if (alp < MinLamda) { alp = 1.0;break; }
                if (alp > 0.5 * lamdaPrev) { alp = 0.5 * lamdaPrev; }
                if (alp < 0.1 * lamdaPrev) { alp = 0.1 * lamdaPrev; }
                double funValAtAlp = OneDimensionalMinFunc(fun, x, Sdir, alp);
                sufficientDecrease = funValAtAlp <= (f0 + sets.WolfeConstants.C1 * alp * f0Derivative);
                alpList.RemoveAt(0);
                funvals.RemoveAt(0);
                alpList.Add(alp);
                funvals.Add(funValAtAlp);
                if (alp < MinLamda) { break; }
                if (count > 10) break;
            }
            if (double.IsNaN(alp) || double.IsInfinity(alp)) { return 0.001; }
            //if (alp < 0.1) { alp = 0.1; }
            //if (alp > 0.5) { alp = 0.5; }
            return alp;
        }
        internal static (double, Vector) mue(Matrix H, Vector g, ref double u, double lamda, double delPrev, ref bool newtonStepUsed, bool firstCalled)
        {
            Func<Matrix, Vector, Vector> QFunc = new Func<Matrix, Vector, Vector>((m, S) =>
                 {
                     double len = S.Norm();
                     double Phii = len - lamda;
                     var qr = MatrixFactorization.QR(m, QRFactorizationMethod.HouseHolder);
                     Vector PhiSolvee = (qr.R ^ -1) * ~qr.Q * S;
                     double PhiiDerivative = S.Dot(PhiSolvee) / len;
                     Vector res = new Vector(2);
                     res[0] = Phii;
                     res[1] = PhiiDerivative;
                     return res;
                 });
            Matrix eye = Matrix.IdentityMatrix(H.Nrow);
            #region Initialized
            double Uc = 0.0;
            double UL = 0.75 * lamda;
            double Uu = 1.5 * lamda;
            Vector Su = -LinearSolvers.Solve(H, g);
            double Newtonlen = Su.Norm();
            if (Newtonlen <= 1.5 * lamda)
            {
                newtonStepUsed = true;
                Uc = 0.0;
                lamda = lamda < Newtonlen ? lamda : Newtonlen;
                return (Uc, Su);
            }
            else
            {
                newtonStepUsed = false;
                var Qf = QFunc(H, Su);
                if (u > 0.0)
                {

                    Uc = u - ((Qf[0] + delPrev) / lamda) * (((delPrev - lamda) + Qf[0]) / Qf[1]);
                }
                if (firstCalled)
                {
                    firstCalled = false;
                }
                UL = -Qf[0] / Qf[1];
                Uu = g.Norm() / lamda;
            }
            #endregion
            bool exit = false;
            double Norm = Su.Norm();
            double delLowerBound = 0.75 * lamda;
            double delUpperBound = 1.5 * lamda;
            exit = Norm < UL && Norm > Uu;
            while (exit)
            {
                if (Uc < UL || Uc > Uu) { Uc = Sqrt(UL * Uu) > 1e-3 * Uu ? Sqrt(UL * Uu) : 1e-3 * Uu; }
                Matrix Ueye = Uc * eye;
                Matrix Hmod = H + Ueye;
                Su = -LinearSolvers.Solve(Hmod, g);
                Vector Qvec = QFunc(Hmod, Su);
                double stepLenght = Su.Norm();
                exit = (Norm >= UL && Norm <= Uu) || (Uu - UL <= 0.0);
                if (exit) { return (Uc, Su); }
                Uc = Uc - Qvec[0] / Qvec[1];
                UL = UL > Uc - (Qvec[0] / Qvec[1]) ? UL : Uc - (Qvec[0] / Qvec[1]);

                double UuCurrent = g.Norm() / lamda;
                if (UL < Uc) { UL = Uc; }
                if (Qvec[0] < 0.0)
                {
                    Uu = Uc;
                    Uc = Uc - (stepLenght / lamda) * (Qvec[0] / Qvec[1]);
                }
            }
            return (Uc, Su);
        }
        private static double Landa(Matrix Q, Matrix R, Vector S, double TrRadius, double LamdaPrev)
        {
            var Qres = CalculatePhiandPhiDer( Q, R, S, TrRadius);
            double term1 = S.Norm() / TrRadius;
            double term2 = Qres.Phi / Qres.PhiDerivative;
            double currentLamda = LamdaPrev - (term1 * term2);
            return currentLamda;
        }
        private static (double Phi, double PhiDerivative) CalculatePhiandPhiDer(Matrix Q, Matrix R, Vector S, double TrRadius)
        {
            double len = S.Norm();
            double Phii = len - TrRadius;
            //var qr = MatrixFactorization.QR(m, QRFactorizationMethod.HouseHolder);
            Vector Rhs= ~Q * S;
            Vector PhiSolvee = LinearSolvers.BackSubstition(R,Rhs);
            Vector wwwwww = (R ^ -1) * Rhs;
            double PhiiDerivative =(-1)* S.Dot(PhiSolvee) / len;
            return (Phii, PhiiDerivative);
        }
        internal static(Vector Sdir,double lan,double Tr) Hook(Matrix H,Vector g,Vector Sn, double TrRadius,ref double lamdaPrev,ref double Qprev,ref double QderPrev,bool firstHook,OptimizationAndSolverSettings settings)
        {
            Func<Matrix, double, Matrix> PDM = new Func<Matrix, double, Matrix>((M, L) =>
                 {
                     Matrix eye = Matrix.IdentityMatrix(M.Nrow);
                     Matrix Hmodified = M + L * eye;
                     return Hmodified;
                 });
            Func<Matrix, (Matrix Q, Matrix R)> QrFac = new Func<Matrix, (Matrix, Matrix)>(m => MatrixFactorization.QR(m, QRFactorizationMethod.HouseHolder));
            Vector S = Sn.Copy();
            double hi = 1.5;double low = 0.75;
            double lamda = lamdaPrev;
            double Trr = 0.0;
            double newtonLen = Sn.Norm();
            bool newtonStepTaken = newtonLen <= hi * TrRadius;
            if (newtonStepTaken)
            {
                lamda = 0.0;
                lamdaPrev = lamda;
                Trr = TrRadius<=newtonLen?TrRadius:newtonLen;
                return (S, lamda, Trr);
            }
            Matrix Hmod = PDM(H, lamda);
            var qr = QrFac(Hmod);
            if (lamda > 0.0)
            {
                lamda = Landa(qr.Q, qr.R, S, TrRadius, lamdaPrev);//initial 
            }
            double LamdaLow = 0.0;
            double lamdaUp = 0.0;
            if (firstHook)
            {
                Matrix HmodIni = PDM(H, 0.0);
                var qrIni = QrFac(HmodIni);
                var phiRes = CalculatePhiandPhiDer(qrIni.Q, qrIni.R, S, TrRadius);
                 LamdaLow = -(phiRes.Phi / phiRes.PhiDerivative);
                 lamdaUp = g.Norm() / TrRadius;
            }
            else
            {
                LamdaLow = -(Qprev / QderPrev);
                lamdaUp = g.Norm() / TrRadius;
            }
           
            bool exitCondition = false;
            int iter = 0;
            do
            {
                iter++;
                if (lamda < LamdaLow || lamda > lamdaUp) { lamda = Sqrt(LamdaLow * lamdaUp) > 1e-3 * lamdaUp ? Sqrt(LamdaLow * lamdaUp) : 1e-3 * lamdaUp; }
                 Hmod = PDM(H,lamda);
                Hmod = settings.MatrixModificationStrategy.SPD(Hmod);
                S = -LinearSolvers.Solve(Hmod, g);
                double stepLen = S.Norm();
                lamdaPrev = lamda;
                exitCondition = stepLen >= low * TrRadius && stepLen <= hi * TrRadius;
                if (!exitCondition)
                {
                    qr = QrFac(Hmod);
                  var  phiRes = CalculatePhiandPhiDer(qr.Q, qr.R, S, TrRadius);
                    double muo = lamda - (phiRes.Phi / phiRes.PhiDerivative);
                    LamdaLow = LamdaLow >= muo ? LamdaLow : muo;
                    if (phiRes.Phi < 0.0) { lamdaUp = lamda; }
                    lamda = Landa(qr.Q, qr.R, S, TrRadius, lamda);
                    lamdaPrev = lamda;
                    Qprev = phiRes.Phi;
                    QderPrev = phiRes.PhiDerivative;
                }
                if (iter > 10) break;

            } while (!exitCondition);
            return (S, lamda, TrRadius);
        }
        private static (Vector x,double fx, bool callBack) TrustRegionUpdate(Func<Vector ,double>fun, Vector xc,Vector sdir, Matrix H, Vector g,ref double TrR,ref Vector xprev,ref double fxPrev, bool newtonLenTaken,OptimizationAndSolverSettings setting)
        {
            Func<Vector,double,double, Vector, bool> AcceptanceTest = new Func<Vector, double, double, Vector, bool>((xn,fx,fxnect, xnn) =>
              {
                  double curve = setting.WolfeConstants.C1 * g.Dot(xnn);
                  bool accept = fxnect <= fx + curve;
                  return accept;
              });
            double steplen = sdir.Norm();
            Vector xnext = xc + sdir;
            Vector xdiff = xnext - xc;
            double fxc = fun(xc);
            double fnext = fun(xnext);
            double deltaf = fnext - fxc;
            double initialSlope = g.Dot(sdir);
            bool returnBack = false;
            if (!AcceptanceTest(xnext, fxc, fnext,xdiff))
            {
                Vector sScaled = sdir / xnext.Max<double>();
                double relLength = sScaled.Max<double>();
                if (relLength < setting.Tol)
                {
                    returnBack = false;
                    return (xc, fxPrev, returnBack);
                }
                else
                {
                    returnBack = true;
                    double dterm = (-initialSlope * steplen) / (2 * (deltaf - initialSlope));
                    if (dterm < 0.1 * TrR)
                    {
                        TrR = 0.1 * TrR;
                    }
                    else if (dterm > 0.5 * TrR)
                    {
                        TrR = 0.5 * TrR;
                    }
                    else
                    {
                        TrR = dterm;
                    }
                    return (xnext, fnext, returnBack);
                }
            }
            else
            {
                double shs = (double) (~sdir * (H * sdir));
                double deltafPred = initialSlope + shs;
                bool deltafMinusdeltafPred = Abs(deltafPred - deltaf) <= 0.1 * Abs(deltaf);
                bool deltaflessThanInitialslope = deltaf <= initialSlope;
                bool trradiusOkay = TrR <= 0.99 * setting.MaxStep;
                if((deltafMinusdeltafPred ||deltaflessThanInitialslope)&&trradiusOkay)
                {
                    returnBack = true;
                    xprev = xnext;
                    fxPrev = fnext;
                    TrR = (2.0 * TrR)<setting.MaxStep ?2.0*TrR:setting.MaxStep;
                    return (xnext, fnext, returnBack);
                }
                else
                {
                    returnBack = false;
                    if (deltaf >= 0.1 * deltafPred) { TrR = TrR / 2; }
                    else if (deltaf <= 0.75 * deltafPred) { TrR = 2 * TrR; }
                    else { TrR = 1.0 * TrR; }
                    return (xnext, fnext, returnBack);
                }

            }
           
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="Sn"></param>
        /// <param name="g"></param>
        /// <param name="X"></param>
        /// <param name="H"></param>
        /// <param name="Trradius"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector TrustRegionDriver(Func<Vector,double>fun,Vector Sn,Vector g,Vector X, Matrix H,ref double Trradius, OptimizationAndSolverSettings setting)
        {

            bool goBack = true;
            Vector xvec = null!;
            bool firstHook = true;
            double Lamda = 0.0;
            double LamdaPrev = 0.0;
            double Fx = 0.0;
            double Fxprev = 0.0;
            double Qprev = 0.0;
            double QderPrev = 0.0;
            do
            {
                var hookres = Hook(H, g, Sn, Trradius,ref LamdaPrev,ref Qprev,ref QderPrev, firstHook,setting);
                firstHook = false;
                 Lamda = hookres.lan;
                Vector sdir = hookres.Sdir;
                double Trr = hookres.Tr;
                Trradius = Trr;
                var TrUpdate = TrustRegionUpdate(fun, X, sdir, H, g, ref Trradius, ref X, ref Fxprev, true,setting);
                 goBack = TrUpdate.callBack;
                Fx = TrUpdate.fx;
                Fxprev = Fx;
                xvec = TrUpdate.x;
                Sn = sdir.Copy();
            } while (goBack);
            return xvec;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="UpperLimit"></param>
        /// <param name="LowerLimit"></param>
        /// <param name="NumberOfIntervals"></param>
        /// <returns></returns>
        public static Vector LinSpace(double UpperLimit,double LowerLimit,double NumberOfIntervals)
        {
            double Ul = UpperLimit > LowerLimit ? UpperLimit : LowerLimit;
            double Ll = LowerLimit < UpperLimit ? LowerLimit : UpperLimit;
            int nn = (int)Floor(NumberOfIntervals);
            Vector res = new Vector(nn)
            {
                [0] = Ll
            };
            double h = (Ul - Ll) / (nn - 1);
            for (int i = 1; i < res.Length; i++)
            {
                res[i] = Ll + (h * i);
            }

            return ~res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="UpperLimit"></param>
        /// <param name="LowerLimit"></param>
        /// <param name="NumberOfIntervals"></param>
        /// <returns></returns>
        public static Vector LogSpace(double UpperLimit, double LowerLimit, double NumberOfIntervals)
        {
            double Ul = UpperLimit > LowerLimit ? UpperLimit : LowerLimit;
            double Ll = LowerLimit < UpperLimit ? LowerLimit : UpperLimit;
            int nn = (int)Floor(NumberOfIntervals);
            Vector res = new Vector(nn)
            {
                [0] =Pow(10, Ll)
            };
            double h = (Ul - Ll) / (nn - 1);
            for (int i = 1; i < res.Length; i++)
            {
                res[i] =Pow(10, Ll + (h * i));
            }

            return ~res;
        }
        
    }
}
