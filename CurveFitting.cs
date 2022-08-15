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
    public enum LMUpdateStrategy
    {
        /// <summary>
        /// 
        /// </summary>
        LM,
        /// <summary>
        /// 
        /// </summary>
        Quadratic,
        /// <summary>
        /// 
        /// </summary>
        Nielsen,
        /// <summary>
        /// 
        /// </summary>
        TrustRegion
    }
    /// <summary>
    /// 
    /// </summary>
   public class CurveFitting
    {
        #region Linear Least Square
        /// <summary>
        /// No Normal Equation is formed, instead the algorithm performs QR decomposition of the Vandermonde matrix using the Householder method as default. 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="degree"></param>
        /// <param name="qrsch"></param>
        /// <returns></returns>
        public static Vector PolyFit(Vector x,Vector y,int degree=2,QRFactorizationMethod qrsch=QRFactorizationMethod.HouseHolder)
        {
            Matrix VanderMondeMat = Matrix.VandermondeMatrix(x, degree);
            var Ls = LinearSolvers.LeastSquare(VanderMondeMat, y, qrsch);
            Vector coeffs = Ls.OptimalVector;
            Vector ypred = new Vector();
            for (int i = 0; i < x.Length; i++)
            {
                double yp = Polynomial.polyEval(coeffs, x[i]);
                ypred.Add(yp);
            }
            var errorAna = ErrorAnalysis(x, y, ypred, coeffs);
            return coeffs;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="degree"></param>
        /// <param name="colToModify"></param>
        /// <param name="ModificationCallBack"></param>
        /// <returns></returns>
        public static Vector PolyFit(Vector x, Vector y, int degree,int[]colToModify, List<Func<Vector>>ModificationCallBack =null!)
        {
            Matrix VanderMondeMat = Matrix.VandermondeMatrix(x, degree);
            VanderMondeMat.Modify(colToModify, ModificationCallBack);
            var Ls = LinearSolvers.LeastSquare(VanderMondeMat, y, QRFactorizationMethod.HouseHolder);
            Vector coeffs = Ls.OptimalVector;
            Vector ypred = new Vector();
            for (int i = 0; i < x.Length; i++)
            {
                double yp = Polynomial.polyEval(coeffs, x[i]);
                ypred.Add(yp);
            }
            var errorAna = ErrorAnalysis(x, y,ypred, coeffs);
            return coeffs;
        }
        /// <summary>
        /// This version generates the normal matrix from the vondermonde matrix and set up the normal equation.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        public static Vector PolyFit2(Vector x, Vector y, int degree=2)
        {
            Matrix VanderMondeMat = Matrix.VandermondeMatrix(x, degree);
            Matrix Normal = ~VanderMondeMat * VanderMondeMat;
            Vector NormalVec = ~VanderMondeMat * y;
            Vector coeffs = LinearSolvers.Solve(Normal, NormalVec);
            Vector ypred = new Vector();
            for (int i = 0; i < x.Length; i++)
            {
                double yp = Polynomial.polyEval(coeffs, x[i]);
                ypred.Add(yp);
            }
            var errorAna = ErrorAnalysis(x, y, ypred, coeffs);
            return coeffs;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="ypred"></param>
        /// <param name="coeffs"></param>
        /// <returns></returns>
        public  static (double,double RSquared,double correlation)ErrorAnalysis(Vector x,Vector y,Vector ypred, Vector coeffs)
        {
            double ymean = y.Average();
            Vector ydev = (y - ymean)^2;
            double sumYdev = ydev.Sum();
           
            int n = 0;
            if (coeffs.Length == 2) {n= x.Length - 2; }
            else { n = x.Length - (coeffs.Length); }
            Vector ydiff = y - ypred;
            Vector R_squared = ydiff^2;
            double Sy = Sqrt(sumYdev / (y.Length - 1));
            double SumR_sqr = R_squared.Sum();
            double sRe = Sqrt(SumR_sqr / n);
            double coeffOfDet = (sumYdev - SumR_sqr) / sumYdev;
            double corr = Sqrt(coeffOfDet);
            return (sRe, coeffOfDet, corr);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="numDegree"></param>
        /// <param name="denomDegree"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static (Vector coeffs,Vector predictedVals)RationalFit(Vector x,Vector y,int numDegree,int denomDegree)
        {
            if ((numDegree + denomDegree) > x.Length) { throw new Exception("The sum of the degree of the polynomials must be less than or equal to the number of data point to be fitted."); }
            Matrix ModifiedVandermontMat = new Matrix();
            for (int i = 0; i <= numDegree; i++)
            {
                Vector C = x.Pow(i);
                ModifiedVandermontMat.Add(C);
            }
            for (int i = 1; i <= denomDegree ; i++)
            {
                Vector C = y * x.Pow(i);
                ModifiedVandermontMat.Add(-C);
            }
            var Ls = LinearSolvers.LeastSquare(ModifiedVandermontMat, y, QRFactorizationMethod.HouseHolder);
            Vector coeffs = Ls.OptimalVector;
           
            Vector ypred = new Vector();
            for (int i = 0; i < x.Length; i++)
            {
                double numterm = 0.0;
                double denterm = 0.0;
                for (int j = 0; j < coeffs.Length; j++)
                {
                    if (j <= numDegree)
                    {
                        numterm += coeffs[j] * Pow(x[i], j);
                    }
                    else
                    {
                        denterm += coeffs[j] * Pow(x[i], j-numDegree);
                    }
                }
                double yval = numterm /(1.0+ denterm);
                ypred.Add(yval);
            }
            var errorAna = ErrorAnalysis(x, y,ypred, coeffs);
            return (coeffs, ypred);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="EstP"></param>
        /// <returns></returns>
        public static (Vector coeff,Vector FuncValues) ExponetialFit(Vector x, Vector y,Vector EstP)
        {
           
            Vector Trimy = new Vector();
            Vector trimx = new Vector();
            if (y.Any(d => d <= 0))
            {
                for (int i = 0; i < y.Length; i++)
                {
                    if (y[i] > 0)
                    {
                        trimx.Add( x[i]);
                        Trimy.Add(y[i]);
                    }
                }
            }
            else
            {
                trimx = x.Copy();
                Trimy = y.Copy();
            }
            Vector lny =(Vector) Trimy.Select(item => Log(item));
            Vector coeffs = PolyFit(trimx, lny, 1);
            coeffs[0] = Exp(coeffs[0]);
            Vector res = new Vector(EstP.Length);
            for (int i = 0; i < EstP.Length; i++)
            {
                res[i] = coeffs[0] * Exp(coeffs[1] * EstP[i]);
            }

            return (coeffs,res);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="EstP"></param>
        /// <returns></returns>
        public static (Vector coeff, Vector FuncValues) PowerFit(Vector x, Vector y, Vector EstP)
        {
            Vector Trimy = new Vector();
            Vector trimx = new Vector();
            if (y.Any(d => d <= 0))
            {
                for (int i = 0; i < y.Length; i++)
                {
                    if (y[i] > 0)
                    {
                        trimx.Add( x[i]);
                        Trimy.Add(y[i]);
                    }
                   
                }
            }
            if (x.Any(d => d <= 0))
            {
                for (int i = 0; i < y.Length; i++)
                {
                    if (x[i] > 0)
                    {
                        trimx.Add(x[i]);
                        Trimy.Add(y[i]);
                    }
                }
            }
            if(!(y.Any(d => d <= 0)) && !( x.Any(d => d <= 0)))
            {
                trimx = x.Copy();
                Trimy = y.Copy();
            }
            
            Vector logy =(Vector) Trimy.Select(item => Log10(item));
            Vector logx =(Vector) trimx.Select(item => Log10(item));
            Vector coeffs = PolyFit(logx, logy, 1);
            coeffs[0] = Pow(10,coeffs[0]);
            Vector res = new Vector(EstP.Length);
            for (int i = 0; i < EstP.Length; i++)
            {
                res[i] = coeffs[0] * Pow(EstP[i],coeffs[1]);
            }
            return (coeffs,res);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="IndependentVariables"></param>
        /// <param name="dependenrVariable"></param>
        /// <param name="EstPoint"></param>
        /// <param name="includeProductterm"></param>
        /// <returns></returns>
        public static (Vector coeff,Vector FuncValues)  MultiLinearRegression(List<Vector>IndependentVariables,Vector dependenrVariable,List<Vector> EstPoint,bool includeProductterm=false)
        {
            Matrix vander = Matrix.VandermondeMatrix(IndependentVariables,includeProductterm);
            var Ls = LinearSolvers.LeastSquare(vander, dependenrVariable, QRFactorizationMethod.HouseHolder);
            Vector coeffs = Ls.OptimalVector;
            Vector res = Vector.Zeros(EstPoint[0].Length);
            Vector prodT = Vector.Ones(EstPoint[0].Length);
            for (int i = 0; i < EstPoint.Count; i++)
            {
                if (includeProductterm) { prodT *= EstPoint[i]; }
              res+=  coeffs[i+1] * EstPoint[i];
            }
            res +=  coeffs[0];
            if (includeProductterm) {res+= coeffs.Last() * prodT; }
            return (coeffs,res);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="xlist"></param>
        /// <param name="yMat"></param>
        /// <param name="FirstIndepVar"></param>
        /// <param name="secondIndVar"></param>
        /// <returns></returns>
        public static (Vector coeff, Matrix funcval) QuadraticBivariateFit(Matrix xlist, Matrix yMat,Vector FirstIndepVar,Vector secondIndVar)
        {
            Matrix vander = new Matrix();
            vander.Add(Vector.Ones(yMat.Ncol*yMat.Nrow));
            Vector C1 = new Vector(); Vector C2 = new Vector(); Vector C3 = new Vector();
            Vector C4 = new Vector(); Vector C5 = new Vector();
            for (int i = 0; i < xlist[0].Length; i++)
            {
                for (int j = 0; j < xlist[0].Length; j++)
                {
                    C1.Add(xlist[0,i]);C2.Add(xlist[1,j]);C3.Add(xlist[0,i] * xlist[1,j]);
                    C4.Add(Pow(xlist[0,i], 2));C5.Add(Pow(xlist[1,j],2));
                }
            }
            vander.Add(C1); vander.Add(C2); vander.Add(C3); vander.Add(C4); vander.Add(C5);
            Vector y = Matrix.vec(yMat);
            var Ls = LinearSolvers.LeastSquare(vander, y, QRFactorizationMethod.HouseHolder);
            Vector coeffs = Ls.OptimalVector;
            Matrix res = new Matrix(FirstIndepVar.Length, FirstIndepVar.Length);
            for (int i = 0; i < FirstIndepVar.Length; i++)
            {
                for (int j = 0; j < secondIndVar.Length; j++)
                {
                    res[j,i] = coeffs[0] + coeffs[1] * FirstIndepVar[i] + coeffs[2] * secondIndVar[j] + coeffs[3] * FirstIndepVar[i] * secondIndVar[j] +
                        coeffs[4] * Pow(FirstIndepVar[i], 2) + coeffs[5] * Pow(secondIndVar[j], 2);
                }
            }
           
            return (coeffs,res);
        }
        #endregion
        #region Non Linear Least Square
        /// <summary>
        /// Buggy do not use
        /// </summary>
        /// <param name="fun"></param>
        /// <param name="InitialGuess"></param>
        /// <param name="setting"></param>
        /// <returns></returns>
        public static Vector GaussNewton(Func<Vector,Vector>fun,Vector InitialGuess,OptimizationAndSolverSettings setting=null!)
        {
            if (setting == null) { setting = new OptimizationAndSolverSettings(); }
            Vector x = InitialGuess.Copy();
            Vector r = fun(x);
            Func<Vector, double> objFunc = new Func<Vector, double>((k) => 0.5 * k.Dot(k));
            Matrix J = BMath.Jacobian(fun, x, setting.Delta);
            Matrix H = ~J *J;
            Vector g = ~J * r;
            double alpha = 1.0;
            int iter = 0;
            bool exitCriteria = true;
            do
            {
                iter++;
                Vector Dsearch = LinearSolvers.Solve(H, -g);
                if (BMath.WolfeCondition(objFunc, x, Dsearch, 0.5, setting.WolfeConstants.C1, setting.WolfeConstants.C2, setting.WolfeConstants.IncludeCurvatureCondition))
                {
                    alpha = BMath.InterpolationLineSearch(objFunc, x, Dsearch, setting);
                }
                Vector ScaledDsearch = alpha * Dsearch;
                x = x + ScaledDsearch;
                r = fun(x);
                double N = r.Norm();
                exitCriteria = N >= setting.Tol;
                J = BMath.Jacobian(fun, x, setting.Delta);
                H = ~J * J;
                g = ~J * r;
                if (iter > setting.MaxIteration) break;
            } while (exitCriteria);

            return x;
        }
        /// <summary>
        /// buggy do not use
        /// </summary>
        /// <param name="ypredicted"></param>
        /// <param name="ymeasured"></param>
        /// <param name="P0"></param>
        /// <param name="settings"></param>
        /// <returns></returns>
        public static LMResults LevenbergMarquart(Func<Vector,Vector>ypredicted,Vector ymeasured,Vector P0,LMSetting settings)
        {
            LMResults result = new LMResults();
            Vector Pcurrent = P0.Copy();
            Vector Pnext = P0.Copy();
            Matrix weight = settings.Weight;
            int Np = Pcurrent.Length;
            int Ndp = ymeasured.Length;
            int DOF = Ndp - Np + 1;
            if (weight == null!)
            {
                double val =  (double)(~ymeasured * ymeasured);
                double diaVal = 1.0 / val;
                weight = Matrix.Diagonal(Ndp, diaVal);
            }
            Func<Vector, Vector> rFunc = new Func<Vector, Vector>(parameters =>
                {
                    Vector y = ypredicted(parameters);
                    Vector res = ymeasured - y;
                    return res;
                });
            Func<Vector, double> CostFunc = new Func<Vector, double>(resVector =>
                {
                    //Vector resVector = rFunc(parameters);
                    double res = (double)(~resVector *(weight* resVector));
                    return res;
                });

            Matrix jac = BMath.Jacobian(rFunc, Pcurrent);
            Matrix Hessian = ~jac * (weight * jac);
            // MATRIX eye = MATRIX.IdentityMatrix(Hessian.Nrow);
            Matrix DtD = Matrix.IdentityMatrix(Hessian.Nrow); 
            double maxDiaVal = 0.0;
            double lamda = settings.Dampingfactor;
            switch (settings.LMUpdateMode)
            {
                case LMUpdateStrategy.LM:
                    lamda = settings.Dampingfactor;
                    
                    break;
                case LMUpdateStrategy.Quadratic:
                     maxDiaVal = Hessian.GetDiagonalElements().Max<double>();
                    lamda = lamda * maxDiaVal;
                    break;
                case LMUpdateStrategy.Nielsen:
                    
                    maxDiaVal = Hessian.GetDiagonalElements().Max<double>();
                    lamda = lamda * maxDiaVal;
                    break;
                case LMUpdateStrategy.TrustRegion:
                    return LevenbergMarquartTR(ypredicted, ymeasured, P0, settings);
                default:
                    break;
            }
            Vector residual = rFunc(Pcurrent); 
            double currentCost = CostFunc(residual);
            double nextCost = currentCost;
            List<double> Costs = new List<double>();
            List<Vector> deltapList = new List<Vector>();
            List<Vector> geodesicList = new List<Vector>();
            List<bool> geodesicUpdateAcceptedList = new List<bool>();
            bool exitcriteria = true;
            int counter = 0;
            double alp = 0.0;
            double nu = settings.Nu;
            Matrix HessianModified = null!;
            Vector gradient = ~jac * (weight * residual); 
            Vector residualnext = residual.Copy();
            Vector ypredC = ypredicted(Pcurrent);
            Vector ypredNext = ypredC.Copy();
            double rho = 0.0;
            Vector deltaP = null!;
            Vector acc = null!;
            double actualReduction = 0.0;
            double predictedReduction = 0.0;
            bool geodesicUpdateAccepted = false;
            bool gradientConverged = false;
            bool ParametersConvreged = false;
            bool CostConverged = false;
            result.Costs = new Vector();
            do
            {
                counter++;
                result.Costs.Add(currentCost);
                List<Vector> ContourList = new List<Vector>();
                for (int i = 0; i < Pcurrent.Length; i++)
                {
                    Vector v = new Vector();
                    ContourList.Add(v);
                }
                for (int i = 0; i < Pcurrent.Length; i++)
                {
                    ContourList[i].Add(Pcurrent[i]);
                }
                result.ContourData.C.Add(currentCost);
                result.ContourData.P.Add(ContourList);
                switch (settings.LMUpdateMode)
                {
                    case LMUpdateStrategy.LM:
                        HessianModified = Hessian + (lamda * DtD);
                        break;
                    case LMUpdateStrategy.Quadratic:
                        for (int i = 0; i < DtD.Nrow; i++)
                        {
                            if (Abs(Hessian[i, i]) > Abs(DtD[i, i])) { DtD[i, i] = Abs(Hessian[i, i]); }
                        }
                        HessianModified = Hessian + lamda * DtD;
                        break;
                    case LMUpdateStrategy.Nielsen:
                        for (int i = 0; i < DtD.Nrow; i++)
                        {
                            if (Abs(Hessian[i, i]) > Abs(DtD[i, i])) { DtD[i, i] = Abs(Hessian[i, i]); }
                        }
                        HessianModified = Hessian + lamda * DtD;
                        break;
                    default:
                        break;
                }
                deltaP = -LinearSolvers.Solve(HessianModified, gradient);
                if (settings.ApplyGeodesicAcceleration)
                {
                    Vector SecondOrderDir = SecondOrderDirectionalDerivative(Pcurrent, rFunc, jac, residual, deltaP,settings.SecondOrderDirectionDerivativeDelta);
                    Vector accGradient = ~jac * SecondOrderDir;
                    acc = -LinearSolvers.Solve(HessianModified, accGradient);
                    acc = 0.5 * acc;
                    double geoNorm = acc.Norm();
                    double deltapNorm = deltaP.Norm();
                    double normRation = (2.0 * geoNorm) / deltapNorm;
                    geodesicUpdateAccepted = normRation <= settings.GeodesicFactor;
                    if (geodesicUpdateAccepted)
                    {
                        geodesicList.Add(acc);
                        deltaP = deltaP + acc;
                    }
                }
                Pnext = Pcurrent + deltaP;
                if (settings.ApplyParameterBounds)
                {
                    for (int i = 0; i < Pnext.Length; i++)
                    {
                        if (Pnext[i] < settings.Lbound) { Pnext[i] = settings.Lbound; }
                        if (Pnext[i] > settings.Ubound) { Pnext[i] = settings.Ubound; }
                    }
                }
                residualnext = rFunc(Pnext);
                nextCost = CostFunc(residualnext);
                ypredNext = ypredicted(Pnext);
                if (settings.LMUpdateMode == LMUpdateStrategy.Quadratic)
                {
                     alp = (double)(~gradient * deltaP) / (((nextCost - currentCost) / 2.0) + (2.0 * ((double)(~gradient * deltaP))));
                    deltaP = alp * deltaP;
                    Pnext = Pcurrent + deltaP;
                    if (settings.ApplyParameterBounds)
                    {
                        for (int i = 0; i < Pnext.Length; i++)
                        {
                            if (Pnext[i] < settings.Lbound) { Pnext[i] = settings.Lbound; }
                            if (Pnext[i] > settings.Ubound) { Pnext[i] = settings.Ubound; }
                        }
                    }
                    residualnext = rFunc(Pnext);
                    nextCost = CostFunc(residualnext);
                    ypredNext = ypredicted(Pnext);
                }
                actualReduction = currentCost - nextCost;
                if (settings.LMUpdateMode == LMUpdateStrategy.LM)
                {
                    predictedReduction = (double)(~deltaP * (lamda * deltaP + gradient));
                    //predictedReduction = (double)(~deltaP * (lamda *( MATRIX.Diagonal(Hessian) * deltaP) + gradient));
                    rho = actualReduction / predictedReduction;
                }
                else
                {
                    predictedReduction = (double)(~deltaP * (lamda * (Matrix.Diagonal(Hessian) * deltaP) + gradient));
                    //predictedReduction = (double)(~deltaP * (lamda * deltaP + gradient));
                    rho = actualReduction / predictedReduction;
                }
                if (rho > settings.StepAcceptanceTol)
                {
                    AcceptStep(settings, ref Pcurrent, Pnext, weight, rFunc, ref jac, out Hessian, ref lamda, ref residual, currentCost, Costs, deltapList, counter, alp, ref nu, out gradient, residualnext, ref ypredC, ypredNext, rho, deltaP);
                    geodesicUpdateAcceptedList.Add(geodesicUpdateAccepted);
                }
                else
                {
                    if (settings.AllowUphillSteps&& deltapList.Count!=0)
                    {
                        if(settings.ApplyGeodesicAcceleration&& geodesicUpdateAccepted)
                        {
                            deltaP = deltaP - acc;
                           
                        }
                        Vector deltaPold = deltapList.Last();

                        double theta = deltaP.CosTheta(deltaPold);
                        double costPrev = Costs.Last();
                        bool acceptUphillStep = Pow((1.0 - theta), 2) * nextCost <= costPrev;
                        if (acceptUphillStep)
                        {
                            AcceptStep(settings, ref Pcurrent, Pnext, weight, rFunc, ref jac, out Hessian, ref lamda, ref residual, currentCost, Costs, deltapList, counter, alp, ref nu, out gradient, residualnext, ref ypredC, ypredNext, rho, deltaP);
                           
                        }
                        else
                        {
                            RejectStep(settings, ref lamda, currentCost, nextCost, alp, ref nu, HessianModified);
                        }

                    }
                    else
                    {
                        RejectStep(settings, ref lamda, currentCost, nextCost, alp, ref nu, HessianModified);
                    }
                }
                currentCost = nextCost;
                 gradientConverged = gradient.Norm() < settings.GradTol;
                 ParametersConvreged = (deltaP / (Pcurrent.Absolute() + 1e-12)).Norm() < settings.ParameterTol;
                 CostConverged = nextCost / (DOF) < settings.CostTol;
                exitcriteria = gradientConverged || ParametersConvreged || CostConverged;
                if (counter > settings.MaxIter) break;
            }
            while (!exitcriteria);
            double CostErrorCritierionvv =  currentCost / DOF;
            double CostErrorCritierion =(double)(~residual*(weight*residual)) / DOF;
            Matrix ParameterCov = Matrix.Inverse(Hessian);
            Vector AsymptoticStandardParameterErrors = Matrix.GetDiagonalElements(ParameterCov).Absolute().Pow(0.5);
            Matrix Corr = Matrix.TwDiv(ParameterCov, Matrix.OuterProduct(AsymptoticStandardParameterErrors, ~AsymptoticStandardParameterErrors));
            Matrix hhh = jac * ParameterCov * ~jac;
            Vector AsymptoticStandardErrorOfFit = new Vector();
            for (int i = 0; i < Ndp; i++)
            {
                Vector jr = jac[i, "row"];
                AsymptoticStandardErrorOfFit.Add((double)(jr * (ParameterCov * ~jr)));
            }
            AsymptoticStandardErrorOfFit = AsymptoticStandardErrorOfFit.Pow(0.5);
            Matrix mat = null!;
            Vector AsymptoticStandardPredictionError = null!;
            if (settings.Weight == null! || settings.Weight == Matrix.IdentityMatrix(Ndp))
            {
                double variance =(double) (~residual * residual) / DOF;
                mat = Matrix.Diagonal(Ndp, (1.0 / variance));
                AsymptoticStandardPredictionError = Matrix.GetDiagonalElements(mat+hhh).Absolute().Pow(0.5);
            }
            else
            {
                mat = weight.Copy();
                AsymptoticStandardPredictionError = Matrix.GetDiagonalElements(mat+hhh).Absolute().Pow(0.5);
            }
            //collect results
           
            result.OptimumParameters = Pcurrent.Copy();
            result.CorrelationMatrix = Corr.Copy();
            result.Jacobian = jac.Copy();
            result.Gradient = gradient.Copy();
            result.NumberOfIteration = counter;
            result.ConvergenceAchieved = exitcriteria;
            result.ReducedErrorCritirion = CostErrorCritierion;
            result.StandardParameterError = AsymptoticStandardParameterErrors;
            result.StandardErrorOfFit = AsymptoticStandardErrorOfFit;
            result.PercentParameterFit = (AsymptoticStandardParameterErrors / Pcurrent) * 100;
            if (exitcriteria)
            {
                result.Info = $"Convergence Achieved after {counter} iteration, dou to the following criteria been met:" +
                    $"Gradient Norm is {gradientConverged}, Parameter Norm is {ParametersConvreged} " +
                    $"Objective Function Criterion is {CostConverged}";
            }
            else
            {
                result.Info = $"Fail to converged because the maximum number {settings.MaxIter} of iteration is exceeded";
            }
           
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ypredicted"></param>
        /// <param name="ymeasured"></param>
        /// <param name="P0"></param>
        /// <param name="settings"></param>
        /// <returns></returns>
        public static LMResults LevenbergMarquartTR(Func<Vector, Vector> ypredicted, Vector ymeasured, Vector P0, LMSetting settings)
        {
            LMResults result = new LMResults();
            Vector Pcurrent = P0.Copy();
            Vector Pnext = P0.Copy();
            Matrix weight = settings.Weight;
            int Np = Pcurrent.Length;
            int Ndp = ymeasured.Length;
            int DOF = Ndp - Np + 1;
            if (weight == null!)
            {
                double val = (double)(~ymeasured * ymeasured);
                double diaVal = 1.0 / val;
                weight = Matrix.Diagonal(Ndp, diaVal);
            }
            Func<Vector, Vector> rFunc = new Func<Vector, Vector>(parameters =>
            {
                Vector y = ypredicted(parameters);
                Vector res = ymeasured - y;
                return res;
            });
            Func<Vector, double> CostFunc = new Func<Vector, double>(resVector =>
            {
                double res = (double)(~resVector * (weight * resVector));
                return res;
            });

            Matrix jac = BMath.Jacobian(rFunc, Pcurrent);
            Matrix Hessian = ~jac * (weight * jac);
            Matrix DtD = Matrix.IdentityMatrix(Hessian.Nrow);
            double maxDiaVal = Hessian.GetDiagonalElements().Max<double>();
            double lamda = settings.Dampingfactor*maxDiaVal;
            Vector residual = rFunc(Pcurrent);
            double currentCost = CostFunc(residual);
            double nextCost = currentCost;
            List<double> Costs = new List<double>();
            List<Vector> deltapList = new List<Vector>();
            List<Vector> geodesicList = new List<Vector>();
            List<bool> geodesicUpdateAcceptedList = new List<bool>();
            bool exitcriteria = true;
            int counter = 0;
            Matrix HessianModified = null!;
            Vector gradient = ~jac * (weight * residual);
            Vector residualnext = residual.Copy();
            double rho = 0.0;
            Vector deltaP = null!;
            Vector acc = null!;
            double nu = settings.Nu;
            double actualReduction = 0.0;
            double predictedReduction = 0.0;
            bool geodesicUpdateAccepted = false;
            bool gradientConverged = false;
            bool ParametersConvreged = false;
            bool CostConverged = false;
            result.Costs = new Vector();
            double Tr = settings.Trustregion;
            bool firstH = true; double Qprev = 0.0; double QderPrev = 0.0;
            
            do
            {
                counter++;
                result.Costs.Add(currentCost);
                List<Vector> ContourList = new List<Vector>();
                for (int i = 0; i < Pcurrent.Length; i++)
                {
                    Vector v = new Vector();
                    ContourList.Add(v);
                }
                for (int i = 0; i < Pcurrent.Length; i++)
                {
                    ContourList[i].Add(Pcurrent[i]);
                }
                result.ContourData.C.Add(currentCost);
                result.ContourData.P.Add(ContourList);
                for (int i = 0; i < DtD.Nrow; i++)
                {
                    if (Abs(Hessian[i, i]) > Abs(DtD[i, i])) { DtD[i, i] = Abs(Hessian[i, i]); }
                }
                HessianModified = Hessian + lamda * DtD;
                deltaP = -LinearSolvers.Solve(HessianModified, gradient);
                if (settings.ApplyGeodesicAcceleration)
                {
                    Vector SecondOrderDir = SecondOrderDirectionalDerivative(Pcurrent, rFunc, jac, residual, deltaP, settings.SecondOrderDirectionDerivativeDelta);
                    Vector accGradient = ~jac * SecondOrderDir;
                    acc = -LinearSolvers.Solve(HessianModified, accGradient);
                    acc = 0.5 * acc;
                    double geoNorm = acc.Norm();
                    double deltapNorm = deltaP.Norm();
                    double normRation = (2.0 * geoNorm) / deltapNorm;
                    geodesicUpdateAccepted = normRation <= settings.GeodesicFactor;
                    if (geodesicUpdateAccepted)
                    {
                        geodesicList.Add(acc);
                        deltaP = deltaP + acc;
                    }
                }
               
                var HookRes = BMath.Hook(Hessian, gradient, deltaP, Tr, ref lamda, ref Qprev, ref QderPrev, firstH, new OptimizationAndSolverSettings());
                deltaP = HookRes.Sdir;
                lamda = HookRes.lan;
                Tr = HookRes.Tr;
                Pnext = Pcurrent + deltaP;
                if (settings.ApplyParameterBounds)
                {
                    for (int i = 0; i < Pnext.Length; i++)
                    {
                        if (Pnext[i] < settings.Lbound) { Pnext[i] = settings.Lbound; }
                        if (Pnext[i] > settings.Ubound) { Pnext[i] = settings.Ubound; }
                    }
                }

                residualnext = rFunc(Pnext);
                nextCost = CostFunc(residualnext);
                actualReduction = currentCost - nextCost;
                Vector cxc = residual + jac * Pcurrent;
                 predictedReduction = cxc.Dot(cxc);
                double rNorm = residual.Norm();
                double rnextNorm = residualnext.Norm();
                Vector jp = jac * deltaP;
                double jpNorm = jp.Norm();
                Vector Dp = DtD * deltaP;
                double dpNorm = Dp.Norm();
                double numerator = 1.0 - Pow((rnextNorm / rNorm), 2);
                double numerator2 = (1.0 - nextCost/currentCost);
                double T1 = Pow((jpNorm / rNorm), 2);
                double T2 = 2 * Pow(((Sqrt(lamda) * dpNorm) / rNorm), 2);
                double denom = T1 + T2;
                rho = numerator / denom;

              double rho2 = actualReduction / (currentCost + predictedReduction);
               
                if (rho > settings.StepAcceptanceTol)
                {
                    AcceptStepTR(settings, ref Pcurrent, Pnext, weight, rFunc, ref jac, out Hessian, ref lamda, ref residual, currentCost, Costs, deltapList, counter, out gradient, residualnext, deltaP);
                    geodesicUpdateAcceptedList.Add(geodesicUpdateAccepted);
                    //bool newtonStepTaken = false;
                    if (lamda == 0)
                    {
                        //newtonStepTaken = true;
                        lamda = settings.Dampingfactor ;
                        lamda = new Vector(new double[] { lamda / settings.LamdaDownFac, 1e-7 }).Max<double>();
                    }
                    #region former trust region update scheme
                    //double temp = 0.0;
                    //if (nextCost > currentCost) { temp = 0.5; }
                    //else if (nextCost < 0.1 * currentCost) { temp = 0.1; }
                    //else
                    //{
                    //    double tau = (double)(~deltaP * ~jac * residual) / residual.Dot(residual);
                    //    temp = (0.5 * tau) / (tau + 0.5 * (1.0 - Pow((residualnext.Norm() / residual.Norm()), 2)));
                    //}
                    //if (temp < 0.1) { temp = 0.1; }
                    //if (temp > 0.5) { temp = 0.5; }
                    //Tr = Tr - temp * Tr;
                    //if (rho > 0.75)
                    //{
                    //    double norm = (DtD * Pcurrent).Norm();
                    //    Tr = 2 * norm;
                    //}
                    //if (rho > 0.25 && rho < 0.75 && !newtonStepTaken)// to be reviewd
                    //{
                    //    lamda = settings.Dampingfactor * maxDiaVal;
                    //    Tr = settings.Trustregion;
                    //}
                    //else if (rho > 0.25 && rho < 0.75 && newtonStepTaken)
                    //{
                    //    Tr = settings.Trustregion;
                    //}
                    #endregion
                }
                else
                {
                        lamda = nu* settings.Dampingfactor ;
                        nu = 2 * nu;
                        Tr =  0.5*Tr;
                }
                currentCost = nextCost;
                gradientConverged = gradient.Norm() < settings.GradTol;
                ParametersConvreged = (deltaP / (Pcurrent.Absolute() + 1e-12)).Norm() < settings.ParameterTol;
                CostConverged = nextCost / (DOF) < settings.CostTol;
                exitcriteria = gradientConverged || ParametersConvreged || CostConverged;
                if (counter > settings.MaxIter) break;
            }
            while (!exitcriteria);
            double CostErrorCritierionvv = currentCost / DOF;
            double CostErrorCritierion = (double)(~residual * (weight * residual)) / DOF;
            Matrix ParameterCov = Matrix.Inverse(Hessian);
            Vector AsymptoticStandardParameterErrors = Matrix.GetDiagonalElements(ParameterCov).Absolute().Pow(0.5);
            Matrix Corr = Matrix.TwDiv(ParameterCov, Matrix.OuterProduct(AsymptoticStandardParameterErrors, ~AsymptoticStandardParameterErrors));
            Matrix hhh = jac * ParameterCov * ~jac;
            Vector AsymptoticStandardErrorOfFit = new Vector();
            for (int i = 0; i < Ndp; i++)
            {
                Vector jr = jac[i, "row"];
                AsymptoticStandardErrorOfFit.Add((double)(jr * (ParameterCov * ~jr)));
            }
            AsymptoticStandardErrorOfFit = AsymptoticStandardErrorOfFit.Pow(0.5);
            Matrix mat = null!;
            Vector AsymptoticStandardPredictionError = null!;
            if (settings.Weight == null! || settings.Weight == Matrix.IdentityMatrix(Ndp))
            {
                double variance = (double)(~residual * residual) / DOF;
                mat = Matrix.Diagonal(Ndp, (1.0 / variance));
                AsymptoticStandardPredictionError = Matrix.GetDiagonalElements(mat + hhh).Absolute().Pow(0.5);
            }
            else
            {
                mat = weight.Copy();
                AsymptoticStandardPredictionError = Matrix.GetDiagonalElements(mat + hhh).Absolute().Pow(0.5);
            }
            //collect results

            result.OptimumParameters = Pcurrent.Copy();
            result.CorrelationMatrix = Corr.Copy();
            result.Jacobian = jac.Copy();
            result.Gradient = gradient.Copy();
            result.NumberOfIteration = counter;
            result.ConvergenceAchieved = exitcriteria;
            result.ReducedErrorCritirion = CostErrorCritierion;
            result.StandardParameterError = AsymptoticStandardParameterErrors;
            result.StandardErrorOfFit = AsymptoticStandardErrorOfFit;
            result.PercentParameterFit = (AsymptoticStandardParameterErrors / Pcurrent) * 100;
            if (exitcriteria)
            {
                result.Info = $"Convergence Achieved after {counter} iteration, due to the following criteria been met:" +
                    $"Gradient Norm is {gradientConverged}, Parameter Norm is {ParametersConvreged} " +
                    $"Objective Function Criterion is {CostConverged}";
            }
            else
            {
                result.Info = $"Fail to converged after the maximum number {settings.MaxIter} of iteration is exceeded";
            }

            return result;
        }

        private static Vector SecondOrderDirectionalDerivative(Vector Pcurrent, Func<Vector, Vector> rFunc, Matrix jac, Vector residual, Vector deltaP,double h)
        {
            Vector Increment = Pcurrent +h * deltaP;
            Vector Decrement = Pcurrent - h * deltaP;
            Vector residualIncrement = rFunc(Increment);
            Vector residualDecrement = rFunc(Decrement);
            //Vector rdiff = (residualInr - residual) / h;
            //Vector rprime = rdiff - (jac * deltaP);
            //rprime = (2.0 / 0.1) * rprime;
            Vector rprime = ((residualIncrement - 2*residual)+residualDecrement) / Pow(h,2);
            return rprime;
        }

        private static void RejectStep(LMSetting settings, ref double lamda, double currentCost, double nextCost, double alp, ref double nu, Matrix H =null!)
        {
            switch (settings.LMUpdateMode)
            {
                case LMUpdateStrategy.LM:
                    lamda = new Vector(new double[] { lamda * settings.LamdaUpFac, 1e7 }).Min<double>();
                    break;
                case LMUpdateStrategy.Quadratic:
                    double condest = H.Condest();
                    if (H.IsPositiveDefinite()&&condest<1e8)
                    {
                        lamda = lamda / settings.DecayRate;
                    }
                    else
                    {
                       if(H.IsDiagonallyDominant())
                        {
                            lamda = H.GetDiagonalElements().Min<double>();
                        }
                        else
                        {
                            AddingMultipleOfIdentityMatrix modH = new AddingMultipleOfIdentityMatrix();
                            Matrix HHH = modH.SPD(H);
                            lamda = HHH.GetDiagonalElements().Min<double>();
                        }

                    }
                    
                    //lamda = lamda +Abs( (nextCost - currentCost) / (2 * alp));
                    break;
                case LMUpdateStrategy.Nielsen:
                    lamda = lamda * nu;
                    nu = 2 * nu;
                    break;
                default:
                    break;
            }
        }

        private static void AcceptStep(LMSetting settings, ref Vector Pcurrent, Vector Pnext, Matrix weight, Func<Vector, Vector> rFunc, ref Matrix jac, out Matrix JtWJ, ref double lamda, ref Vector residual, double currentCost, List<double> Costs, List<Vector> deltapList, int counter, double alp, ref double nu, out Vector Jtr, Vector residualnext, ref Vector ypredC, Vector ypredNext, double rho, Vector deltaP)
        {

            Vector deltaR = residual - residualnext;
            Vector delP = Pcurrent - Pnext;

            //Vector deltaR = residualnext - residual;
            //Vector delP = Pnext - Pcurrent;
            double hh = Matrix.FrobeniusNorm(jac);
            Matrix jacOld = jac.Copy();
            if (settings.RankOneUpdateSettings.ApplyBrodenUpdate)
            {
                if (counter % settings.RankOneUpdateSettings.UpdateAfter != 0)
                {
                    Matrix jacFd = BMath.Jacobian(rFunc, Pcurrent);
                    jac = jac + (Matrix.OuterProduct((deltaR - (jac * delP)), ~delP)) / (double)(~delP * delP);//Broyden rank one update
                    bool vcvcv = jac == jacFd;
                    double gg = Matrix.FrobeniusNorm(jac);
                    double ttt = Matrix.FrobeniusNorm(jacFd);
                    double uu = Abs(gg - hh);
                    if (uu > 0.1)
                    {
                        jac = BMath.Jacobian(rFunc, Pcurrent);
                    }
                }
                else
                {
                    jac = BMath.Jacobian(rFunc, Pcurrent);
                }
            }
            else
            {
                jac = BMath.Jacobian(rFunc, Pcurrent);
            }
            switch (settings.LMUpdateMode)
            {
                case LMUpdateStrategy.LM:
                    lamda = new Vector(new double[] { lamda / settings.LamdaDownFac, 1e-7 }).Max<double>();
                    break;
                case LMUpdateStrategy.Quadratic:
                    double wip = Pnext.Dot(Pcurrent) / (Pnext.Norm() * Pcurrent.Norm());
                    lamda = lamda * Pow(settings.DecayRate, wip);
                    //lamda = new Vector(new double[] { lamda / (1 + alp), 1e-7 }).Max<double>();
                    break;
                case LMUpdateStrategy.Nielsen:
                    lamda = lamda * new Vector(new double[] { 1.0 / 3.0, 1 - Pow((2 * rho - 1), 3) }).Max<double>();
                    nu = settings.Nu;
                    break;
                default:
                    break;
            }
            Pcurrent = Pnext.Copy();
            JtWJ = ~jac * (weight * jac);
            residual = residualnext.Copy();// rFunc(Pcurrent);
            Jtr = ~jac * (weight * residual);
            ypredC = ypredNext.Copy();
            deltapList.Add(deltaP);
            Costs.Add(currentCost);
        }
        private static void AcceptStepTR(LMSetting settings, ref Vector Pcurrent, Vector Pnext, Matrix weight, Func<Vector, Vector> rFunc, ref Matrix jac, out Matrix JtWJ, ref double lamda, ref Vector residual, double currentCost, List<double> Costs, List<Vector> deltapList, int counter, out Vector Jtr, Vector residualnext, Vector deltaP)
        {

            Vector deltaR = residual - residualnext;
            Vector delP = Pcurrent - Pnext;

            //Vector deltaR = residualnext - residual;
            ////Vector delP = Pnext - Pcurrent;
            //double hh = MATRIX.FrobeniusNorm(jac);
            //MATRIX jacOld = jac.Copy();
            if (settings.RankOneUpdateSettings.ApplyBrodenUpdate)
            {
                if (counter % settings.RankOneUpdateSettings.UpdateAfter != 0)
                {
                    jac = jac + (Matrix.OuterProduct((deltaR - (jac * delP)), ~delP)) / (double)(~delP * delP);//Broyden rank one update
                    //double gg = MATRIX.FrobeniusNorm(jac);
                    //double uu = Abs(gg - hh);
                    //if (uu > 0.1)
                    //{
                    //    jac = BMath.Jacobian(rFunc, Pcurrent);

                    //}
                }
                else
                {
                    jac = BMath.Jacobian(rFunc, Pcurrent);
                }
            }
            else
            {
                jac = BMath.Jacobian(rFunc, Pcurrent);
            }

            Pcurrent = Pnext.Copy();
            JtWJ = ~jac * (weight * jac);
            residual = residualnext.Copy();// rFunc(Pcurrent);
            Jtr = ~jac * (weight * residual);
            deltapList.Add(deltaP);
            Costs.Add(currentCost);
        }
        private static double UpdateTrustRegion(double rho,double Tr,ref double  lambda,Vector deltap, Matrix D, Matrix jtJ, double Cc,double Cn)
        {
            double Pnorm = deltap.Dot(D * deltap);
            double Term = 0.0;
            if (rho>0.25)
            {
                if (lambda > 0.0 && rho < 0.75)
                {
                    Term = 1.0;
                }
                else
                {
                    Term = 2.0 * Pnorm / Tr;
                }

            }
            else
            {
                double actRed = (1.0 - Cn / Cc);
                double T1 = deltap.Dot(jtJ * deltap) / Cc;
                //double T2 = 0.5*lambda * deltap.Dot(D * deltap) / Cc;
                double T2 = Sqrt(lambda) * deltap.Dot(D * deltap) / Cc;
                double gamma = -1.0 * (T1 + T2);
                if (actRed >= 0)
                {
                    Term = 0.5;
                }
                else
                {
                    Term = (0.5 * gamma) / (gamma + 0.5 * actRed);
                }
                if (0.1 * Cn >= Cc||Term<0.1) { Term = 0.1; }
            }
            Vector vec = new double[] { Tr, 10 * Pnorm };
            double UpdatedTr = Term * vec.Max<double>();
            if (UpdatedTr < 1e-2) { UpdatedTr = 1e-2; }
            lambda = lambda / Term;
            return UpdatedTr;
        }
        #endregion

    }

    /// <summary>
    /// 
    /// </summary>
    public class CurveFitResults
    {
        private double standardError;
        private double correlation;
        private double coefficientOfDetermination;
        Vector parameters;
        double predictedValue;
        /// <summary>
        /// 
        /// </summary>
        public double StandardError { get => standardError; set => standardError = value; }
        /// <summary>
        /// 
        /// </summary>
        public double Correlation { get => correlation; set => correlation = value; }
        /// <summary>
        /// 
        /// </summary>
        public double CoefficientOfDetermination { get => coefficientOfDetermination; set => coefficientOfDetermination = value; }
        /// <summary>
        /// 
        /// </summary>
        public Vector Parameters { get => parameters; set => parameters = value; }
        /// <summary>
        /// 
        /// </summary>
        public double PredictedValue { get => predictedValue; set => predictedValue = value; }
    }
    /// <summary>
    /// 
    /// </summary>
    public class LMSetting
    {
        /// <summary>
        /// 
        /// </summary>
        public double Dampingfactor { get; set; } = 0.001;
        /// <summary>
        /// 
        /// </summary>
        public double Trustregion { get; set; } = 1.0;
        /// <summary>
        /// 
        /// </summary>
        public double GradTol { get; set; } = new OptimizationAndSolverSettings().Eps;
        /// <summary>
        /// 
        /// </summary>
        public double ParameterTol { get; set; } = new OptimizationAndSolverSettings().Eps;
        /// <summary>
        /// 
        /// </summary>
        public double CostTol { get; set; } = new OptimizationAndSolverSettings().Eps;
        /// <summary>
        /// 
        /// </summary>
        public double StepAcceptanceTol { get; set; } = 0.1;
        /// <summary>
        /// 
        /// </summary>
        public LMUpdateStrategy LMUpdateMode { get; set; } = LMUpdateStrategy.LM;
        /// <summary>
        /// 
        /// </summary>
        public double LamdaDownFac { get; set; } = 9.0;
        /// <summary>
        /// 
        /// </summary>
        public double LamdaUpFac { get; set; } = 11.0;
        /// <summary>
        /// 
        /// </summary>
        public double Nu { get; set; } = 2.0;
        /// <summary>
        /// 
        /// </summary>
        public double Ubound { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public double Lbound { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public int MaxIter { get; set; } = 500;
        /// <summary>
        /// 
        /// </summary>
       public Matrix Weight { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public bool ApplyGeodesicAcceleration { get; set; } = false;
        /// <summary>
        /// 
        /// </summary>
        public bool ApplyParameterBounds { get; set; } = false;
        /// <summary>
        /// 
        /// </summary>
        public bool AllowUphillSteps { get; set; } = false;
        /// <summary>
        /// 
        /// </summary>
        public double GeodesicFactor
        {
            get => geodesicFactor;
            set
            {
                if (value > 1.0) { value = 1.0; }
                if (value < 0.0) { value = 0.0; }
                geodesicFactor = value;
            }
        }
            

        private double geodesicFactor=0.75;
        /// <summary>
        /// 
        /// </summary>
        public double SecondOrderDirectionDerivativeDelta { get; set; } = 0.1;
        /// <summary>
        /// 
        /// </summary>
        public (bool ApplyBrodenUpdate,int UpdateAfter) RankOneUpdateSettings { get; set; } = (false,3);
        /// <summary>
        /// 
        /// </summary>
        public double DecayRate { get; set; } = 0.1;
       
    }
    /// <summary>
    /// 
    /// </summary>
    public class LMResults
    {
        /// <summary>
        /// 
        /// </summary>
        public Vector OptimumParameters { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Matrix CorrelationMatrix { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Vector StandardParameterError { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Vector StandardErrorOfFit { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Matrix Jacobian { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Vector Gradient { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Vector PercentParameterFit { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public bool ConvergenceAchieved { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public int NumberOfIteration { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public string Info { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public double ReducedErrorCritirion { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public Vector Costs { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public (Vector C, List<List<Vector>> P) ContourData { get; set; } = (new Vector(), new List<List<Vector>>());
    }
}
