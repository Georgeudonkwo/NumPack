using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumSharp.ODE
{
    /// <summary>
    /// 
    /// </summary>
   public  class SingleSteps
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix Euler(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy = results[i - 1, 1, 0, "row"];
                Vector k = h * model(meshSize[i - 1], yy);
                results[i, 1, 0, "row"] = yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix RK4(Func<double,Vector,Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy =~results[i - 1, 1, 0, "row"];
                double x =meshSize[i - 1];
                Vector k0 = h * model(x, yy);
                Vector k1 = h * model(x+h/2.0, yy+k0/2.0);
                Vector k2 = h * model(x + h / 2.0, yy + k1 / 2.0);
                Vector k3 = h * model(x + h , yy + k2);
                Vector k = 1.0 / 6.0 * (k0 + 2 * k1 + 2 * k2 + k3);
                results[i, 1, 0, "row"] = ~yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix RK2ModifiedEuler(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy = ~results[i - 1, 1, 0, "row"];
                double x = meshSize[i - 1];
                Vector k0 = h * model(x, yy);
                Vector k1 = h * model(x + h / 2.0, yy + k0*(h / 2.0));
                Vector k = h* k1 ;
                results[i, 1, 0, "row"] = ~yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix RK2Ralston(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy = ~results[i - 1, 1, 0, "row"];
                double x = meshSize[i - 1];
                Vector k0 = h * model(x, yy);
                Vector k1 = h * model(x +(3.0* h / 4.0), yy + k0*(3.0*h / 4.0));
                Vector k = h / 3.0 * (k0 + 2 * k1 );
                results[i, 1, 0, "row"] = ~yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix RKHeun(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy = ~results[i - 1, 1, 0, "row"];
                double x = meshSize[i - 1];
                Vector k0 = h * model(x, yy);
                Vector k1 = h * model(x + h , yy + h*k0);
                Vector k = h / 2.0 * (k0 + k1 );
                results[i, 1, 0, "row"] = ~yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <returns></returns>
        public static Matrix RK3(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();// new MATRIX(meshSize.Length,InitialVal.Length);
            results.Add(~meshSize);
            for (int i = 0; i < InitialVal.Length; i++)
            {
                Vector res = new Vector(meshSize.Length);
                res[0] = InitialVal[i];
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            for (int i = 1; i < meshSize.Length; i++)
            {
                Vector yy = ~results[i - 1, 1, 0, "row"];
                double x = meshSize[i - 1];
                Vector k0 = h * model(x, yy);
                Vector k1 = h * model(x + h / 2.0, yy + k0*(h / 2.0));
                Vector k2 = h * model(x + h , yy - k0*h+2*h*k1);
                Vector k = h / 6.0 * (k0 + 4 * k1 +  k2 );
                results[i, 1, 0, "row"] = ~yy + ~k;
            }
            return results;
        }
        /// <summary>
        /// not working
        /// </summary>
        /// <param name="model"></param>
        /// <param name="InitialVal"></param>
        /// <param name="Mesh"></param>
        /// <param name="tol"></param>
        /// <param name="hmin"></param>
        /// <param name="hmax"></param>
        /// <returns></returns>
        public static Matrix RK45(Func<double, Vector, Vector> model, Vector InitialVal, (double StartPoint, double Endpoint, int numberOfPoints) Mesh,double tol,double hmin,double hmax)
        {
            Vector meshSize = BMath.LinSpace(Mesh.Endpoint, Mesh.StartPoint, Mesh.numberOfPoints);
            Matrix results = new Matrix();
            for (int k = 0; k < InitialVal.Length + 1; k++)
            {
                Vector res = Vector.Zeros(InitialVal.Length);
                results.Add(res);
            }
            double h = meshSize[1] - meshSize[0];
            int i = 0;
            double x = Mesh.StartPoint;
            bool converge = true;
            Vector Prevy = null!;
            do
            {
                x = x + i * h;
                if (i == 0)
                {
                    results[0].Add(x);
                    for (int j = 1; j < results.Ncol; j++)
                    {
                        results[j].Add(InitialVal[j-1]);
                    }
                    i++;
                    continue;
                }
                else
                {
                    Prevy = ~results[i - 1, 1, 0, "row"];

                    Vector k1 =  model(x, Prevy);
                    Vector k2 =  model(x + h / 4.0, Prevy + k1*(h / 4.0));
                    Vector k3 =  model(x + (3*h / 8.0), Prevy + k1*(3.0*h / 32.0)+k2*(9.0*h/32.0));
                    Vector k4 =  model(x + (12*h/13.0), Prevy + k1*(1932.0*h/2197.0)-k2*(7200.0*h/2197.0)+k3*(7296.0*h/2197.0));
                    Vector k5 = model(x + h, Prevy + k1 * (439.0 * h / 216.0) - 8.0 * h * k2 + k3 * (3680.0 * h / 513.0) - k4 * (845.0 * h / 4104));
                    Vector k6 = model(x + h / 2.0, Prevy - k1 * (8.0 * h / 27.0) + 2.0 * h * k2 - k3 * (3544.0 * h / 2565.0) + k4 * (1859.0 * h / 4104.0) - k5 * (11.0 * h / 40.0));
                    Vector err = h * ((k1 / 360.0) + (128 * k3 / 4275.0) - (2197.0 * k4 / 75240.0) + (k5 / 50.0) + (2.0 * k6 / 55.0));
                    double errNorm = err.Norm_Infinity();
                    double hnext = 0.84 * Math.Pow((tol * h / errNorm), 1.0 / 4.0);
                    converge = errNorm < tol;

                }
            } while (converge);
            return results;
        }
    }
}
