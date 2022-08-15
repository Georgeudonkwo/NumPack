using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using static System.Math;

namespace NumSharp
{
    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class Vector : IEnumerable<double>, ICloneable, IEquatable<Vector>
    {

        #region constructors
        //public Vector(double[]vec)
        //{
        //    this.ColVector = new double[vec.Length, 1];
        //    for (int i = 0; i < vec.Count(); i++)
        //    {
        //        this.ColVector[i, 0] = vec[i];
        //    }
        //    Length = this.ColVector.GetLength(0);
        //    IsColvector = true;
        //    this.RowVector = null!;
        //}
        /// <summary>
        /// 
        /// </summary>
        public Vector(params double[] data)
        {
            this.ColVector = new double[data.Length, 1];
            for (int i = 0; i < data.Count(); i++)
            {
                this.ColVector[i, 0] = data[i];
            }
            Length = this.ColVector.GetLength(0);
            IsColvector = true;
            this.RowVector = null!;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        /// <param name="fun"></param>
        public Vector(IEnumerable<double> vec, Func<double, double> fun = null!)
        {
            List<double> arr = new List<double>();
            if (fun != null)
            {
                arr = vec.Select(fun).ToList();
            }
            else
            {
                arr = vec.Select(v => v).ToList();
            }

            this.ColVector = new double[vec.Count(), 1];
            for (int i = 0; i < vec.Count(); i++)
            {
                this.ColVector[i, 0] = arr[i];
            }
            Length = this.ColVector.GetLength(0);
            IsColvector = true;
            this.RowVector = null!;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public Vector(List<double> vec)
        {
            this.ColVector = new double[vec.Count, 1];
            for (int i = 0; i < vec.Count(); i++)
            {
                this.ColVector[i, 0] = vec[i];
            }
            Length = this.ColVector.GetLength(0);
            IsColvector = true;
            this.RowVector = null!;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public Vector(double vec)
        {
            this.ColVector = new double[1, 1];
            this.ColVector[0, 0] = vec;
            Length = this.ColVector.GetLength(0);
            IsColvector = true;
            this.RowVector = null!;
        }
        /// <summary>
        /// Create an empty vector
        /// You can use the <see cref="Add" />Method to add elements
        /// </summary>
        public Vector()
        {
            this.colVec = new double[0, 0];
            Length = -1;
            IsColvector = true;
            this.RowVector = null!;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim"></param>
        public Vector(int dim)
        {
            this.colVec = new double[dim, 1];
            Length = dim;
            IsColvector = true;
            this.RowVector = null!;
        }
        //public Vector (ColVec Cv)
        //{
        //    Length = Cv.Count();
        //    this.ColVector = new double[Length, 1];
        //    for (int i = 0; i < Cv.Count(); i++)
        //    {
        //        this.ColVector[i, 0] = Cv[i];
        //    }
        //    this.IsColvector = true;
        //    this.RowVector = null;
        //}
        //public Vector(RowVec Rv)
        //{
        //    Length = Rv.Count();
        //    this.rowvector = new double[1, Length];
        //    for (int i = 0; i < Rv.Count(); i++)
        //    {
        //        this.rowvector[0, i] = Rv[i];
        //    }
        //    this.IsColvector = false;
        //    this.ColVector = null;
        //}
        #endregion
        #region operators

        //public static implicit operator ColVec(Vector V)
        //{
        //    double[] arr = new double[V.Length];
        //    for (int i = 0; i < arr.Length; i++)
        //    {
        //        arr[i] = V[i];
        //    }
        //    return new ColVec(arr);
        //}
        //public static implicit operator Vector(ColVec V)
        //{
        //    return new Vector(V);
        //}

        //public static implicit operator RowVec(Vector V)
        //{
        //    double[] arr = new double[V.Length];
        //    for (int i = 0; i < arr.Length; i++)
        //    {
        //        arr[i] = V[i];
        //    }
        //    return new RowVec(arr);
        //}
        //public static implicit operator Vector(RowVec V)
        //{
        //    return new Vector(V);
        //}
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>

        public static explicit operator double[](Vector vec)
        {
            double[] res = new double[vec.Length];
            for (int i = 0; i < res.Count(); i++)
            {
                res[i] = vec[i];
            }
            return res;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="scaler"></param>
        public static implicit operator Vector(double scaler)
        {
            Vector res = new Vector(scaler);
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static explicit operator double(Vector vec)
        {
            double res = vec[0];
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static implicit operator List<double>(Vector vec)
        {
            List<double> res = new List<double>();
            for (int i = 0; i < vec.Length; i++)
            {
                res.Add(vec[i]);
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static implicit operator Vector(List<double> vec)
        {
            Vector res = new Vector();
            res.ColVector = new double[vec.Count, 1];
            for (int i = 0; i < res.ColVector.GetLength(0); i++)
            {
                res.ColVector[i, 0] = vec[i];
            }
            res.Length = vec.Count;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static implicit operator Vector(double[] vec)
        {
            Vector res = new Vector();
            res.ColVector = new double[vec.Length, 1];
            for (int i = 0; i < res.ColVector.GetLength(0); i++)
            {
                res.ColVector[i, 0] = vec[i];
            }
            res.Length = vec.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static explicit operator Vector(double[,] vec)
        {
            if (vec.GetLength(1) > 1) { throw new Exception("Matrix must be a column vector"); }
            Vector res = new Vector();
            res.ColVector = new double[vec.Length, 1];
            for (int i = 0; i < res.ColVector.GetLength(0); i++)
            {
                res.ColVector[i, 0] = vec[i, 0];
            }
            res.Length = vec.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static implicit operator double[,](Vector vec)
        {
            double[,] res = new double[vec.Length, 1];
            for (int i = 0; i < vec.ColVector.GetLength(0); i++)
            {
                res[i, 0] = vec[i];
            }
            return res;
        }
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool operator ==(Vector v1, Vector v2)
        {
            bool isequal = false;
            isequal = v1.Equals(v2);
            return isequal;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool operator !=(Vector v1, Vector v2)
        {

            bool isNotequal = v1 == v2;// true;
            if (isNotequal == true) { return false; } else { return true; }
            #region
            //if (ReferenceEquals(v2, null))
            //{
            //    if (ReferenceEquals(v1, null)) { return false; }
            //    else { return true; }
            //}

            //    if (v1.Length == v2.Length) { isNotequal = false; }
            //    for (int i = 0; i < v1.Length; i++)
            //    {
            //        bool isValueEqual = v1[i] == v2[i];
            //        bool isSignEqual = Sign(v1[i]) == Sign(v2[i]);
            //        if (isValueEqual && isSignEqual)
            //        { isNotequal = false; }
            //        else { isNotequal = true; break; }
            //    }
            //return isNotequal;
            #endregion
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static Vector operator +(Vector v1, Vector V2)
        {
            return Addition(v1, V2);
        }
        /// <summary>
        /// transpose row vector to column vector and vice versa
        /// </summary>
        /// <param name="V"></param>
        /// <returns></returns>
        public static Vector operator ~(Vector V)
        {
            if (V.RowVector == null)
            {
                return V.GetRowVector();
            }
            else
            {
                double[] C = new double[V.rowvector.GetLength(1)];
                for (int i = 0; i < V.RowVector.GetLength(1); i++)
                {
                    C[i] = V.RowVector[0, i];
                }
                return new Vector(C);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector operator +(Vector v1, double Scalar)
        {
            return Addition(v1, Scalar);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="v1"></param>
        /// <returns></returns>
        public static Vector operator +(double Scalar, Vector v1)
        {
            return Addition(Scalar, v1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static Vector operator -(Vector v1, Vector V2)
        {
            return Subtraction(v1, V2);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector operator -(Vector v1, double Scalar)
        {
            return Subtraction(v1, Scalar);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="v1"></param>
        /// <returns></returns>
        public static Vector operator -(double Scalar, Vector v1)
        {
            return Subtraction(Scalar, v1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <returns></returns>
        public static Vector operator -(Vector v1)
        {
            return Multiplication(-1, v1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static Vector operator *(Vector v1, Vector v2)
        {
            return Multiplication(v1, v2);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector operator *(Vector v1, double Scalar)
        {
            return Multiplication(v1, Scalar);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="v1"></param>
        /// <returns></returns>
        public static Vector operator *(double Scalar, Vector v1)
        {
            return Multiplication(Scalar, v1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static Vector operator /(Vector v1, Vector V2)
        {
            return TwDivision(v1, V2);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector operator /(Vector v1, double Scalar)
        {
            return TwDivision(v1, Scalar);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="v1"></param>
        /// <returns></returns>
        public static Vector operator /(double Scalar, Vector v1)
        {
            return TwDivision(Scalar, v1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <param name="exponent"></param>
        /// <returns></returns>
        public static Vector operator ^(Vector V, double exponent)
        {
            return Pow(V, exponent);
        }


        #endregion
        #region Indexers
        /// <summary>
        /// 
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        public double this[int i]
        {
            get
            {
                if (this.ColVector != null) { return this.ColVector[i, 0]; }
                else { return this.RowVector[0, i]; }
            }
            set
            {
                if (this.ColVector != null) { this.ColVector[i, 0] = value; }
                else { this.RowVector[0, i] = value; }
            }

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="i"></param>
        /// <param name="vectorType"></param>
        /// <returns></returns>
        public Vector this[int i, string vectorType = "col"]
        {
            get
            {
                if (vectorType.Equals("col", StringComparison.InvariantCultureIgnoreCase))
                {
                    List<double> v = new List<double>();
                    for (int k = i; k < this.colVec.GetLength(0); k++)
                    {
                        v.Add(this.colVec[k, 0]);
                    }
                    return new Vector(v);
                }
                else
                {
                    List<double> v = new List<double>();
                    //for (int k = i; k < this.rowvector.GetLength(0); k++)
                    for (int k = i; k < this.rowvector.GetLength(1); k++)
                    {
                        v.Add(this.rowvector[0, k]);
                    }

                    return new Vector(v).GetRowVector();
                }
            }
            set
            {
                if (vectorType.Equals("col", StringComparison.InvariantCultureIgnoreCase))
                {
                    for (int k = 0; k < value.colVec.Length; k++)
                    {
                        this.colVec[k + i, 0] = value[k];
                    }
                }
                else
                {
                    for (int k = 0; k < value.rowvector.Length; k++)
                    {
                        this.colVec[0, k + i] = value[k];
                    }
                }
            }
        }

        #endregion
        #region Methods
        /// <summary>
        /// calculate the mean value of the vector components
        /// </summary>
        /// <returns></returns>
        public double Mean()
        {
            return this.Average();
        }
        /// <summary>
        /// calculate the variance of the vector components
        /// </summary>
        /// <returns></returns>
        public double Variance()
        {
            double mean = Mean();
            Vector v = this - mean;
            v = v ^ 2;
            double res = v.Sum() / (Length - 1);
            return res;
        }
        /// <summary>
        /// calculate the standard deviation of the vector components
        /// </summary>
        /// <returns></returns>
        public double STD()
        {
            double v = Variance();
            return Math.Sqrt(v);
        }
        /// <summary>
        /// normalizes the vector,such that the vector length is unity
        /// </summary>
        /// <returns></returns>
        public Vector Normalized()
        {
            Vector v = this;
            v = this / Norm();
            return v;
        }
        /// <summary>
        /// concatenate the vectors,essentially flattening them into a single vector
        /// </summary>
        /// <param name="vec"></param>
        /// <returns></returns>
        public static Vector Concat(params Vector[] vec)
        {
            Vector res = new Vector();
            for (int k = 0; k < vec.Length; k++)
            {
                for (int i = 0; i < vec[k].Length; i++)
                {
                    res.Add(vec[k][i]);
                }
            }
            return res;

        }
        /// <summary>
        /// Insert the value at the top of the vector, at index 0
        /// </summary>
        /// <param name="val"></param>
        /// <returns></returns>
        public Vector Prepend(double val)
        {
            return Insert(val, 0);
        }
        /// <summary>
        /// Add item to the end of the vector <see cref="Vector.Add(double)"/>
        /// </summary>
        /// <param name="val"></param>
        /// <returns></returns>
        public Vector Append(double val)
        {
            return Insert(val, this.Length);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="fun"></param>
        /// <returns></returns>
        public T Get<T>(Func<T> fun)
        {
            return fun();
        }
        /// <summary>
        /// 
        /// </summary>
        public double Transform(Func<double, double> fun)
        {

            return fun(this[0]);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="v"></param>
        /// <param name="fun"></param>
        /// <returns></returns>
        public static T Get<T>(Vector v, Func<Vector, T> fun)
        {
            return fun(v);
        }
        /// <summary>
        /// return true if any component of the vector contains NAN of infinity
        /// </summary>
        /// <returns></returns>
        public bool NotAVector()
        {
            if (this.Contains(double.NaN) || this.Contains(double.PositiveInfinity) || this.Contains(double.NegativeInfinity)) return true;
            return false;
        }
        /// <summary>
        /// return the index of the vector component or -1 if the component is not in the vector
        /// </summary>
        /// <param name="val"></param>
        /// <returns></returns>
        public int IndexOf(double val)
        {
            for (int i = 0; i < this.Length; i++)
            {
                if (val == this[i]) return i;
            }
            return -1;
        }
        /// <summary>
        /// return the row vector
        /// </summary>
        /// <returns></returns>
        public Vector GetRowVector()
        {

            if (this.RowVector != null) { return this.Copy(); }
            Vector? Vcopy = this!.Copy()as Vector;
            Vector res = new Vector();
            res.ColVector = null!;
            res.RowVector = new double[1, Vcopy!.ColVector.GetLength(0)];
            for (int i = 0; i < Vcopy.Length; i++)
            {
                res.RowVector[0, i] = Vcopy[i];
            }
            res.Length = Vcopy.Length;
            res.IsColvector = false; ;
            return res;
        }
        /// <summary>
        /// return the row vector of V
        /// </summary>
        /// <param name="V"></param>
        /// <returns></returns>
        public static Vector GetRowVector(Vector V)
        {
            if (V.RowVector != null) { return V.Copy(); ; }
            Vector res = new Vector();
            res.ColVector = null!;
            res.RowVector = new double[1, V.Length];
            for (int i = 0; i < V.Length; i++)
            {

                res.RowVector[0, i] = V[i];
            }
            res.Length = V.Length;
           res. IsColvector = false; ;
            return res;
        }
        /// <summary>
        /// return a column vector
        /// </summary>
        /// <param name="vect"></param>
        public  void CopyTo(Vector vect)
        {
             vect = new Vector();
            vect.ColVector = new double[vect.Length, 1];
            vect.Length = vect.Length;
            for (int i = 0; i < this. Length; i++)
            {
                vect.ColVector[i, 0] = this.ColVector[i, 0];
            }
            vect.Length = this.Length;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Vector Copy()
        {
            Vector res = new Vector();
            if (this.ColVector != null)
            {
                res.ColVector = new double[this.Length, 1];
                res.Length = this.Length;
                res.IsColvector = true;
                res.RowVector = null!;
                for (int i = 0; i < this.ColVector.GetLength(0); i++)
                {
                    res.ColVector[i, 0] = this.ColVector[i, 0];
                }
               
            }
            else
            {
                res.RowVector = new double[1, this.Length];
                res.Length = this.Length;
                for (int i = 0; i < this.RowVector.GetLength(1); i++)
                {
                    res.RowVector[0, i] = this.RowVector[0, i];
                }
                res.ColVector = null!;
                res.IsColvector = false;
                res.Length = this.Length;
            }
            return res;
        }
        /// <summary>
        /// Calculates the inner product of two vectors
        /// </summary>
        /// <param name="V1"> A ROW VECTOR</param>
        /// <param name="V2"> A Column vector</param>
        /// <returns> returns a double</returns>
        public static double Dot(Vector V1,Vector V2)
        {
            double res = 0.0;
            for (int i = 0; i < V1.RowVector.GetLength(0); i++)
            {
                for (int j = 0; j < V2.ColVector.GetLength(0); j++)
                {
                    res +=V1.RowVector[0, j] * V2.ColVector[j, 0];
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <returns></returns>
        public  double Dot(Vector V)
        {

            double res = 0.0;
            Vector RV = this.GetRowVector();
            for (int i = 0; i < RV.rowvector.GetLength(0); i++)
            {
                for (int j = 0; j < V.ColVector.GetLength(0); j++)
                {
                    res += RV.rowvector[0, j] * V.ColVector[j, 0];
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public  double Max()
        {
            Vector Cv = this.Copy();
            if (Cv.ColVector == null)
            {
                Cv = ~Cv;
                return Max();
            }
            bool exit = true;
            int iter = 0;
            do
            {
                for (int i = 0; i < Cv.Length - 1; i++)
                {
                    double T = Cv.ColVector[i, 0];
                    double Tnext = Cv.ColVector[i + 1, 0];
                    if (T < Tnext)
                    {
                        Swap(ref Cv.ColVector[i, 0], ref Cv.ColVector[i + 1, 0]);
                    }
                   
                }
                exit = iter < Cv.Length - 1;//  Cv.vec[0, 0] < Cv.vec[ 1, 0];
                iter++;

            } while (exit);
          
          return Cv.ColVector[0,0];
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Min()
        {
           
            Vector Cv = this.Copy();
            if (Cv.ColVector == null)
            {
                Cv = ~Cv;
                return Min();
            }
            bool exit = true;
            int iter = 0;
            do
            {
                for (int i = 0; i < Cv.Length - 1; i++)
                {
                    double T = Cv.ColVector[i, 0];
                    double Tnext = Cv.ColVector[i + 1, 0];
                    if (T > Tnext)
                    {
                        Swap(ref Cv.ColVector[i, 0], ref Cv.ColVector[i + 1, 0]);
                    }

                }
               
                exit =iter<Cv.Length-1;
                iter++;
            } while (exit);

            return Cv.ColVector[0, 0];
        }
        /// <summary>
        /// exchange two items
        /// </summary>
        /// <param name="Item1"></param>
        /// <param name="Item2"></param>
        public void Swap(ref double Item1,ref double Item2)
        {
            
            double temp = Item1;
            Item1 = Item2;
            Item2 = temp;
        }
        /// <summary>
        /// perform vector additrion
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns> returns a column vector</returns>
        /// <exception cref="Exception"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Vector Addition(Vector V1,Vector V2)
        {
            if (V1.ColVector == null && V2.ColVector == null)
            {
                return ~Addition(~V1, ~V2);
            }
            if ((V1.ColVector == null && V1.ColVector != null) || (V1.RowVector == null && V1.rowvector != null)) { throw new Exception("Both vectors must be either column or row vectors"); }
            if (V1.Length != V2.Length) { throw new InvalidOperationException($"{V1.Length} is not the same as {V2.Length}"); }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i,0] = V1.ColVector![i, 0] + V2.ColVector[i, 0];
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector Addition(Vector V1, double Scalar)
        {
            if (V1.ColVector == null)
            {
                return ~Addition(~V1, Scalar);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] + Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="V1"></param>
        /// <returns></returns>
        public static Vector Addition(double Scalar, Vector V1)
        {
            if (V1.ColVector == null)
            {
                return ~Addition(Scalar,~V1);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] + Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// subtract two vectors
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Vector Subtraction(Vector V1, Vector V2)
        {
            if (V1.ColVector == null && V2.ColVector == null)
            {
                return ~Subtraction(~V1, ~V2);
            }
            if ((V1.ColVector == null && V1.ColVector != null) || (V1.RowVector == null && V1.rowvector != null)) { throw new Exception("Both vectors must be either column or row vectors"); }
            if (V1.Length != V2.Length) { throw new InvalidOperationException($"{V1.Length} is not the same as {V2.Length}"); }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector![i, 0] - V2.ColVector[i, 0];
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector Subtraction(Vector V1,double Scalar)
        {
            if (V1.ColVector == null)
            {
                return ~Subtraction(~V1, Scalar);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] - Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="V1"></param>
        /// <returns></returns>
        public static Vector Subtraction(double Scalar,Vector V1)
        {
            if (V1.ColVector == null)
            {
                return ~Subtraction( Scalar, ~V1);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = Scalar- V1.ColVector[i, 0]  ;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Vector Multiplication(Vector V1, Vector V2)
        {
            //if (V1.ColVector == null && V2.ColVector == null)
            //{
            //    return ~Multiplication(~V1, ~V2);
            //}
           
            if ((V1.ColVector == null && V1.ColVector != null) || (V1.RowVector == null && V1.rowvector != null)) { throw new Exception("Both vectors must be either column or row vectors"); }
            if (V1.Length != V2.Length) { throw new InvalidOperationException($"{V1.Length} is not the same as {V2.Length}"); }
            Vector res = new Vector();
            if (V1.isColvector&&V2.isColvector)
            {
                res.ColVector = new double[V1.Length, 1];
                for (int i = 0; i < V1.Length; i++)
                {
                    res.ColVector[i, 0] = V1.ColVector![i, 0] * V2.ColVector[i, 0];
                }
                res.Length = V1.Length;
            }
            else if (!V1.isColvector && V2.isColvector)
            {
                res.ColVector = new double[1, 1];
                for (int i = 0; i < V1.Length; i++)
                {
                    res.ColVector[0, 0] += V1.RowVector![0, i] * V2.ColVector[i, 0];
                }

                res.Length = 1;// V1.Length;
            }
            else if (!V1.isColvector && !V2.isColvector)
            {
                res.RowVector = new double[1, V1.Length];
                for (int i = 0; i < V1.Length; i++)
                {
                    res.RowVector[0, i] = V1.RowVector![0, i] * V2.RowVector[0, i];
                }
                res.Length = V1.Length;
            }
            else if (V1.isColvector && !V2.isColvector)
            {
                throw new Exception("This operation will return a matrix, please use the OuterProduct method in the MATRIX class");
            }
            else { throw new Exception("incompactible vectors"); }
                return res;
        }
        /// <summary>
        /// scalar vector multiplication
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector Multiplication(Vector V1, double  Scalar)
        {
            if (V1.ColVector == null)
            {
                return ~Multiplication(~V1, Scalar);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] * Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="V1"></param>
        /// <returns></returns>
        public static Vector Multiplication(double Scalar,Vector V1 )
        {
            if (V1.ColVector == null)
            {
                return ~Multiplication(Scalar, ~V1);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] * Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// termwise division
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Vector TwDivision(Vector V1, Vector V2)
        {
            if (V1.ColVector == null && V2.ColVector == null)
            {
                return ~TwDivision(~V1, ~V2);
            }
            if((V1.ColVector==null&&V1.ColVector!=null)|| (V1.RowVector == null && V1.rowvector!= null)){ throw new Exception("Both vectors must be either column or row vectors"); }
            if (V1.Length != V2.Length) { throw new InvalidOperationException($"{V1.Length} is not the same as {V2.Length}"); }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector![i, 0] / V2.ColVector[i, 0];
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="Scalar"></param>
        /// <returns></returns>
        public static Vector TwDivision(Vector V1, double Scalar)
        {
            if (V1.ColVector == null)
            {
                return ~TwDivision( ~V1, Scalar);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = V1.ColVector[i, 0] / Scalar;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scalar"></param>
        /// <param name="V1"></param>
        /// <returns></returns>
        public static Vector TwDivision(double Scalar, Vector V1)
        {
            if (V1.ColVector == null)
            {
                return ~TwDivision(Scalar, ~V1);
            }
            Vector res = new Vector();
            res.ColVector = new double[V1.Length, 1];
            for (int i = 0; i < V1.Length; i++)
            {
                res.ColVector[i, 0] = Scalar/ V1.ColVector[i, 0] ;
            }
            res.Length = V1.Length;
            return res;
        }
        /// <summary>
        /// return the first nth entries
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public Vector FirstNth(int n)
        {
            Vector res = new Vector();
            for (int i = 0; i < n; i++)
            {
                res.Add(this[i]);
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="val"></param>
        /// <param name="index"></param>
        public  void Replace(double val,int index)
        {
            for (int i = 0; i < this.Length; i++)
            {
                if (i == index) { this[i] = val; }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="val"></param>
        /// <param name="newval"></param>
        public void Replace(double val ,double newval)
        {
            for (int i = 0; i < this.Length; i++)
            {
                if (this[i] == val) { this[i] = newval; }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public Vector LastNth(int n)
        {
            Vector res = new Vector();
            for (int i = n; i < this.Length; i++)
            {
                res.Add(this[i]);
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <returns></returns>
        public Vector ggg(Func<double,double> fun)
        {
            var res = new Vector();
            for (int i = 0; i < this.Count(); i++)
            {
                res.Add(fun(this[i]));
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="exponent"></param>
        /// <returns></returns>
        public Vector Pow(double exponent)
        {
            Vector res = this.Copy();
            if (res.isColvector)
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.ColVector[i, 0] = Math.Pow(res.ColVector[i, 0], exponent);
                }
            }
            else
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.RowVector[0,i] = Math.Pow(res.RowVector[0, i], exponent);
                }
            }
            
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <param name="exponent"></param>
        /// <returns></returns>
        public static Vector Pow(Vector V, double exponent)
        {
            Vector res = V.Copy();
            if (res.IsColvector)
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.ColVector[i, 0] = Math.Pow(res.ColVector[i, 0], exponent);
                }
            }
            else
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.RowVector[0, i] = Math.Pow(res.RowVector[0, i], exponent);
                }
            }
           
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim"></param>
        /// <returns></returns>
        public static Vector Ones(int dim)
        {
            Vector one = new Vector(dim);
            for (int i = 0; i < one.ColVector.Length; i++)
            {
                one[i] = 1.0;
            }
            return one;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim"></param>
        /// <returns></returns>
        public static Vector Zeros(int dim)
        {
            Vector one = new Vector(dim);
            for (int i = 0; i < one.ColVector.Length; i++)
            {
                one[i] = 0.0;
            }
            return one;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Norm()
        {
            double N = 0.0;
            Vector V = this.Copy();
            double maxVal = V.Max<double>(d=>Abs(d)); // Max();
          V = V / maxVal;
            V = V.Absolute().Pow(2);
            double S = V.Sum();
            N = maxVal * Math.Pow(S,1.0/2.0);
            return N;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Vector Absolute()
        {
            Vector res = this.Copy();
            if (res.ColVector != null)
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.ColVector[i, 0] = Math.Abs(res.ColVector[i, 0]);
                }
            }
            else
            {
                for (int i = 0; i < res.Length; i++)
                {
                    res.RowVector[0, i] = Math.Abs(res.RowVector[0, i]);
                }
            }
            
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Norm_1()
        {
            double res = 0.0;
            if (this.isColvector)
            {
                for (int i = 0; i < this.Length; i++)
                {
                    res += Abs(this.ColVector[i, 0]);
                }
            }
            else
            {
                for (int i = 0; i < this.Length; i++)
                {
                    res += Abs(this.RowVector[0, i]);
                }
            }
            
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="infinity"></param>
        /// <returns></returns>
        public double Norm_Infinity(int infinity)
        {
            double N = 0.0;
           
            Vector V = this.Copy();
            double maxVal = V.Max<double>();
            V = V / maxVal;
            V = V.Absolute().Pow(infinity);
            double S = V.Sum();
            N = maxVal * Math.Pow(S, 1.0 /infinity);
            return N;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Norm_Infinity()
        {
          
            Vector V = this.Copy();
            double N = Abs(V[0]);
            for (int i = 1; i < V.Length; i++)
            {
                if (Abs( V[i]) > N) { N =Abs(V[i]); }
            }
            return N;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Val"></param>
        /// <param name="position"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        public Vector Insert(double Val,int position)
        {
            Vector cVec =(Vector) this.Clone();
            double[,]cmat = null!;
           
            if (IsColvector)
            {
                cmat = new double[cVec.Length + 1, 1];
                if (position > cmat.GetLength(0)) { throw new InvalidOperationException($"invalid index {position} can not be greater than{cmat.GetLength(0)}"); }
                cmat[position, 0] = Val;
                bool met = false;
                for (int i = 0; i < cVec.Length; i++)
                {
                    if (i == position)
                    {
                        cmat[i + 1, 0] = cVec.colVec[i, 0];
                        met = true;
                    }
                    else
                    {
                        if (met)
                        {
                            cmat[i + 1, 0] = cVec.colVec[i, 0];
                        }
                        else
                        {
                            cmat[i, 0] = cVec.colVec[i, 0];
                        }

                    }

                }
            }
            else
            {
                cmat = new double[ 1, cVec.Length + 1];
                if (position > cmat.GetLength(1)) { throw new InvalidOperationException($"invalid index {position} can not be greater than{cmat.GetLength(0)}"); }
                cmat[ 0, position] = Val;
                bool met = false;
                for (int i = 0; i < cVec.Length; i++)
                {
                    if (i == position)
                    {
                        cmat[0, i + 1] = cVec.RowVector[0, i];
                        met = true;
                    }
                    else
                    {
                        if (met)
                        {
                            cmat[ 0, i + 1] = cVec.RowVector[0, i];
                        }
                        else
                        {
                            cmat[0, i ] = cVec.RowVector[0, i];
                        }

                    }

                }
            }
            
            Vector res = this;
            res.Length = cmat.GetLength(0);
            res.colVec = cmat;
            //var Ro = res.GetRowVector();
            //res.RowVector = Ro.RowVector;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="position"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public Vector Remove( int position)
        {
            if (position < 0) { throw new Exception($"{position} is not a valid index"); }
            Vector cVec = (Vector)this.Clone();
            double[,] cmat = null!;
            if (IsColvector)
            {
                 cmat = new double[cVec.Length - 1, 1];
                if (position > cmat.GetLength(0)) { throw new InvalidOperationException($"invalid index {position} can not be greater than{cmat.GetLength(0)}"); }
                bool met = false;
                for (int i = 0; i < cVec.Length - 1; i++)
                {
                    if (i == position)
                    {
                        cmat[i, 0] = cVec.colVec[i + 1, 0];
                        met = true;
                    }
                    else
                    {
                        if (met)
                        {
                            cmat[i, 0] = cVec.colVec[i + 1, 0];
                        }
                        else
                        {
                            cmat[i, 0] = cVec.colVec[i, 0];
                        }

                    }

                }
            }
            else
            {
                 cmat = new double[ 1, cVec.Length - 1];
                if (position > cmat.GetLength(1)) { throw new InvalidOperationException($"invalid index {position} can not be greater than{cmat.GetLength(0)}"); }
                bool met = false;
                for (int i = 0; i < cVec.Length - 1; i++)
                {
                    if (i == position)
                    {
                        cmat[0, i] = cVec.RowVector[0, i + 1];
                        met = true;
                    }
                    else
                    {
                        if (met)
                        {
                            cmat[0, i] = cVec.RowVector[0, i + i];
                        }
                        else
                        {
                            cmat[0, i] = cVec.RowVector[0, i];
                        }
                    }
                }
            }
           
            Vector res = this;
            res.Length = res.isColvector ? cmat.GetLength(0) : cmat.GetLength(1);
            if(isColvector)
            {
                res.ColVector=cmat;
            }
            else
            {
                res.RowVector= cmat;
            }
            
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public Vector Remove(double value)
        {
            if (!this.Contains(value)) { throw new Exception($"{value} is not an item of the {this} Vector"); }
            int index = this.IndexOf(value);
            return Remove(index);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <returns></returns>
        public double CosTheta(Vector V)
        {
            return this.Dot(V) / (this.Norm() * V.Norm());
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static double CosTheta(Vector V1, Vector V2)
        {
            return V1.Dot(V2) / (V1.Norm() * V2.Norm());
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="unit"></param>
        /// <returns></returns>
        public static double Angle(Vector v1, Vector v2, string unit = "deg")
        {
            if (unit.Equals("deg", StringComparison.InvariantCultureIgnoreCase))
            {
                double theta = CosTheta(v1, v2);
                if (theta > 1.0) { theta = 1.0; }
                if (theta < -1) { theta = -1; }
                return Acos(theta) * (180 / PI);
            }
            else
            {
                double theta = CosTheta(v1, v2);
                if (theta > 1.0) { theta = 1.0; }
                if (theta < -1) { theta = -1; }
                return Acos(theta);
            }

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <param name="unit"></param>
        /// <returns></returns>
        public double Angle(Vector v,string unit="deg")
        {
            if (unit.Equals("deg", StringComparison.InvariantCultureIgnoreCase))
            {
                double theta = CosTheta(v);
                if (theta > 1.0) { theta = 1.0; }
                if (theta < -1) { theta = -1; }
                return Acos(theta) * (180/PI);
            }
            else
            {
                double theta = CosTheta(v);
                if (theta > 1.0) { theta = 1.0; }
                if (theta <-1) { theta = -1; }
                return Acos(theta);
            }
            
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Vector Reverse()
        {
            Vector Vc = this.Copy();
            List<double> Rvc = Vc.Reverse<double>().ToList();
            for (int i = 0; i < Rvc.Count; i++)
            {
                Vc[i] = Rvc[i];
            }
            return Vc;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Vector cumSum()
        {
            Vector  vec = new Vector() { this[0]};
            for (int i = 1; i < this.Length; i++)
            {
                vec.Add(vec[i-1]+this[i]);
            }
            return vec;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Vector CumProd()
        {
            Vector prod = new Vector() { this[0] };
            for (int i = 1; i < this.Length; i++)
            {
                prod.Add(prod[i-1]*this[i]);
            }
            return prod;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Val"></param>
        /// <returns></returns>
        public Vector Add(double Val)
        {
            Vector cVec = (Vector)this.Clone();
            if (cVec.IsColvector)
            {
                if (cVec.Length == -1) { cVec.Length = 0; }
                double[,] cmat = new double[cVec.Length + 1, 1];

                cmat[cVec.Length, 0] = Val;
                for (int i = 0; i < cVec.Length; i++)
                {
                    cmat[i, 0] = cVec.ColVector[i, 0];
                }
                Vector res = this;
                res.Length = cmat.GetLength(0);
                res.colVec = cmat;
                res.RowVector = null!;
                return res;
            }
            else
            {
                if (cVec.Length == -1) { cVec.Length = 0; }
                double[,] cmat = new double[ 1, cVec.Length + 1];

                cmat[0, cVec.Length] = Val;
                for (int i = 0; i < cVec.Length; i++)
                {
                    cmat[0, i] = cVec.RowVector[0, i];
                }
                Vector res = this;
                res.Length = cmat.GetLength(1);
                res.colVec = null!;
                res.RowVector = cmat;
                return res;
            }
            
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fun"></param>
        /// <returns></returns>
        public  Vector Transform(  Func<Vector, Vector> fun)
        {
            return fun(this);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Log(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Log(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Log10(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Log10(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Sin(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Sin(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Cos(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Cos(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Exp(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Exp(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Sqrt(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Sqrt(c[i]);
            }
            return c;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static Vector Tan(Vector v)
        {
            Vector c = v.Copy();
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Math.Tan(c[i]);
            }
            return c;
        }
        
        #endregion
        #region implemented interfaces
        /// <summary>
        /// An enumerator for iterating through the vector
        /// </summary>
        /// <returns></returns>
        public IEnumerator<double> GetEnumerator()
        {
            if (ColVector?.Length > 0)
            {
                for (int i = 0; i < ColVector.GetLength(0); i++)
                {
                    yield return ColVector[i, 0];
                }
            }
            else
            {
                for (int i = 0; i < RowVector.GetLength(1); i++)
                {
                    yield return RowVector[0, i];
                }
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            if (ColVector?.Length>0)
            {
                for (int i = 0; i < ColVector.GetLength(0); i++)
                {
                    yield return ColVector[i, 0];
                }
            }
            else
            {
                for (int i = 0; i < RowVector.GetLength(1); i++)
                {
                    yield return RowVector[0, 1];
                }
            }
           
        }
        /// <summary>
        /// shallow copy of the vector object
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            //using (MemoryStream str = new MemoryStream())
            //{
                //BinaryFormatter formater = new BinaryFormatter();
                //formater.Serialize(str, this);
                //str.Position = 0;
                //Vector Cvector = (Vector)formater.Deserialize(str);
                
                Vector Cvector = new Vector();
                for (int i = 0; i < this.Length; i++)
                {
                    Cvector.Add(this[i]);
                }
                return Cvector;
            //}
        }
        /// <summary>
        /// equality check 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Vector? other)
        {
            return this.SequenceEqual(other!);
            //bool isequal = false;
            //if (ReferenceEquals(this, other)) { return true; }
            //if (this.All((v) => double.IsNaN(v) && other!.All((c) => double.IsNaN(c)))) return true;
            //if (other!.Contains(double.NaN) || this.Contains(double.NaN)) { return false; }
            //if (other!.Any((x)=>double.IsInfinity(x)) || this.Any((x) => double.IsInfinity(x))) { return false; }
            //if (this.Length == other!.Length) { isequal = true; }
            //    for (int i = 0; i < this.Length; i++)
            //    {
            //      bool isValueEqual = (Abs(this[i] - other[i])) < 1e-6;
            //      bool isSignEqual = Sign(this[i]) == Sign(other[i]);
            //        if (isValueEqual && isSignEqual)
            //        { isequal = true; }
            //        else { isequal = false; break; }
            //    }
            //    return isequal;
        }
        /// <summary>
        /// implements equality semantic
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object? obj)
        {
            Vector V = (Vector)obj!;
            return Equals(V);
        }
        /// <summary>
        /// Calculate the hash code of the object
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            int res = 0;
            for (int i = 0; i < this.Length; i++)
            {
                res += this[i].GetHashCode() + 17;
            }
            return res;
        }
        #endregion
        double[,] colVec;
        double[,] rowvector;
        bool isColvector = true;
        /// <summary>
        /// vector length, denote number of components in the vector
        /// </summary>
        public int Length { get; private set; }
        /// <summary>
        /// column vector
        /// </summary>
        public double[,] ColVector { get => colVec; private set => colVec = value; }
        /// <summary>
        /// Row vector
        /// </summary>
        public double[,] RowVector { get => rowvector; set => rowvector = value; }
        /// <summary>
        /// Return true if the vector is a row vector
        /// </summary>
        public bool IsColvector { get => isColvector; set => isColvector = value; }
    }
}
