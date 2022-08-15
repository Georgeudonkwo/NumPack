using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Data;
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
    public enum SwapState { 
        /// <summary>
        /// 
        /// </summary>
        RowWise,
        /// <summary>
        /// 
        /// </summary>
        ColumnWise,
        /// <summary>
        /// 
        /// </summary>
        BothRowAndColumnWise}
    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    public class Matrix : IEnumerable<Vector>, ICloneable, IEquatable<Matrix>
    {
        #region Constructors
        /// <summary>
        /// 
        /// </summary>
        /// <param name="array"></param>
        public Matrix(double[,] array)
        {
            mat = new double[array.GetLength(0), array.GetLength(1)];
            Nrow = array.GetLength(0);
            Ncol = array.GetLength(1);
            for (int i = 0; i < Nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    mat[i, j] = array[i, j];

                }
            }

        }
        /// <summary>
        /// copy constructor
        /// </summary>
        /// <param name="A"></param>
        public Matrix(Matrix A)
        {
            mat = new double[A.Nrow, A.Ncol];
            Nrow = A.Nrow;
            Ncol = A.Ncol;
            for (int i = 0; i < Nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    mat[i, j] = A[i, j];
                }
            }

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public Matrix(Vector vec)
        {
            if (vec.ColVector?.GetLength(0) > 1)
            {
                Nrow = vec.ColVector.GetLength(0);
                Ncol = 1;
                mat = new double[Nrow, Ncol];
                for (int i = 0; i < Nrow; i++)
                {
                    mat[i, 0] = vec.ColVector[i, 0];

                }
            }
            else
            {
                Nrow = 1;
                Ncol = vec.RowVector.GetLength(1);
                mat = new double[Nrow, Ncol];
                for (int i = 0; i < Nrow; i++)
                {
                    mat[0, i] = vec.RowVector[0, i];
                }
            }

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ListArray"></param>
        /// <exception cref="InvalidOperationException"></exception>
        public Matrix(List<List<double>> ListArray)
        {
            List<int> NumListItems = new List<int>();
            NumListItems.Add(ListArray[0].Count);
            Ncol = ListArray.Count;
            Nrow = ListArray[0].Count;
            mat = new double[Nrow, Ncol];
            for (int i = 0; i < Ncol; i++)
            {
                NumListItems.Add(ListArray[i].Count);
                if (NumListItems[0] != NumListItems[i]) { throw new InvalidOperationException(); }
                for (int j = 0; j < Nrow; j++)
                {
                    mat[j, i] = ListArray[i][j];
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ListArray"></param>
        /// <param name="Orientation"></param>
        /// <exception cref="InvalidOperationException"></exception>
        public Matrix(List<List<double>> ListArray, string Orientation = "row")
        {
            List<int> NumListItems = new List<int>();
            NumListItems.Add(ListArray.Count);
            Nrow = ListArray.Count;
            Ncol = ListArray[0].Count;
            mat = new double[Nrow, Ncol];
            for (int i = 0; i < Nrow; i++)
            {
                NumListItems.Add(ListArray[i].Count);
                if (NumListItems[0] != NumListItems[i]) { throw new InvalidOperationException(); }
                for (int j = 0; j < Ncol; j++)
                {
                    mat[i, j] = ListArray[i][j];
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="scaler"></param>
        public Matrix(double scaler)
        {
            Nrow = 1;
            Ncol = 1;
            mat = new double[1, 1];
            mat[0, 0] = scaler;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        internal Matrix(int row, int col)
        {
            Nrow = row;
            Ncol = col;
            mat = new double[row, col];
        }
        /// <summary>
        /// 
        /// </summary>
        public Matrix()
        {
            Nrow = -1;
            Ncol = -1;
            mat = new double[0, 0];
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="table"></param>
        public Matrix(DataTable table)
        {
            Ncol = table.Columns.Count;
            Nrow = table.Rows.Count;
            mat = new double[Nrow, Ncol];
            for (int i = 0; i < Nrow; i++)
            {
                for (int j = 0; j < Ncol; j++)
                {
                    mat[i, j] = (double)table.Rows[i][j];
                }
            }
        }
        
        #endregion
        #region operators
        /// <summary>
        /// 
        /// </summary>
        /// <param name="array"></param>
        public static implicit operator Matrix(double[,] array)
        {
            return new Matrix(array);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="data"></param>
        public static implicit operator DataTable(Matrix data)
        {
            DataTable table = new DataTable();
            for (int i = 0; i < data.Ncol; i++)
            {
                DataColumn C = new DataColumn();
                C.DataType = typeof(double);
                table.Columns.Add(C);
            }
            for (int i = 0; i < data.Nrow; i++)
            {
                DataRow row = table.NewRow();
                for (int j = 0; j < data.Ncol; j++)
                {
                    row[j] = data[i, j];
                }
                table.Rows.Add(row);
            }
            return table;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="data"></param>
        public static implicit operator Matrix(DataTable data)
        {
            return new Matrix(data);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="MatObj"></param>
        public static implicit operator double[,](Matrix MatObj)
        {
            double[,] arr = new double[MatObj.Nrow, MatObj.Ncol];
            for (int i = 0; i < arr.GetLength(0); i++)
            {
                for (int j = 0; j < arr.GetLength(1); j++)
                {
                    arr[i, j] = MatObj.mat[i, j];
                }
            }
            return arr;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ListArr"></param>
        public static implicit operator Matrix(List<List<double>> ListArr)
        {
            return new Matrix(ListArr);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="matObj"></param>
        public static explicit operator List<List<double>>(Matrix matObj)
        {
            List<List<double>> listarr = new List<List<double>>();
            for (int i = 0; i < matObj.Nrow; i++)
            {
                List<double> innerList = new List<double>();
                listarr.Add(innerList);
                for (int j = 0; j < matObj.Ncol; j++)
                {
                    listarr[i].Add(matObj.mat[i,j]);
                }
            }
            return listarr;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vec"></param>
        public static implicit operator Matrix(Vector vec)
        {
            return new Matrix(vec);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        public static explicit operator Vector (Matrix A)
        {
            if (A.mat.GetLength(1) > 1&& A.mat.GetLength(0)>1 ){ throw new InvalidOperationException(); }
            if (A.mat.GetLength(0) > 1 && A.mat.GetLength(1) == 1|| A.mat.GetLength(0) == 1 && A.mat.GetLength(1) == 1)
            {
                double[] arr = new double[A.mat.GetLength(0)];
                for (int i = 0; i < arr.Length; i++)
                {
                    arr[i] = A.mat[i, 0];
                }
                Vector v = new Vector(arr);
                return v;
            }
            else
            {
                double[] arr = new double[A.mat.GetLength(1)];
                for (int i = 0; i < arr.Length; i++)
                {
                    arr[i] = A.mat[0, i];
                }
                Vector v = new Vector(arr);
                Vector Rv = v.GetRowVector();
                return Rv;
            }
           
        }
      /// <summary>
      /// 
      /// </summary>
      /// <param name="A"></param>
      /// <param name="exponent"></param>
      /// <returns></returns>
        public static Matrix operator ^(Matrix A, double exponent)
        {
            if (exponent == -1) { return Inverse(A); }
            else { return Power(A, exponent); }
        }
        /// <summary>
        /// Matrix Transpose 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix operator ~(Matrix A)
        {
            return Transpose(A);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            return MatrixMult(A,B);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Scaler"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator *(double Scaler, Matrix B)
        {
            return ScalarMatrixMult(Scaler, B);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector operator *(Matrix A ,Vector b)
        {
            return MatrixVectorMult(A, b);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            return MatrixAddition(A, B);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            return MatrixSubtraction (A, B);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix operator -(Matrix A)
        {
            Matrix B = new Matrix();
            foreach (var item in A)
            {
                B.Add(-item);
            }
            return  B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="scaler"></param>
        /// <returns></returns>
        public static Matrix operator /(Matrix A, double scaler)
        {
            return MatrixScalerDivision( A, scaler);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static bool operator==(Matrix A, Matrix B)
        {
            bool isequal = false;
            if (ReferenceEquals(A, B)) { return true; }
            else if(!ReferenceEquals(A, B)) { return false;}
            else
            {
                if ((A.Nrow == B.Nrow) && (A.Ncol == B.Ncol)) { isequal = true; }
                else { return false; }
                for (int i = 0; i < A.Nrow; i++)
                {
                    if (A[i, "col"] == B[i, "col"]) { isequal = true; }
                    else { isequal = false; break; }
                }
            }
            return isequal;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static bool operator !=(Matrix A, Matrix B)
        {
            return A == B;
            //bool isNotequal = true;
            //if (ReferenceEquals(A, B)) { return false; }
            //else if (!ReferenceEquals(A, B)) { return true; }
            //else
            //{
            //    if ((A.Nrow == B.Nrow) && (A.Ncol == B.Ncol)) { isNotequal = false; }
            //    else { return true; }
            //    for (int i = 0; i < A.Nrow; i++)
            //    {
            //        if (A[i, "col"] == B[i, "col"]) { isNotequal = false; }
            //        else { isNotequal = true; break; }
            //    }
            //}
            //return isNotequal;
        }
        #endregion
        #region Indexers
        /// <summary>
        /// 
        /// </summary>
        /// <param name="index"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public Vector this[int index,string type = "col"]
        {
           
            get
            {
                if (type.Equals("col", StringComparison.CurrentCultureIgnoreCase))
                {
                    double[] colArr = new double[this.mat.GetLength(0)];
                    for (int i = 0; i < this.mat.GetLength(0); i++)
                    {
                        colArr[i] = this.mat[i, index];
                    }
                    return new Vector(colArr);
                }
                else
                {
                    double[] rowarr = new double[this.mat.GetLength(1)];
                    for (int i = 0; i < this.mat.GetLength(1); i++)
                    {
                        rowarr[i] = this.mat[index, i];
                    }
                    Vector colv = new Vector(rowarr);
                    return colv.GetRowVector();
                }
            }
            set
            {
                if (type.Equals("col", StringComparison.CurrentCultureIgnoreCase))
                {
                    for (int i = 0; i < this.mat.GetLength(0); i++)
                    {
                        this.mat[i, index] = value[i];
                    }
                  
                }
                else
                {
                    for (int i = 0; i < this.mat.GetLength(1); i++)
                    {
                        this.mat[index, i] = value[i];
                    }
                }
            }
        }


        /// <summary>
        /// Get a column or row vector 
        /// </summary>
        /// <param name="col_Row_Index"></param>
        /// <param name="startP">The starting index of the respective column or row</param>
        /// <param name="endp"> The ending index of the respective column or row</param>
        /// <param name="VectorType"> the orientation; column or row wise</param>
        /// <returns></returns>
        public Vector this[int col_Row_Index, int startP, int endp=0, string VectorType = "col"]
        {
            get
            {
              
                int iter = 0;
                if (VectorType.Equals("col", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (endp == 0) { endp = this.mat.GetLength(0) - 1; }
                    if (endp < startP) { throw new Exception($"starting index {startP} is greater than end index {endp}"); }
                    if (col_Row_Index > this.mat.GetLength(1) || endp > this.mat.GetLength(0)) { throw new Exception(); }
                    double[] arr = new double[endp - startP + 1];
                    for (int i = 0; i < this.mat.GetLength(0); i++)
                    {
                        if (iter == arr.Length) break;
                        iter++;
                        arr[iter - 1] = this.mat[startP + i, col_Row_Index];

                    }
                    return new Vector(arr);
                }

                else
                {
                    if (endp == 0) { endp = this.mat.GetLength(1) - 1; }
                    if (col_Row_Index > this.mat.GetLength(0) || endp > this.mat.GetLength(1)) { throw new Exception(); }
                    double[] arr = new double[endp - startP + 1];
                    for (int i = 0; i < this.mat.GetLength(1); i++)
                    {
                        if (iter == arr.Length) break;
                        iter++;
                        arr[iter - 1] = this.mat[col_Row_Index, startP + i];

                    }
                    return new Vector(arr).GetRowVector();
                }
            }
            set
            {
                int iter = 0;
                if (VectorType.Equals("col", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (endp == 0) { endp = this.mat.GetLength(0) - 1; }
                    if (col_Row_Index > this.mat.GetLength(0) || endp > this.mat.GetLength(1)) { throw new Exception(); }
                    double[] arr = new double[endp - startP + 1];
                    for (int i = 0; i < value.ColVector.GetLength(0); i++)
                    {
                        if (iter == arr.Length) break;
                        iter++;
                        this[ startP+i, col_Row_Index] = value[ i];
                    }
                  
                }

                else
                {
                    if (endp == 0) { endp = this.mat.GetLength(1) - 1; }
                    if (col_Row_Index > this.mat.GetLength(0) || endp > this.mat.GetLength(1)) { throw new Exception(); }
                    double[] arr = new double[endp - startP + 1];
                    for (int i = 0; i < value.RowVector.GetLength(1); i++)
                    {
                        if (iter == arr.Length) break;
                        iter++;
                        this[col_Row_Index,startP+i] =value[ i];

                    }
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="rowIndex"></param>
        /// <param name="colIndex"></param>
        /// <returns></returns>
        public double this[int rowIndex,int colIndex]
        {
            get => this.mat[rowIndex, colIndex];
            set => this.mat[rowIndex, colIndex] = value;
        }
        #endregion
        #region Methods
        #region Basic Matrix Operations
        /// <summary>
        /// 
        /// </summary>
        /// <param name="colOrRowsToModify"></param>
        /// <param name="funs"></param>
        /// <returns></returns>
        public Matrix Modify(int[] colOrRowsToModify,List<Func<Vector>>funs)
        {

            foreach (var i in colOrRowsToModify)
            {
                this[i] = this[i].Get(funs[i]);
            }
            return this;
        }
        /// <summary>
        /// modify the mattrix using the transform function
        /// </summary>
        /// <param name="colOrRowToModify">the column or row to modify</param>
        /// <param name="fun">the modifying function</param>
        /// <returns></returns>
        public Matrix Modify(int colOrRowToModify, Func<Vector> fun)
        { 
            this[colOrRowToModify] = this[colOrRowToModify].Get(fun);
            return this;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="exponent"></param>
        /// <returns></returns>
        public Matrix Power(double exponent)
        {
            if (exponent == 1) { return this.Copy(); }
            else if (exponent == 2) { return this * this; }
            else
            {
                Matrix A = this * this;
                for (int i = 3; i <= exponent; i++)
                {
                  A=  A * this;
                }
                return A;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="exponent"></param>
        /// <returns></returns>
        public static Matrix Power(Matrix A, double exponent)
        {
            if (exponent == 1) { return A.Copy(); }
            else if (exponent == 2) { return A * A; }
            else
            {
                Matrix B = A * A;
                for (int i = 3; i <= exponent; i++)
                {
                    B = B * A;
                }
                return B;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="fun"></param>
        /// <returns></returns>
        public IEnumerable<T> Select<T>(Func<Vector , T> fun)
        {
            Vector v = null!;
            T? res = default(T);
            if (this.orientation.Equals("col", StringComparison.InvariantCultureIgnoreCase))
            {
                for (int i = 0; i < this.Ncol; i++)
                {
                    v = this[i, this.Orientation];
                    res = fun(v);
                    yield return res;
                }
            }
            else
            {
                for (int i = 0; i < this.Nrow; i++)
                {
                    v = this[i, this.Orientation];
                    res = fun(v);
                    yield return res;
                }
            }
            
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public Matrix Apend(Vector v)
        {
            Matrix C = this.Copy();
            Matrix newMatrix = null!;
            if (v.IsColvector)
            {
                newMatrix = new Matrix(C.Nrow, C.Ncol + 1);
                InsectSubMatrix(newMatrix, C, 0, 0);
                newMatrix = InsertVectorColumnwise(newMatrix, v, C.Ncol);
            }
            else
            {
                newMatrix = new Matrix(C.Nrow + 1, C.Ncol);
                InsectSubMatrix(newMatrix, C, 0, 0);
                newMatrix = InsertVectorRowwise(newMatrix, v, C.Nrow);
            }
            return newMatrix;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public Matrix PrePend(Vector v)
        {
            Matrix C = this.Copy();
            Matrix newMatrix = null!;
            if (v.IsColvector)
            {
                newMatrix = new Matrix(C.Nrow, C.Ncol + 1);
                InsectSubMatrix(newMatrix, C, 0, 1);
                newMatrix = InsertVectorColumnwise(newMatrix, v, 0);
            }
            else
            {
                newMatrix = new Matrix(C.Nrow+1,C.Ncol);
                InsectSubMatrix(newMatrix, C, 1, 0);
                newMatrix = InsertVectorRowwise(newMatrix, v, 0);
            }
            return newMatrix;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public Matrix Add(Vector V)
        {
            if (this.Ncol != -1)
            {
                if(V.IsColvector&&(V.Length != this.mat.GetLength(0))) { throw new Exception(); }
                if (!V.IsColvector && (V.Length != this.mat.GetLength(1))) { throw new Exception(); }
            }
            if (V.IsColvector)
            {
                return Addc(V);
            }
            else
            {
                return Addr(V);
            }
        }
        
        private Matrix Addc(Vector V)
        {

            Matrix A = this.Copy();
            double[,] cmat = null!;
            if (A.Ncol == -1 || A.Nrow == -1)
            {
                A.Nrow = 0;A.Ncol = 0;
                
                    cmat = new double[V.Length, A.Ncol + 1];
                    for (int k = 0; k < V.Length; k++)
                    {
                        cmat[k, 0] = V[k];
                    }
               
                
            }
            else
            {
                double[,] M = A.mat;
                    cmat = new double[V.Length, A.Ncol + 1];
                    for (int i = 0; i < cmat.GetLength(1); i++)
                    {
                        for (int k = 0; k < V.Length; k++)
                        {
                            if (i < M.GetLength(1))
                            {
                                cmat[k, i] = M[k, i];
                            }
                            else
                            {
                                cmat[k, i] = V[k];
                            }
                        }
                    }
            }
            A = this;
            A.mat = cmat;
            Nrow = A.mat.GetLength(0);
            Ncol = A.mat.GetLength(1);
            A.orientation = "col";
            return A;
        }
        private Matrix Addr(Vector V)
        {

            Matrix A = this.Copy();
            double[,] cmat = null!;
            if (A.Ncol == -1 || A.Nrow == -1)
            {
                A.Nrow = 0; A.ncol = 0;

                cmat = new double[ A.Ncol + 1, V.Length];
                for (int k = 0; k < V.Length; k++)
                {
                    cmat[0,k] = V[k];
                }
            }
            else
            {
                //for (int i = 0; i < A.Nrow; i++)
                //{
                //    Vector cr = A[i, "row"];
                //    bool allZeros = cr == Vector.Zeros(A.Ncol);
                //    if (allZeros)
                //    {
                //        this[i, "row"] = V;
                //        return this;
                //    }
                //}
                double[,] M = A.mat;
                //M = A.Mat;
                cmat = new double[A.Nrow + 1,V.Length];
                for (int i = 0; i < cmat.GetLength(0); i++)
                {
                    for (int k = 0; k < V.Length; k++)
                    {
                        if (i < M.GetLength(0))
                        {
                            cmat[i, k] = M[i, k];
                        }
                        else
                        {
                            cmat[i, k] = V[k];
                        }
                    }
                }
            }
            A = this;
            A.mat = cmat;
            Nrow = A.mat.GetLength(0);
            Ncol = A.mat.GetLength(1);
            A.Orientation = "row";
            return A;
        }
        /// <summary>
        /// normalizes each column of the matrix
        /// </summary>
        /// <returns></returns>
        public Matrix Normalized()
        {
            for (int i = 0; i < this.Ncol; i++)
            {
                this[i] = this[i].Normalized();
            }
            return this;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="swapfrom"></param>
        /// <param name="swapTo"></param>
        /// <param name="SwapMode"></param>
        public void swap(int swapfrom, int swapTo, SwapState SwapMode= SwapState.ColumnWise)
        {
            string originalorientation = this.orientation;
            if (SwapMode == SwapState.ColumnWise)
            {
                this.orientation = "col";
                Vector c1 = this[swapfrom, orientation];
                Vector c2 = this[swapTo, orientation];
                this[swapfrom, orientation] = c2;
                this[swapTo, orientation] = c1;
                this.orientation = originalorientation;
            }
            else if (SwapMode == SwapState.RowWise)
            {
                this.orientation = "row";
                Vector r1 = this[swapfrom, orientation];
                Vector r2 = this[swapTo, orientation];
                this[swapfrom, orientation] = r2;
                this[swapTo, orientation] = r1;
                this.orientation = originalorientation;
            }
            else
            {
                this.orientation = "row";
                Vector r1 = this[swapfrom, orientation];
                Vector r2 = this[swapTo, orientation];
                this[swapfrom, orientation] = r2;
                this[swapTo, orientation] = r1;
                this.orientation = originalorientation;
                this.orientation = "col";
                Vector c1 = this[swapfrom, orientation];
                Vector c2 = this[swapTo, orientation];
                this[swapfrom, orientation] = c2;
                this[swapTo, orientation] = c1;
                this.orientation = originalorientation;
            }
        }
     /// <summary>
     /// 
     /// </summary>
     /// <param name="colv"></param>
     /// <returns></returns>
        public static int MaxValIndex(Vector colv)
        {
            int index = 0;
            double maxva = colv[0];
            for (int i = 0; i < colv.ColVector.GetLength(0) - 1; i++)
            {
                if (colv[i] > maxva) { maxva = colv[i]; index = i; }
            }
            return index;
        }
        /// <summary>
        /// Get the maximum off diagonal element
        /// </summary>
        /// <returns></returns>
        public double MaxOffDiagonalElement()
        {
            Vector v = new Vector();
            for (int i = 0; i < this.nrow; i++)
            {
                for (int j = 0; j < this.Ncol; j++)
                {
                    if(i!=j)v.Add(this[i,j]);
                }
            }
            //for (int i = 0; i < this.Ncol; i++)
            //{
            //    Vector Ci = this[i];
            //    double maxv = Ci.Where((d => d != Ci[i])).Max<double>();
            //    v.Add(maxv);
            //}
            return v.Max<double>();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix Transpose(Matrix A)
        {
            Matrix Trans = new Matrix(A.Ncol, A.Nrow);
            Trans.Nrow = A.Ncol;
            Trans.Ncol = A.Nrow;
            for (int i = 0; i < Trans.Nrow; i++)
            {
                for (int j = 0; j < Trans.Ncol; j++)
                {
                    Trans.mat [i, j] = A.mat[j, i];
                }
            }
            Trans.Orientation = A.Orientation.Equals("col", StringComparison.InvariantCultureIgnoreCase) ? Trans.Orientation = "row" : Trans.Orientation = "col";
            return Trans;
        }
        /// <summary>
        /// calculate the trace of the matrix
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static double Trace(Matrix A)
        {
            double res = 0.0;
            Matrix d = Matrix.Diagonal(A);
            for (int i = 0; i < d.Nrow; i++)
            {
                Vector rvec = d[i, "row"];
                res +=rvec.Sum();
            }
            return res;
        }
        /// <summary>
        /// calculate determinant
        /// </summary>
        /// <returns></returns>
        public  double Det()
        {
            if (this.Nrow != this.Ncol) throw new InvalidOperationException("matrix must be square");
            if (this.mat.GetLength(0) == 2 && this.mat.GetLength(1) == 2) { return twoByTwoMatrixDet(); }
         var lup=MatrixFactorization.LUP(this);
           Vector rVect = GetDiagonalElements(lup.U);
            double res = 1.0;
            for (int i = 0; i < rVect.ColVector.Length; i++)
            {
                res *= rVect[i];
            }
            int col = this.mat. GetLength(1);
            res = Pow(-1, lup.NumberOfRowExchange) * res;
          
            return res;
        }
        private double twoByTwoMatrixDet()
        {
            double d = (this[0, 0] * this[1, 1]) - (this[0, 1] * this[1, 0]);
            return d;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static double Det(Matrix A)
        {
            if (A.mat.GetLength(0) == 2 && A.mat.GetLength(1) == 2) { return twoByTwoMatrixDet(A); }
            var lup = MatrixFactorization.LUP(A);
            Vector rVect = GetDiagonalElements(lup.U);
            double res = 1.0;
            for (int i = 0; i < rVect.ColVector.Length; i++)
            {
                res *= rVect[i];
            }
            int col = A.mat.GetLength(1);
            res = Pow(-1, lup.NumberOfRowExchange) * res;

            return res;
        }
        private static double twoByTwoMatrixDet(Matrix A)
        {
            double d = (A[0, 0] * A[1, 1]) - (A[0, 1] * A[1, 0]);
            return d;
        }
        private static Matrix twoByTwoMatrixInverse(Matrix A)
        {
            double D = Det(A);
            Matrix LRflip = A.FlipLR();
            Matrix Udflip = LRflip.FlipUD();
            Udflip[0, 1] = -Udflip[0, 1];
            Udflip[1,0] = -Udflip[1, 0];
            Matrix inv = (1.0 / D) * Udflip;
            return inv;
        }
        private static Matrix DiagonalMatrixInverse(Matrix A)
        {
            Matrix Acopy = A.Copy();
            for (int i = 0; i < Acopy.mat.GetLength(0); i++)
            {
                for (int j = 0; j < Acopy.mat.GetLength(1); j++)
                {
                    if (i == j)
                    {
                        Acopy[i, j] = 1.0 / Acopy[i, j];
                    }
                }
            }
            return Acopy;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix Inverse(Matrix A)
        {
            if (A.mat.GetLength(0) == 2 && A.mat.GetLength(1) == 2) { return twoByTwoMatrixInverse(A); }
            if (A.iSLowerTriangular()) { return InverseOfLowerTriangularMatrix(A); }
            if (A.iSUpperTriangular()) { return InverseOfUpperTriangularMatrix(A); }
            if (A.IsDiagonal()) { return DiagonalMatrixInverse(A); }
            var svd = MatrixFactorization.SVD(A);
            Matrix leftVec = svd.LeftSingularVec;
            Matrix rightVec = svd.RightSingularVec;
            Matrix sigma = svd.SingularValues;
            Matrix sigmaTrans = sigma;//^ -1;
            for (int i = 0; i < sigmaTrans.mat.GetLength(0); i++)
            {
                for (int j = 0; j < sigmaTrans.mat.GetLength(1); j++)
                {
                    if (i == j)
                    {
                        if (sigma[i, j] <= 1e-6) { sigmaTrans[i, j] = sigma[i, j]; }
                        else { sigmaTrans[i, j] = Pow(sigma[i, j], -1); }
                    }
                  
                }
            }
            if (A.mat.GetLength(0) > A.mat.GetLength(1))
            {
                Matrix Pinv = rightVec * Transpose(sigmaTrans) * Transpose(leftVec);
                return Pinv;
            }

            Matrix inv = rightVec * sigmaTrans * Transpose(leftVec); 
            return inv;
        }
         static Matrix InverseOfUpperTriangularMatrix(Matrix U)
        {
            Matrix V = new Matrix(new double[U.mat.GetLength(0), U.mat.GetLength(1)]);
            for (int k = U.mat.GetLength(0) - 1; k >= 0; k--)
            {
                V[k, k] = Pow(U[k, k], -1);
                for (int i = k - 1; i >= 0; i--)
                {
                    double sum = 0.0;
                    for (int j = i + 1; j <= k; j++)
                    {
                        sum += U[i, j] * V[j, k];
                    }
                    V[i, k] = -Pow(U[i, i], -1) * sum;
                }
            }
            return V;
        }
         static Matrix InverseOfLowerTriangularMatrix(Matrix L)
        {
            Matrix u = ~L;
            Matrix upp = InverseOfUpperTriangularMatrix(u);
            Matrix lower = ~upp;
            return lower;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
       public bool iSLowerTriangular()
        {
            bool iSlowerT = false;
            double sum = 0.0;
            for (int i = 0; i < this.mat.GetLength(0)-1; i++)
            {
                Vector v= this[i, i + 1, 0, "row"];
              sum += v.Sum();
            }
            if (sum == 0.0) { iSlowerT = true; }
            return iSlowerT;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public bool iSUpperTriangular()
        {
            bool iSuuperT = false;
            double sum = 0.0;
            for (int i = 0; i < this.mat.GetLength(0) - 1; i++)
            {
                sum += this[i, i + 1, 0, "col"].Sum();
            }
            if (sum == 0.0) { iSuuperT = true; }
            return iSuuperT;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Condest()
        {
            try
            {
                double a1 = Matrix.MatrixNorm_1(this);
                int n = this.mat.GetLength(0);
                Matrix eye = IdentityMatrix(n);
                Vector e = eye[0, "col"];
                Vector v = LinearSolvers.SolveQR(this, e);
                double gamma = v.Norm_1();
                Vector epsilon = new Vector(v, (i) => Sign(i));
                Vector x = LinearSolvers.SolveQR(Transpose(this), epsilon);
                if (x.Any((item) => double.IsInfinity(item) || double.IsNaN(item))) { throw new Exception(); }
                bool exit = true;
                int k = 0;
                do
                {
                    int signMaxz = Sign(x.ColVector[0, 0]); double maxz = x[0];
                    for (int i = 0; i < x.ColVector.Length; i++)
                    {
                        if (Abs(x.ColVector[i, 0]) > Abs(maxz))
                        {
                            signMaxz = Sign(x.ColVector[i, 0]);
                            maxz = Abs(x.ColVector[i, 0]);
                        }
                    }
                    int j = x.Get<int>(() => x.IndexOf(signMaxz * x.Max(t => Abs(t))));
                    Vector ej = eye[j, "col"];
                    v = LinearSolvers.SolveQR(this, ej);
                    Vector epsilonj = new Vector(v, (i) => Sign(i));
                    bool isEpsilonEqual = epsilon.Equals(epsilonj);
                    double gamaaNew = gamma;
                    bool CestLtPrevEst = gamma <= gamaaNew;
                    gamma = v.Norm_1();
                    if (isEpsilonEqual || CestLtPrevEst)
                    {
                        x = new Vector(n);
                        for (int i = 0; i < n; i++)
                        {
                            x[i] = Pow(-1, i + 2) * (1.0 + (i / (n - 1)));
                        }
                        x = LinearSolvers.SolveQR(Transpose(this), x);
                        if ((2.0 * x.Norm_1() / (3.0 * n)) > gamma)
                        {
                            v = x.Copy();
                            gamma = 2.0 * x.Norm_1() / (3.0 * n);
                            return gamma * a1;
                        }
                        else
                        {
                            return gamma * a1;
                        }
                    }
                    else
                    {
                        epsilon = new Vector(v, (i) => Sign(i));
                        x = LinearSolvers.SolveQR(Transpose(this), epsilon);
                        k++;
                        if (k >= 2)
                        {
                            //exit = x.Norm_Infinity(100) == x[j] || k > n + 1;
                            exit = x.Norm_1() == Abs(x[j]) || k > n + 1;
                        }

                    }

                } while (exit);


                double res = a1 * gamma;
                return res;
            }
            catch (Exception)
            {

                return double.MaxValue;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double Condest2()
        {
            double a1 = Matrix.MatrixNorm_1(this);
            int n = this.mat.GetLength(0);
            Vector b = new Vector(n);
            double rho = 0.0;
            for (int i = 0; i < b.ColVector.Length; i++)
            {
                b[i] = 1.0 / n;
            }
            do
            {
                Vector x = LinearSolvers.SolveQR(this, b);
                if (x.Norm_1() <= rho)
                {
                    return rho*a1;
                }
                else { rho = x.Norm_1(); }
                Vector y = new Vector(x, (i) => Sign(i));
                Vector z = LinearSolvers.SolveQR(Transpose(this), y);
                int signMaxz = Sign(z.ColVector[0,0]);double maxz = x[0];
                for (int i = 0; i < z.ColVector.Length; i++)
                {
                    if (Abs(z.ColVector[i, 0]) > Abs(maxz))
                    {
                        signMaxz = Sign(z.ColVector[i, 0]);
                        maxz = Abs(z.ColVector[i, 0]);
                    }
                }
                int j = z.Get<int>(() =>z.IndexOf(signMaxz* z.Max(t=>Abs(t))));
                if ( Abs( z[j]) < z.Dot(b))
                {
                    return rho*a1;
                }
                else
                {
                    for (int i = 0; i < b.ColVector.Length; i++)
                    {
                        if (i == j) { b[i] = 1.0; }
                        else { b[i] = 0.0; }
                    }
                }


            } while (true);
          
        }
        /// <summary>
        /// return the condition number of the matrix
        /// </summary>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        public double Cond()
        {
            double Cn = 0.0;
            var svd = MatrixFactorization.SVD(this);
            Vector diag = svd.SingularValues.GetDiagonalElements();
            double maxVal = diag.Max();
            double minVal = diag.Min();
            if (minVal >1e-7)
            {
                Cn = maxVal / minVal;
            }
            else { throw new InvalidOperationException("The Matrix is singular"); }
            return Cn;
        }
        /// <summary>
        /// return the rank of the matrix
        /// </summary>
        /// <returns></returns>
        public int Rank()
        {
            Vector SingVal = MatrixFactorization.SVD(this).SingularValues.GetDiagonalElements();
            int res = SingVal.Where(s => s > 1e-6).Count();
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public Matrix Remove(Vector v)
        {
            Matrix M = new Matrix();
            if (v.IsColvector)
            {
                for (int i = 0; i < this.Ncol; i++)
                {
                    Vector vec = this[i];
                    if (vec == v) { continue; }
                    else { M.Add(vec); }
                }
            }
            else
            {
                for (int i = 0; i < this.Nrow; i++)
                {
                    Vector vec = this[i,"row"];
                    if (vec == v) { continue; }
                    else { M.Add(vec); }
                }
            }
               
            return  M;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ColOrRowIndex"></param>
        /// <param name="vectorType"></param>
        /// <returns></returns>
        public Matrix Remove(int ColOrRowIndex,string vectorType="col")
        {
            Vector vec = this[ColOrRowIndex, vectorType]; 
            return Remove(vec);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public bool isillConditioned()
        {
            bool IsillCond = false;
            double condNumber = this.Cond();
            if (condNumber > 1e12) { IsillCond = true; }
            return IsillCond;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static bool IsSymmetric(Matrix A)
        {
            try
            {
                if (A.mat.GetLength(0) != A.mat.GetLength(1)) { return false; }
                bool sym = true;
                for (int i = 0; i < A.Nrow; i++)
                {
                    if (!sym) break;
                    for (int j = 0; j < A.Ncol; j++)
                    {
                        bool cond1 = Abs(A[i, j] - A[j, i]) < 1e-6;
                        bool cond2 = Sign(A[i, j]) == Sign(A[j, i]);
                        if (cond1 && cond2) { sym = true; }
                        else { sym = false; break; }
                    }

                }
                return sym;
            }
            catch (Exception)
            {

                throw;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public  bool IsDiagonal()
        {
            bool res = false;
            for (int i = 0; i < this.mat.GetLength(0); i++)
            {
                for (int j = 0; j < this.mat.GetLength(0); j++)
                {
                    if (i != j)
                    {
                        if (this[i, j] != 0.0) { res = false;return res; }
                        else { res = true; }
                    }
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public bool IsDiagonallyDominant()
        {
            bool res = false;
            for (int i = 0; i < this.mat.GetLength(0); i++)
            {
                double diff = Abs( this[i, i]) - this[i, "row"].Where((val, j) => j != i).Select(v=>Abs(v)).Sum();
                if (diff >= 0) { res = true; }
                else { res = false;break; }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="Startcol"></param>
        /// <returns></returns>
        public static Matrix SubMatrix(Matrix A, int Startcol)
        {
            int row = A.mat.GetLength(0);
            int col = A.mat.GetLength(1);
            Matrix B = new Matrix(row - Startcol, col - Startcol);
            int k = 0;

            for (int i = 0; i < row; i++)
            {
                if (i >= Startcol)
                {
                    int o = 0;
                    for (int j = 0; j < col; j++)
                    {
                        if (j >= Startcol)
                        {
                            B[k, o] = A[i, j];
                            o++;
                        }
                    }
                    k++;
                }
            }
            return B;
        }
       /// <summary>
       /// 
       /// </summary>
       /// <param name="A"></param>
       /// <param name="startRow"></param>
       /// <param name="Startcol"></param>
       /// <returns></returns>
        public static Matrix SubMatrix(Matrix A, int startRow, int Startcol)
        {
            int row = A.mat.GetLength(0);
            int col = A.mat.GetLength(1);
            Matrix B = new Matrix(row - startRow, col - Startcol);
            int k = 0;

            for (int i = 0; i < row; i++)
            {
                if (i >= startRow)
                {
                    int o = 0;
                    for (int j = 0; j < col; j++)
                    {
                        if (j >= Startcol)
                        {
                            B[k, o] = A[i, j];
                            o++;
                        }
                    }
                    k++;
                }
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="startRow"></param>
        /// <param name="Startcol"></param>
        /// <param name="endrow"></param>
        /// <param name="endcol"></param>
        /// <returns></returns>
        public static Matrix SubMatrix(Matrix A, int startRow, int Startcol,int endrow,int endcol)
        {
            int row = A.mat.GetLength(0);
            int col = A.mat.GetLength(1);
            Matrix B = new Matrix(endrow - startRow+1, endcol - Startcol+1);
            int k = 0;

            for (int i = 0; i <= endrow; i++)
            {
                if (i >= startRow)
                {
                    int o = 0;
                    for (int j = Startcol; j <= endcol; j++)
                    {
                            B[k, o] = A[i, j];
                        o++;
                    }
                    k++;
                    
                }
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="V"></param>
        /// <param name="colindex"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix InsertVectorColumnwise(Matrix A, Vector V, int colindex)
        {
            if(colindex > A.mat.GetLength(1) - 1) { throw new Exception($"{colindex} is greater than  {A.mat.GetLength(1) - 1}, the upper column bound of the matrix"); }
            if (V.ColVector.GetLength(0) > A.mat.GetLength(0)) { throw new Exception("Inconsistent dimensions"); }

            for (int i = 0; i < A.mat.GetLength(0); i++)
            {
                A[i, colindex] = V[i];
            }
            return A;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="V"></param>
        /// <param name="rowindec"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix InsertVectorRowwise(Matrix A,Vector V, int rowindec)
        {
            if (rowindec > A.mat.GetLength(0) - 1) { throw new Exception($"{rowindec} is greater than the {A.mat.GetLength(0) - 1}, the upper row bound of the matrix"); }
            if (V.RowVector.GetLength(1) > A.mat.GetLength(1)) { throw new Exception("Inconsistent dimensions"); }
            for (int i = 0; i < A.mat.GetLength(1); i++)
            {
                A[rowindec, i] = V[i];
            }
            return A;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="V"></param>
        /// <param name="colindex"></param>
        /// <param name="startIndex"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix InsertVectorColumnwise(Matrix A, Vector V, int colindex,int startIndex)
        {
            if (A.mat.GetLength(0) - startIndex + 1 > V.ColVector.GetLength(0)) { throw new Exception(); }
            for (int i = 0; i <V.ColVector.GetLength(0); i++)
            {
                A[i+startIndex, colindex] = V[i];
            }
            return A;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="V"></param>
        /// <param name="rowindec"></param>
        /// <param name="startIndex"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix InsertVectorRowwise(Matrix A, Vector V, int rowindec,int startIndex)
        {
            if (A.mat.GetLength(1) - startIndex + 1 > V.RowVector.GetLength(1)) { throw new Exception(); }
            for (int i = 0; i < V.RowVector.GetLength(1) ; i++)
            {
                A[rowindec, i+startIndex] = V[i];
            }
            return A;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="colIndex"></param>
        /// <returns></returns>
        public static Vector ExtractColumnVector(Matrix A, int colIndex)
        {
            Vector B = new Vector();
            for (int i = 0; i < A.mat.GetLength(0); i++)
            {

                for (int j = 0; j < A.mat.GetLength(1); j++)
                {
                    if (j == colIndex)
                    {
                        B[i] = A[i, j];
                        break;
                    }
                }
            }
            return B;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="exponent"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        public static Matrix TwPow(Matrix A, double exponent)
        {
            Matrix B = A.Copy();
            for (int i = 0; i < B.mat.GetLength(0); i++)
            {
                for (int j = 0; j < B.mat.GetLength(1); j++)
                {
                    if (Sign(exponent) == -1)
                    {
                        if (B[i, j] != 0.0) { B[i, j] = Math.Pow(B[i, j], exponent); }
                        else if(B[i,j]<0.0){ throw new InvalidOperationException(); }
                    }
                    else
                    {
                        B[i, j] = Math.Pow(B[i, j], exponent);
                    }

                }
            }
            return B;
        }
      /// <summary>
      /// 
      /// </summary>
      /// <returns></returns>
        public Matrix Copy()
        {
            return (Matrix)Clone();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static   Vector MatrixVectorMult(Matrix A, Vector b)
        {
            return (Vector)MatrixMult(A, b);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public  Vector MatrixVectorMult( Vector b)
        {
            return (Vector)MatrixMult(this, b);
        }
        /// <summary>
        /// Calculates the outer product of two vectors
        /// </summary>
        /// <param name="V1"> A column Vector</param>
        /// <param name="V2"> A row vector</param>
        /// <returns></returns>
        public static Matrix OuterProduct(Vector V1,Vector V2)
        {

            Matrix res = new Matrix(V1.Length, V2.Length);
            for (int i = 0; i < V1.Length; i++)
            {
                    for (int k = 0; k < V2.Length; k++)
                    {
                        res[i, k] = V1.ColVector[i,0] * V2.RowVector[0,k];
                    }

            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Matrix FlipLR()
        {
            Matrix Flip = this.Copy();
            int C = this.mat.GetLength(1);
            if (C % 2 == 0)
            {
                for (int i = 0; i < Flip.mat.GetLength(1); i++)
                {
                    Vector v = this[i, "col"];
                    Vector vv = this[C - 1 - i, "col"];
                    Flip[C - 1 - i, "col"] = v;
                    Flip[i, "col"] = vv;
                }
            }
            else
            {

                for (int i = 0; i < Flip.mat.GetLength(1); i++)
                {
                    if(i!=(C/2) )
                    {
                        Vector v = this[i, "col"];
                        Vector vv = this[C - 1 - i, "col"];
                        Flip[C - 1 - i, "col"] = v;
                        Flip[ i, "col"] = vv;
                    }
                   
                }
            }
            return Flip;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Matrix FlipUD()
        {
            Matrix Flip = this.Copy();
            int C = this.mat.GetLength(0);
            if (C % 2 == 0)
            {
                for (int i = 0; i < Flip.mat.GetLength(0); i++)
                {
                    Vector v = this[i, "row"];
                    Vector vv = this[C - 1 - i, "row"];
                    Flip[C - 1 - i, "row"] = v;
                    Flip[i, "row"] = vv;
                }
            }
            else
            {

                for (int i = 0; i < Flip.mat.GetLength(1); i++)
                {
                    if (i != (C / 2))
                    {
                        Vector v = this[i, "row"];
                        Vector vv = this[C - 1 - i, "row"];
                        Flip[C - 1 - i, "row"] = v;
                        Flip[i, "row"] = vv;
                    }

                }
            }
            return Flip;
        }
        /// <summary>
        /// return the kronecker matrix
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Matrix KroneckerProduct(Matrix A, Matrix B)
        {

            int DCcol = A.Ncol* B.Ncol;
            int Drow = A.Nrow * B.Nrow;
            Matrix D = new Matrix(Drow, DCcol);
            int o = 0;
           
            for (int i = 0; i < A.mat.GetLength((0)); i++)
            {
                double scaler = 0.0;
                int p = 0;
                for (int j = 0; j < A.mat.GetLength(1); j++)
                {
                    scaler = A[i, j];
                    Matrix sub = scaler * B;
                    InsectSubMatrix(D,sub,i+o,j+p);
                    p += B.Ncol-1;
                }
                o+= B.Nrow-1 ;
            }
            return D;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="S"></param>
        /// <param name="StartRowIndex"></param>
        /// <param name="startColIndex"></param>
        public static void InsectSubMatrix(Matrix A, Matrix S,int StartRowIndex,int startColIndex)
        {
            Matrix Acopy = A.Copy();
            for (int i = 0; i < S.mat.GetLength(0); i++)
            {
                for (int j = 0; j < S.mat.GetLength(1); j++)
                {
                    A[StartRowIndex + i, startColIndex + j] = S[i, j];
                }
            }
        }
        /// <summary>
        /// termwise or Schur product of two matrices
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix HadamardProduct(Matrix A, Matrix B)
        {
            return TwMult(A, B);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix MatrixMult(Matrix A, Matrix B)
        {

            bool isConformable = A.Ncol == B.Nrow;
            if (!isConformable) throw new Exception("Matrices are not conformable");
            double[,] C = new double[A.Nrow, B.Ncol];
            for (int i = 0; i < A.Nrow; i++)
            {
                for (int j = 0; j < B.Ncol; j++)
                {
                    C[i, j] = Vector.Dot(A[i, "row"], B[j, "col"]);
                }
            }
            return C;
        }
        
         static Matrix TwMult(Matrix A, Matrix B)
        {

            bool isConformable = A.Ncol*A.nrow == B.Nrow*B.ncol;
            if (!isConformable) throw new Exception("Matrices are not conformable");
            double[,] C = new double[A.Nrow, B.Ncol];
            for (int i = 0; i < A.Nrow; i++)
            {
                for (int j = 0; j < B.Ncol; j++)
                {
                    C[i, j] = A[i,j]*B[i,j];
                }
            }
            return C;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix TwDiv(Matrix A, Matrix B)
        {

            bool isConformable = A.Ncol * A.nrow == B.Nrow * B.ncol;
            if (!isConformable) throw new Exception("Matrices are not conformable");
            double[,] C = new double[A.Nrow, B.Ncol];
            for (int i = 0; i < A.Nrow; i++)
            {
                for (int j = 0; j < B.Ncol; j++)
                {
                    C[i, j] = A[i, j] / B[i, j];
                }
            }
            return C;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix MatrixMult(params Matrix[] A)
        {
            Matrix res = MatrixMult(A[0], A[1]);
            for (int i = 2; i < A.Length; i++)
            {
                res = MatrixMult(res, A[i]);
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
       public bool IsPositiveDefinite()
        {
            bool res = false;
            if (!IsSymmetric(this))
            {
                return res;
            }
            else
            {
                try
                {
                  var Chol=  MatrixFactorization.Cholesky(this);
                    res = true;
                }
                catch (Exception)
                {

                    return false;
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="scalar"></param>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix ScalarMatrixMult(double scalar, Matrix A)
        {
            Matrix B = A.Copy();
            for (int i = 0; i < A.mat.GetLength(0); i++)
            {
                for (int j = 0; j < A.mat.GetLength(1); j++)
                {
                    B[i, j] = scalar * A[i, j];
                }
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="scalar"></param>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix ScalarMatrixDiv(double scalar, Matrix A)
        {
            Matrix B = A.Copy();
            for (int i = 0; i < A.mat. GetLength(0); i++)
            {
                for (int j = 0; j < A.mat. GetLength(1); j++)
                {
                    B[i, j] = A[i, j] / scalar;
                }
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        /// <exception cref="InvalidDataException"></exception>
        public static Matrix MatrixAddition(Matrix A, Matrix B)
        {
            if (A.Nrow * A.Ncol != B.Nrow * B.Ncol) { throw new InvalidDataException("Matrices must be the same dimension"); }
            Matrix C = new Matrix(A.mat.GetLength(0), A.mat.GetLength(1));
           
                    for (int j = 0; j < A.mat.GetLength(0); j++)
                    {
                        for (int k = 0; k < A.mat.GetLength(1); k++)
                        {
                            C[j, k] = A[j, k] + B[j, k];
                        }
                    }
            return C;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix MatrixAddition(params Matrix[] A)
        {
            Matrix B = MatrixAddition(A[0], A[1]);
            for (int i = 2; i < A.Length; i++)
            {
                B = MatrixAddition(B, A[i]);
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        /// <exception cref="InvalidDataException"></exception>
        public static Matrix MatrixSubtraction(Matrix A, Matrix B)
        {
            if (A.Nrow * A.Ncol != B.Nrow * B.Ncol) { throw new InvalidDataException("Matrices must be the same dimension"); }
            Matrix C = new Matrix(A.mat.GetLength(0), A.mat.GetLength(1));

            for (int j = 0; j < A.mat.GetLength(0); j++)
            {
                for (int k = 0; k < A.mat.GetLength(1); k++)
                {
                    C[j, k] = A[j, k] - B[j, k];
                }
            }
            return C;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix MatrixSubtraction(params Matrix[] A)
        {
            Matrix B = MatrixSubtraction(A[0], A[1]);

            for (int i = 2; i < A.Length; i++)
            {
                B = MatrixSubtraction(B, A[i]);
            }
            return B;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="scaler"></param>
        /// <returns></returns>
        public static Matrix MatrixScalerDivision(Matrix A, double scaler)
        {

            Matrix C = new Matrix(A.mat.GetLength(0), A.mat.GetLength(1));

            for (int j = 0; j < A.mat.GetLength(0); j++)
            {
                for (int k = 0; k < A.mat.GetLength(1); k++)
                {
                    C[j, k] = A[j, k] / scaler;
                }
            }
            return C;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="PivotRowIndex"></param>
        /// <param name="maxValIndex"></param>
        public void SwapRowVectors(int PivotRowIndex, int maxValIndex)
        {
            //Vector r1 = this[PivotRowIndex,PivotRowIndex,0, "row"];
            //Vector r2 = this[maxValIndex,PivotRowIndex,0, "row"];
            Vector r1 = this[PivotRowIndex, 0, 0, "row"];
            Vector r2 = this[maxValIndex, 0, 0, "row"];
            this[PivotRowIndex,0,0, "row"] = r2;
            this[maxValIndex,0,0, "row"] = r1;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double MaxElement()
        {
            double maxV = Abs( this.mat[0, 0]);
            for (int i = 0; i < this.mat.GetLength(0); i++)
            {
                double M = this[i, "row"].Max(v => Abs(v));
                if (M > maxV) { maxV = M; }
            }
            return maxV;
        }
       
        #endregion
        #region Orthogonization ,Reflectors and Rotators
        internal static Matrix GramSchmidtOthogonalization(Matrix A)
        {
            Matrix B = A.Copy();
            //Matrix res = new Matrix(A.mat.GetLength(0), A.mat.GetLength(1));
            //List<Vector> q = new List<Vector>();
            Matrix q = new Matrix();
            for (int i = 0; i < B.mat.GetLength(1); i++)
            {
                //Vector C = new Vector();
                Vector ColVec = B[i];
                ////for (int j = 0; j < B.mat.GetLength(0); j++)
                ////{
                ////    C.Add(B[j, i]);
                ////}
                q.Add(ColVec);
                for (int j = 0; j <= i - 1; j++)
                {
                    double qc = Vector.Dot(q[j].GetRowVector(), ColVec); 
                    double qq = Vector.Dot(q[j].GetRowVector(), q[j]); ;
                    double fourierConstant = qc / qq;
                   Vector projectedVector = fourierConstant * q[j];
                    Vector orthogonalComplement = ColVec - projectedVector;
                    q.Add(orthogonalComplement);
                    //for (int k = 0; k < B.mat.GetLength(0); k++)
                    //{
                    //    q[i][k] = q[i][k] - sss[k];
                    //}
                }

            }

            //for (int i = 0; i < q.Ncol; i++)
            //{
            //    //q[i] = q[i].Select(item => item / q[i].Norm()).ToList();
            //    q[i] = q[i].Normalized();
            //    //for (int j = 0; j < q[i].Count(); j++)
            //    //{
            //    //    res[j, i] = q[i][j];
            //    //}
            //}

            return q.Normalized();
        }
        internal static Matrix ModifiedGramSchmidtOthogonalization(Matrix A)
        {
            Matrix B = A.Copy();
            //Matrix res = new Matrix(A.mat.GetLength(0), A.mat.GetLength(1));
            //List<Vector> q = new List<Vector>();
            Matrix q = new Matrix();
            for (int i = 0; i < B.mat.GetLength(1); i++)
            {
                //extract the column vector
                Vector ColumnVector = B[i, 0, 0, "col"];
              
                for (int j = 0; j <= i - 1; j++)
                {
                    //calculate the projection along the q[k] direction
                    double dotP = Vector.Dot(q[j].GetRowVector(), ColumnVector);
                    Vector projection = dotP * q[j];
                    ColumnVector = ColumnVector - projection;
                }
                //ColumnVector = ColumnVector / ColumnVector.Norm();
                q.Add(ColumnVector.Normalized());
            }

            //for (int i = 0; i < q.Count; i++)
            //{
            //    for (int j = 0; j < q[i].Count(); j++)
            //    {
            //        res[j, i] = q[i][j];
            //    }
            //}

            return q;
        }
        private static Matrix HouseHolderReflector(Matrix A, Matrix identitiMat, int colIndex)
        {
            Matrix seye = SubMatrix(identitiMat, colIndex);
            Vector e = seye[0,"col"];
            Vector z = A[0,"col"];
            if (z.Norm() == 0.0)
            {
               z= Vector.Ones(z.Length);
                //for (int i = 0; i < z.ColVector. GetLength(0); i++)
                //{
                //    z.ColVector[i, 0] = 1.0;
                //}
            }
            
            Vector u = HouseHolderVector(z, out double gam);
            Matrix uuT = Matrix.OuterProduct(u, u.GetRowVector());
            uuT = gam * uuT;// ScalarMatrixMult(gam, uuT);
            Matrix H = seye - uuT;
            return H;
        }
        private static Vector HouseHolderVector(Vector z)
        {
            if (z[0] == 0.0) { z[0] = z[1]; }
            double znorm = z.Norm();
            double a1 = -Sign(z[0]) * znorm;
           Vector u = new Vector();
            u.Add( Sqrt((a1 - z[0]) / (2 * a1)));
            for (int n = 1; n < z.ColVector.GetLongLength(0); n++)
            {
                double term = z[n] / (-2 * a1 * u[0]);
                u.Add(term);
            }
            return u;
        }
        private static Vector HouseHolderVector(Vector z, out double gamma)
        {
            if (z[0] == 0.0) { z[0] = 1e-6; }
            Vector u = new Vector(z.ColVector.GetLength(0));
            Vector x =z[0,"col"];
            double beta = x.Max(m => Abs(m));
            double tau = 0.0;
            if (beta == 0) { gamma = 0.0; }
            else
            {
                x = x.Select(item => item / beta).ToArray();
                tau = Sqrt(x.Select(s => Pow(s, 2)).Sum());
                if (Sign(x[0]) == -1) { tau = -tau; }
                x[0] = tau + x[0];
                gamma = x[0] / tau;
                for (int n = 1; n < x.Length; n++)
                {
                    x[n] = x[n] / x[0];
                    u[n] = x[n];
                }
                x[0] = 1.0;
                u[0] = x[0];
                tau = beta * tau;
            }
       
            return u;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static (Matrix UpperTriangle, List<Matrix> ListofHouseHolderMatrix) HouseHolder(Matrix A)
        {
            Matrix Hmat = A.Copy();
            Matrix Smat = Hmat.Copy();
            Matrix eye = IdentityMatrix(Hmat.mat.GetLength(0));
            int row = Hmat.mat.GetLength(0);
            int col = Hmat.mat.GetLength(1);
            if (row > col) { col++; }
            List<Matrix> Hlist = new List<Matrix>();
            for (int k = 0; k < col - 1; k++)//revisit the less than or equal sign. this is very important
            {
                double[,] OrigSmart = Smat.Copy();
                Smat = SubMatrix(OrigSmart, k);
                Matrix H = HouseHolderReflector(Smat, eye, k);
                H = Hmodified(H, k);
                Hlist.Add(H);
                Smat = MatrixMult(H, OrigSmart);
            }

            Matrix Hmodified(Matrix Hm, int index)
            {
                Matrix B = IdentityMatrix(Hm.mat.GetLength(0) + index);
                int k = -index;
                for (int i = 0; i < B.mat.GetLength(0); i++)
                {
                    int o = 0;
                    for (int j = 0; j < B.mat.GetLength(1); j++)
                    {
                        if (i >= index && j >= index)
                        {
                            B[i, j] = Hm[k, o];
                            o++;
                        }
                    }
                    k++;
                }
                return B;
            }
            return (Smat, Hlist);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Gmat"></param>
        /// <param name="PivotIndex"></param>
        /// <param name="EliminationIndex"></param>
        /// <returns></returns>
        public static Matrix GivensMatrix(Matrix Gmat, int PivotIndex, int EliminationIndex)
        {
            Matrix G = IdentityMatrix(2);
            double Cc = Gmat.mat[PivotIndex, PivotIndex];
            double Ss = Gmat.mat[EliminationIndex, PivotIndex];
            if (Cc == 0.0 && Ss == 0.0)
            {
                return G;
            }
            double C = Cc / Sqrt((Pow(Cc, 2)) + Pow(Ss, 2));
            double S = Ss / Sqrt((Pow(Cc, 2)) + Pow(Ss, 2));

            G[0, 0] = C;
            G[0, 1] = S;
            G[1, 0] = -S;
            G[1, 1] = C;
            return G;
        }
        #endregion
        #region matrix norms
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static double MatrixNorm_1(Matrix A)
        {
            List<double> colSum = new List<double>();
            for (int i = 0; i < A.mat.GetLength(1); i++)
            {
                Vector v = A[i, 0, 0, "col"];
                colSum.Add(v.Sum(d => Abs(d)));
            }
            return colSum.Max();
          
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static double MatrixNorm_Infinity(Matrix A)
        {
            List<double> rowSum = new List<double>();
            for (int i = 0; i < A.mat.GetLength(0); i++)
            {
                Vector v = A[i, 0, 0, "row"];
                rowSum.Add(v.Sum(d => Abs(d)));
            }
            return rowSum.Max();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static double FrobeniusNorm(Matrix A)
        {
            double Sum = 0.0;
            for (int i = 0; i < A.mat.GetLength(1); i++)
            {
                Sum += A[i].Pow(2).Sum();
            }
            Sum = Sqrt(Sum);
            return Sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Vector GetDiagonalElements(Matrix A)
        {
            Vector res = new Vector();
            for (int i = 0; i < A.mat.GetLength(0); i++)
            {
                for (int j = 0; j < A.mat.GetLength(1); j++)
                {
                    if (i == j) res.Add(A[i, j]);
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public  Vector GetDiagonalElements()
        {
            Vector res = new Vector();
            for (int i = 0; i < this.mat.GetLength(0); i++)
            {
                for (int j = 0; j < this.mat.GetLength(1); j++)
                {
                    if (i == j) res.Add(this[i, j]);
                }
            }
            return res;
        }

        //public static double MatrixNorm_2(double[,] A)
        //{
        //    double maxEigenValue = PowerMethod(MatrixMult(Transpose(A), A)).lamda;
        //    return Sqrt(maxEigenValue);
        //}
        #endregion
        #endregion
        #region Special Matrices
        /// <summary>
        /// 
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns></returns>
        public static Matrix IdentityMatrix(int row, int col)
        {
            double[,] eyeMat = new double[row, col];
            for (int i = 0; i < eyeMat.GetLength(0); i++)
            {
                for (int j = 0; j < eyeMat.GetLength(1); j++)
                {
                    if (i == j) { eyeMat[i, j] = 1.0; }
                }
            }
            Matrix eye = new Matrix(eyeMat);
            eye.Nrow = row;
            eye.Ncol = col;
            return eye;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim"></param>
        /// <returns></returns>
        public static Matrix IdentityMatrix(int dim)
        {
            double[,] eyeMat = new double[dim, dim];
            for (int i = 0; i < eyeMat.GetLength(0); i++)
            {
                for (int j = 0; j < eyeMat.GetLength(1); j++)
                {
                    if (i == j) { eyeMat[i, j] = 1.0; }
                }
            }
            Matrix eye = new Matrix(eyeMat);
            eye.Nrow = dim;
            eye.Ncol = dim;
            return eye;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns></returns>
        public static Matrix Zeros(int row, int col)
        {
            return new Matrix(row, col);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        public static Matrix Diagonal(int dim,double val)
        {
            Matrix D = new Matrix(dim, dim);
            for (int i = 0; i < D.Nrow; i++)
            {
                for (int j = 0; j < D.Ncol; j++)
                {
                    if (i == j)
                    {
                        D[i, j] = val;
                    }
                }
            }
            return D;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix Diagonal(Matrix A)
        {
            Matrix D = new Matrix(A.Nrow, A.Ncol);
            for (int i = 0; i < D.Nrow; i++)
            {
                for (int j = 0; j < D.Ncol; j++)
                {
                    if (i == j)
                    {
                        D[i, j] = A[i, j];
                    }
                }
            }
            return D;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Matrix Diagonal(Vector b)
        {
            if (b.IsColvector)
            {
                Matrix D = new Matrix(b.ColVector.Length, b.ColVector.Length);
                for (int i = 0; i < D.Nrow; i++)
                {
                    for (int j = 0; j < D.Ncol; j++)
                    {
                        if (i == j)
                        {
                            D[i, j] = b.ColVector[i, 0];
                        }
                    }
                }
                return D;
            }
            else
            {
                Matrix D = new Matrix(b.RowVector.Length, b.RowVector.Length);
                for (int i = 0; i < D.Nrow; i++)
                {
                    for (int j = 0; j < D.Ncol; j++)
                    {
                        if (i == j)
                        {
                            D[i, j] = b.RowVector[0, i];
                        }
                    }
                }
                return D;
            }
            
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="rowOrcol"></param>
        /// <returns></returns>
        public static Vector vec(Matrix A,string rowOrcol="col")
        {
            Vector V = new Vector();
            if (rowOrcol.Equals("col", StringComparison.CurrentCultureIgnoreCase))
            {
                for (int i = 0; i < A.mat.GetLength(1); i++)
                {
                    Vector vi = A[i, "col"];
                    for (int j = 0; j < vi.ColVector.GetLength(0); j++)
                    {
                        V.Add(vi[j]);
                    }

                }
            }
            else
            {
                for (int i = 0; i < A.mat.GetLength(0); i++)
                {
                    Vector vi = A[i, "row"];
                    for (int j = 0; j < vi.RowVector.GetLength(1); j++)
                    {
                        V.Add(vi[j]);
                    }

                }
            }
               
            return V;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Xval"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        public static Matrix VandermondeMatrix(Vector Xval, int degree)
        {
            Matrix Vander = new Matrix(Xval.Length, degree+1);
            for (int i = 0; i < Vander.mat.GetLength(0); i++)
            {
                for (int j = 0; j < Vander.mat.GetLength(1); j++)
                {
                        Vander.mat[i, j] = Pow(Xval[i], j);
                }
            }
            return Vander;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="XList"></param>
        /// <param name="includeProductterm"></param>
        /// <returns></returns>
        public static Matrix VandermondeMatrix(List<Vector> XList, bool includeProductterm = false)
        {
            Matrix Vander = new Matrix();
            Vector C1 = Vector.Ones(XList[0].Length);
            Vander.Add(C1);
            Vector Prod = C1;
            for (int i = 0; i < XList.Count; i++)
            {
                if (includeProductterm) { Prod *= XList[i]; }
                Vander.Add(XList[i]);
            }
            if (includeProductterm) { Vander.Add(Prod); }
            return Vander;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="size"></param>
        /// <returns></returns>
        public static Matrix Hilb(double size)
        {
            Matrix res = new Matrix((int)size, (int)size);
            for (int i = 0; i < res.mat.GetLength(0); i++)
            {

                for (int j = 0; j < res.mat.GetLength(1); j++)
                {
                    res[i, j] = size / (size * (j + i + 1));
                }
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public static Matrix Hessenberg(Matrix A)
        {
            if (A.IsDiagonal()) { return A; }
            Matrix? HessMat = A.Clone() as Matrix;
            for (int k = 0; k < A.mat.GetLength(1) - 2; k++)
            {
                Matrix SubMat = SubMatrix(HessMat!, k + 1, k);
                Vector x = SubMat[0, 0, 0, "col"];
                var HouseParams = HouseVectors(x);
                Vector u = HouseParams.Hvec;
                double beta = HouseParams.Beta;
                Matrix outer = beta * OuterProduct(u, u.GetRowVector());
                Matrix eye = IdentityMatrix(outer.mat.GetLength(0));
                outer = eye - outer;
                Matrix HH = outer * SubMat;
                HH = Hmodified(outer, k + 1);
                Matrix Anew = HH * HessMat!;
                HessMat = Anew * Transpose(HH);
            }

            (Vector Hvec, double Beta) HouseVectors(Vector X)
            {
                Vector u = new Vector(X.ColVector.GetLength(0));
                double maxX = X.Max(item => Abs(item));
                Vector XNorm = Vector.TwDivision(X, maxX);
                Matrix id = IdentityMatrix(X.Length);
                Vector idFirstColume = id[0, 0, 0, "col"];
                bool signOfFirstElement = Sign(XNorm[0]) == 1 || XNorm[0] == 0.0 ? true : false;
                if (signOfFirstElement)
                {
                    u = XNorm + (XNorm.Norm() * idFirstColume);
                }
                else
                {
                    u = XNorm - (XNorm.Norm() * idFirstColume);
                }
                double beta = 2.0 / u.Dot(u);
                return (u, beta);
            }

            Matrix Hmodified(Matrix Hm, int index)
            {
                Matrix B = IdentityMatrix(Hm.mat.GetLength(0) + index);
                int k = -index;
                for (int i = 0; i < B.mat.GetLength(0); i++)
                {
                    int o = 0;
                    for (int j = 0; j < B.mat.GetLength(1); j++)
                    {
                        if (i >= index && j >= index)
                        {
                            B[i, j] = Hm[k, o];
                            o++;
                        }
                    }
                    k++;
                }
                return B;
            }
            return HessMat!;
        }
        #endregion
        #region Fields and properties

        private double[,] mat;
        private int nrow;
        private int ncol;
        private string orientation = "col";
        //public double[,]? Mat { get => mat; set => mat = value!; }
        /// <summary>
        /// 
        /// </summary>
        public int Nrow { get => nrow; set => nrow = value; }
        /// <summary>
        /// 
        /// </summary>
        public int Ncol { get => ncol; set => ncol = value; }
        /// <summary>
        /// 
        /// </summary>
        public string Orientation { get => orientation; set => orientation = value; }
        #endregion
        #region Implemented Interfaces
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public IEnumerator<Vector> GetEnumerator()
        {
            if (orientation.Equals("col", StringComparison.CurrentCultureIgnoreCase))
            {
                for (int i = 0; i < this.Ncol; i++)
                {
                        yield return this[i, orientation];
                }
            }
            else
            {
                for (int i = 0; i < this.Nrow; i++)
                {
                    yield return this[i, orientation];
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        IEnumerator IEnumerable.GetEnumerator()
        {
            for (int i = 0; i < this.Nrow; i++)
            {
                if (orientation.Equals("col", StringComparison.CurrentCultureIgnoreCase))
                {
                    yield return this[i, orientation];
                }
                else
                {
                    yield return this[i, orientation];
                }
            }
        }
        /// <summary>
        /// perform a shllow copy of the matrix
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
         
            using (MemoryStream str = new MemoryStream())
            {
                //BinaryFormatter f = new BinaryFormatter();
                //f.Serialize(str, this);
                //str.Position = 0;
                //Matrix Cmatrix = (Matrix)f.Deserialize(str);
                Matrix Cmatrix = null!;
                if (this.Nrow == -1 || this.Ncol == -1) { Cmatrix = new Matrix(); }
                else { Cmatrix = new Matrix(this.Nrow, this.Ncol); }
                 
                for (int i = 0; i < this.nrow; i++)
                {
                    for (int j = 0; j < ncol; j++)
                    {
                        Cmatrix[i,j]=this[i,j];
                    }
                   
                }
                Cmatrix.Nrow = this.Nrow;
                Cmatrix.ncol = this.Ncol;
                Cmatrix.Orientation = this.Orientation;
            
                    return Cmatrix;
            }
         
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Matrix? other)
        {
            bool isequal = false;
            if (ReferenceEquals(this, other)) { return true; }
            else if (!ReferenceEquals(this, other)) { return false; }
            else
            {
                if ((this.Nrow == other.Nrow) && (this.Ncol == other.Ncol)) { isequal = true; }
                else { return false; }
                for (int i = 0; i < this.Nrow; i++)
                {
                    if (this[i, "col"] == other[i, "col"]) { isequal = true; }
                    else { isequal = false; break; }
                }
            }
            return isequal;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object? obj)
        {
            Matrix A = (Matrix)obj!;
            return this.Equals(A);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            List<int> intList = this.GetDiagonalElements().Select(s => (int)Ceiling( s)).ToList();
            int sum = intList.Sum();
            if (sum > 1000) { sum = 999; }
            int bas = this.GetHashCode() * 17*sum;
            int bas2 = (Ncol.GetHashCode() * Nrow.GetHashCode()) * Ncol;
            return bas + bas2; ;
        }
        #endregion
    }
}
