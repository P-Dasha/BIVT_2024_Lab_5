using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data.Common;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static Program;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }

    public int Factorial(int n)
    {
        int answer = 1;
        for (int i = 2; i <= n; i++)
        {
            answer *= i;
        }
        return answer;
    }

    public long Combinations(int n, int k)
    {
        return Factorial(n) / (Factorial(k) * Factorial(n - k));
    }

    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
        if (n < k || k < 0 || n < 0)
            return 0;

        answer = Combinations(n, k);
        // end

        return answer;
    }

    public bool IsTriangle(double a, double b, double c)
    {
        if (a + b > c && a + c > b && b + c > a) 
            return true;

        return false;
    }

    public double GeronArea(double a, double b, double c)
    {
        double p = (a + b + c) / 2;

        return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
    }

    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];

        if (!IsTriangle(a1, b1, c1) || !IsTriangle(a2, b2, c2)) return -1;

        double area1 = GeronArea(a1, b1, c1), area2 = GeronArea(a2, b2, c2);

        if (area1 > area2) answer = 1;
        else if (area1 < area2) answer = 2;
        else answer = 0;

        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }

    public double GetDistance(double v, double a, int t)
    {
        return v * t + a * t * t / 2;
    }

    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here
        double s1 = GetDistance(v1, a1, time), s2 = GetDistance(v2, a2, time);

        if (s1 > s2) 
            answer = 1;
        else if (s1 < s2) 
            answer = 2;
        else 
            answer = 0;

        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        answer = 1;
        while (GetDistance(v1, a1, answer) > GetDistance(v2, a2, answer)) 
            answer++;

        // end

        return answer;
    }
    #endregion

    #region Level 2
    public void FindMaxIndex(int[,] matrix, out int maxI, out int maxJ)
    {
        int max = matrix[0, 0];
        maxI = 0;
        maxJ = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxI = i;
                    maxJ = j;
                }
    }

    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
        int maxI_A, maxJ_A, maxI_B, maxJ_B;

        FindMaxIndex(A, out maxI_A, out maxJ_A);
        FindMaxIndex(B, out maxI_B, out maxJ_B);

        int temp = A[maxI_A, maxJ_A];
        A[maxI_A, maxJ_A] = B[maxI_B, maxJ_B];
        B[maxI_B, maxJ_B] = temp;

        // end
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }

    public int FindDiagonalMaxIndex(int[,] matrix)
    {
        int max = matrix[0, 0], ind = 0;

        for (int i = 1; i < matrix.GetLength(0); i++)
            if (matrix[i, i] > max)
            {
                max = matrix[i, i];
                ind = i;
            }
        return ind;
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        int[,] tempB = new int[4, 5];
        int[,] tempC = new int[5, 6];

        int maxB = FindDiagonalMaxIndex(B);
        int maxC = FindDiagonalMaxIndex(C);

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 5; j++)
            {
                if (i < maxB)
                    tempB[i, j] = B[i, j];
                else
                    tempB[i, j] = B[i + 1, j];
            }

        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 6; j++)
                if (i < maxC)
                    tempC[i, j] = C[i, j];
                else
                    tempC[i, j] = C[i + 1, j];

        B = tempB;
        C = tempC;
        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }

    public int FindMaxInColumn(int[,] matrix, int columnIndex, out int rowIndex)
    {
        int maxValue = int.MinValue;
        int n = matrix.GetLength(0);
        rowIndex = -1;

        for (int i = 0; i < n; i++)
            if (matrix[i, columnIndex] > maxValue)
            {
                maxValue = matrix[i, columnIndex];
                rowIndex = i;
            }

        return rowIndex;
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here
        int maxA, maxB;
        FindMaxInColumn(A, 0, out maxA);
        FindMaxInColumn(B, 0, out maxB);

        for (int j = 0; j < 6; j++)
        {
            int temp = A[maxA, j];
            A[maxA, j] = B[maxB, j];
            B[maxB, j] = temp;
        }

        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }
    
    public int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int answer  = 0, n = matrix.GetLength(1);

        for (int i = 0; i < n; i++)
            if (matrix[rowIndex, i] > 0)
                answer++;

        return answer;
    }
    public int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int answer = 0, n = matrix.GetLength(0);

        for (int i = 0; i < n; i++)
            if (matrix[i, colIndex] > 0)
                answer++;

        return answer;
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        int maxBVal = -1, maxBInd = -1;
        for (int i = 0; i < 4; i++)
            if (CountRowPositive(B, i) > maxBVal)
            {
                maxBVal = CountRowPositive(B, i);
                maxBInd = i;
            }

        int maxCVal = -1, maxCInd = -1;
        for (int j = 0; j < 6; j++)
            if (CountColumnPositive(C, j) > maxCVal)
            {
                maxCVal = CountColumnPositive(C, j);
                maxCInd = j;
            }

        int[,] resultB = new int[5, 5];

        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
            {
                if (i <= maxBInd)
                    resultB[i, j] = B[i, j];
                else if (i == maxBInd + 1)
                    resultB[i, j] = C[j, maxCInd];
                else
                    resultB[i, j] = B[i - 1, j];
            }

        B = resultB;
        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] SumPositiveElementsInColumns(int[,] matrix)
    {
        int[] array = new int[matrix.GetLength(1)];

        for (int j = 0; j < matrix.GetLength(1); j++)
            for (int i = 0; i < matrix.GetLength(0); i++)
                if (matrix[i, j] > 0) 
                    array[j] += matrix[i, j];

        return array;
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here
        answer = new int[A.GetLength(1) + C.GetLength(1)];

        int[] arrA = SumPositiveElementsInColumns(A), arrC = SumPositiveElementsInColumns(C);

        for (int i = 0; i < arrA.Length; i++)
            answer[i] = arrA[i];

        for (int i = arrA.Length; i < arrA.Length + arrC.Length; i++)
            answer[i] = arrC[i - arrA.Length];

        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here
        int maxI_A, maxJ_A, maxI_B, maxJ_B;

        FindMaxIndex(A, out maxI_A, out maxJ_A);
        FindMaxIndex(B, out maxI_B, out maxJ_B);

        int temp = A[maxI_A, maxJ_A];
        A[maxI_A, maxJ_A] = B[maxI_B, maxJ_B];
        B[maxI_B, maxJ_B] = temp;
        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }

    public int[,] RemoveRow(int[,] matrix, int rowIndex)
    {
        int[,] temp = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];

        for (int i = 0; i < temp.GetLength(0); i++)
            for (int j = 0; j < temp.GetLength(1); j++)
                if (i < rowIndex) temp[i, j] = matrix[i, j];
                else if (i >= rowIndex) 
                    temp[i, j] = matrix[i + 1, j];

        return temp;
    }
    public void Task_2_13(ref int[,] matrix)
    {
        // code here
        int max = matrix[0, 0], min = matrix[0, 0], mxRow = 0,  mnRow = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    mxRow = i;
                }

                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    mnRow = i;
                }
            }

        if (mnRow < mxRow)
        {
            matrix = RemoveRow(matrix, mxRow);
            matrix = RemoveRow(matrix, mnRow);
        }
        else if (mnRow > mxRow)
        {
            matrix = RemoveRow(matrix, mnRow);
            matrix = RemoveRow(matrix, mxRow);
        }
        else
            matrix = RemoveRow(matrix, mnRow);

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public double GetAverageWithoutMinMax(int[,] matrix)
    {
        int max = matrix[0, 0], indexMaxI = 0, indexMaxJ = 0, min = matrix[0, 0], indexMinI = 0, indexMinJ = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    indexMaxI = i;
                    indexMaxJ = j;
                }

                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    indexMinI = i;
                    indexMinJ = j;
                }
            }

        double sum = 0;
        int cnt = 0;

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if ((i == indexMaxI && j == indexMaxJ) || (i == indexMinI && j == indexMinJ)) continue;
                sum += matrix[i, j];
                cnt++;
            }

        return sum / cnt;
    }
    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        double[] array = { GetAverageWithoutMinMax(A), GetAverageWithoutMinMax(B), GetAverageWithoutMinMax(C) };

        if (array[0] < array[1] && array[1] < array[2]) 
            answer = 1;
        else if (array[0] > array[1] && array[1] > array[2]) 
            answer = -1;
        else 
            answer = 0;
        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }

    public int[,] SortRowsByMaxElement(int[,] matrix)
    {
        int[] maxs = new int[matrix.GetLength(0)];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            maxs[i] = matrix[i, 0];
            for (int j = 1; j < matrix.GetLength(1); j++)
                if (matrix[i, j] > maxs[i])
                    maxs[i] = matrix[i, j];
        }

        for (int i = 1; i < maxs.Length; i++)
        {
            int k = i - 1;
            while (k >= 0 && maxs[k] < maxs[k + 1])
            {
                (maxs[k], maxs[k + 1]) = (maxs[k + 1], maxs[k]);

                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    int temp = matrix[k, j];
                    matrix[k, j] = matrix[k + 1, j];
                    matrix[k + 1, j] = temp;

                }
                k--;
            }
        }
        return matrix;
    }
    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        A = SortRowsByMaxElement(A);
        B = SortRowsByMaxElement(B);

        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            bool flag = false;
            for (int j = 0; j < matrix.GetLength(1); j++)
                if (matrix[i, j] == 0)
                {
                    flag = true;
                    break;
                }
            if (!flag) continue;
            matrix = RemoveRow(matrix, i);
            i--;
        }
        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public int[] CreateArrayFromMins(int[,] matrix)
    {
        int[] arr = new int[matrix.GetLength(0)];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            arr[i] = matrix[i, i];
            for (int j = i + 1; j < matrix.GetLength(1); j++)
                if (matrix[i, j] < arr[i])
                    arr[i] = matrix[i, j];
        }
        return arr;
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }

    public void MatrixValuesChange(double[,] matrix)
    {
        double[] arr = new double[matrix.GetLength(0) * matrix.GetLength(1)];

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
                arr[i * matrix.GetLength(0) + j] = matrix[i, j];

        for (int i = 1; i < arr.Length; i++)
        {
            int j = i - 1;
            while (j >= 0 && arr[j] < arr[j + 1])
            {
                double temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
                j--;
            }
        }

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                bool flag = false;
                for (int k = 0; k < 5; k++)
                    if (matrix[i, j] == arr[k])
                    {
                        flag = true;
                        matrix[i, j] = (matrix[i, j] > 0) ? matrix[i, j] * 2 : matrix[i, j] / 2;
                        break;
                    }
                if (flag) continue;
                if (matrix[i, j] > 0)
                    matrix[i, j] /= 2;
                else
                    matrix[i, j] *= 2;
            }
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here
        MatrixValuesChange(A);
        MatrixValuesChange(B);
        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }

    public int FindRowWithMaxNegativeCount(int[,] matrix)
    {
        int max = int.MinValue, index = -1;

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int cnt = CountNegativeInRow(matrix, i);
            if (cnt > max)
            {
                max = cnt;
                index = i;
            }
        }
        return index;
    }
    public int CountNegativeInRow(int[,] matrix, int rowIndex)
    {
        int cnt = 0;

        for (int j = 0; j < matrix.GetLength(1); j++)
            if (matrix[rowIndex, j] < 0) 
                cnt++;

        return cnt;
    }
    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here
        maxA = FindRowWithMaxNegativeCount(A);
        maxB = FindRowWithMaxNegativeCount(B);
        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }


    public void FindRowMaxIndex(int[,] matrix, int rowIndex, out int columnIndex)
    {
        int max = matrix[rowIndex, 0];
        columnIndex = 0;

        for (int j = 0; j < matrix.GetLength(1); j++)
            if (matrix[rowIndex, j] > max)
            {
                max = matrix[rowIndex, j];
                columnIndex = j;
            }
    }
    public void ReplaceMaxElementOdd(int[,] matrix, int row, int column)
    {
        matrix[row, column] *= column + 1;
    }
    public void ReplaceMaxElementEven(int[,] matrix, int row, int column)
    {
        matrix[row, column] = 0;
    }
    public void ReplaceMatrixElements(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int index;
            if ((i + 1) % 2 != 0)
            {
                FindRowMaxIndex(matrix, i, out index);
                ReplaceMaxElementOdd(matrix, i, index);
            }
            else
            {
                FindRowMaxIndex(matrix, i, out index);
                ReplaceMaxElementEven(matrix, i, index);
            }
        }
    }
    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here
        ReplaceMatrixElements(A);
        ReplaceMatrixElements(B);
        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3

    public delegate double SumFunction(double x, int i, ref int result);
    public delegate double YFunction(double x);

    public double forSum1(double x, int i, ref int iFact)
    {
        if (i != 0) 
            iFact *= i;

        return Math.Cos(i * x) / iFact;
    }
    
    public double forSum2(double x, int i, ref int sign)
    {
        sign *= -1;
        return sign * Math.Cos(i * x) / (i * i);
    }

    public double Y_1(double x)
    {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }

    public double Y_2(double x)
    {
        return ((x * x) - Math.PI * Math.PI / 3) / 4;
    }

    public double Sum(SumFunction sumFunction, double x, int i)
    {
        double eps = 0.0001, sum = 0;
        int change = 1;
        double element = sumFunction(x, i, ref change);

        while (Math.Abs(element) > eps)
        {
            sum += element;
            element = sumFunction(x, ++i, ref change);
        }
        return sum;
    }
    public void GetSumAndY(SumFunction sFunction, YFunction yFunction, double a, double b, double h, double[,] SumAndY, int startI = 0)
    {
        int steps = (int)Math.Round((b - a) / h) + 1;
        for (int i = 0; i < steps; i++)
        {
            double x = a + i * h;

            double sum = Sum(sFunction, x, startI);
            double y = yFunction(x);

            SumAndY[i, 0] = sum;
            SumAndY[i, 1] = y;
        }
    }
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        firstSumAndY = new double[(int)((b1 - a1) / h1) + 1, 2];

        GetSumAndY(forSum1, Y_1, a1, b1, h1, firstSumAndY, 0);

        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];

        GetSumAndY(forSum2, Y_2, a2, b2, h2, secondSumAndY, 1);
        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }

    public delegate void SwapDirection(double[] array);
    public void SwapLeft(double[] array)
    {
        for (int i = 0; i < array.Length - 1; i += 2)
        {
            double temp = array[i];
            array[i] = array[i + 1];
            array[i + 1] = temp;
        }
    }

    public void SwapRight(double[] array)
    {
        for (int i = array.Length - 1; i > 0; i -= 2)
        {
            double temp = array[i];
            array[i] = array[i - 1];
            array[i - 1] = temp;
        }
    }

    public double GetSum(double[] array, int start, int h)
    {
        double sum = 0;

        for (int i = start; i < array.Length; i += h)
            sum += array[i];

        return sum;
    }

    public double Task_3_3(double[] array)
    {
        double answer = 0;
        SwapDirection swapper = default(SwapDirection);

        // code here
        double avg = 0;

        foreach (double item in array)
            avg += item;

        avg /= array.Length;

        if (array[0] > avg) 
            swapper = SwapLeft;
        else 
            swapper = SwapRight;

        swapper(array);
        answer = GetSum(array, 1, 2);
        // end

        return answer;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public int CountSignFlips(YFunction YFunction, double a, double b, double h)
    {
        double lastY = YFunction(a);
        int cnt = 0;
        for (double i = a + h; i <= b; i += h)
        {
            double curY = YFunction(i);
            if (curY > 0 && lastY < 0 || curY < 0 && lastY > 0) 
                cnt++;
            lastY = curY;
        }
        cnt++;

        return cnt;
    }

    public double Func1(double x)
    {
        return x * x - Math.Sin(x);
    }

    public double Func2(double x)
    {
        return Math.Exp(x) - 1;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here
        func1 = CountSignFlips(Func1, 0, 2, 0.1);
        func2 = CountSignFlips(Func2, -1, 1, 0.2);
        // end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public delegate int CountPositive(int[,] matrix, int index);

    public void InsertColumn(ref int[,] matrixB, int CountRow, int[,] matrixC, int CountColumn)
    {
        int[,] temp = new int[matrixB.GetLength(0) + 1, matrixC.GetLength(0)];

        for (int i = 0; i < temp.GetLength(0); i++)
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                if (i < CountRow) 
                    temp[i, j] = matrixB[i, j];
                else if (i == CountRow) 
                    temp[i, j] = matrixC[j, CountColumn];
                else 
                    temp[i, j] = matrixB[i - 1, j];
            }
        matrixB = temp;
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        CountPositive countPositive;

        int max = int.MinValue, indColumn = -1;
        countPositive = CountRowPositive;

        for (int i = 0; i < B.GetLength(0); i++)
        {
            int cnt = countPositive(B, i);
            if (cnt > max)
            {
                max = cnt;
                indColumn = i;
            }
        }

        max = int.MinValue;
        int indRow = -1;

        countPositive = CountColumnPositive;

        for (int j = 0; j < C.GetLength(1); j++)
        {
            int cnt = countPositive(C, j);
            if (cnt > max)
            {
                max = cnt;
                indRow = j;
            }
        }
        InsertColumn(ref B, indRow + 1, C, indColumn);
        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }

    public delegate void FindElementDelegate(int[,] matrix, out int foundI, out int foundJ);

    public void RemoveRows(ref int[,] matrix, FindElementDelegate findElementDelegate1, FindElementDelegate findElementDelegate2)
    {
        int i1, j1, i2, j2;

        findElementDelegate1(matrix, out i1, out j1);
        findElementDelegate2(matrix, out i2, out j2);

        matrix = RemoveRow(matrix, i1);

        if (i2 < i1)
            matrix = RemoveRow(matrix, i2);
        else if (i2 > i1)
            matrix = RemoveRow(matrix, i2 - 1);
    }

    public void FindMinIndex(int[,] matrix, out int minI, out int minJ)
    {
        minI = -1;
        minJ = -1;
        int minValue = int.MaxValue;

        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < minValue)
                {
                    minValue = matrix[i, j];
                    minI = i;
                    minJ = j;
                }
            }
    }
    public void Task_3_13(ref int[,] matrix)
    {
        // code here
        RemoveRows(ref matrix, FindMaxIndex, FindMinIndex);
        // end
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }

    public delegate void ReplaceMaxElement(int[,] matrix, int i, int index);

    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElementOdd, ReplaceMaxElement replaceMaxElementEven)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int index;
            if ((i + 1) % 2 != 0)
            {
                FindRowMaxIndex(matrix, i, out index);
                replaceMaxElementOdd(matrix, i, index);
            }
            else
            {
                FindRowMaxIndex(matrix, i, out index);
                replaceMaxElementEven(matrix, i, index);
            }
        }
    }
    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion

    #region bonus part

    public delegate void MatrixConverter(double[,] matrix);
    public void ToUpperTriangular(double[,] matrix)
    {
        for (int j = 0; j < matrix.GetLength(1); j++)
            for (int i = j + 1; i < matrix.GetLength(0); i++)
            {
                double k = - (matrix[i, j] / matrix[j, j]);
                matrix[i, j] = 0;
                for (int n = j + 1; n < matrix.GetLength(1); n++)
                    matrix[i, n] += matrix[j, n] * k;
            }
    }
    public void ToLowerTriangular(double[,] matrix)
    {
        for (int j = matrix.GetLength(1) - 1; j >= 0; j--)
            for (int i = j - 1; i >= 0; i--)
            {
                double k = - (matrix[i, j] / matrix[j, j]);
                matrix[i, j] = 0;
                for (int n = j - 1; n >= 0; n--)
                    matrix[i, n] += matrix[j, n] * k;
            }
    }
    public void ToLeftDiagonal(double[,] matrix)
    {
        ToUpperTriangular(matrix);
        ToLowerTriangular(matrix);
    }

    public void ToRightDiagonal(double[,] matrix)
    {
        ToLowerTriangular(matrix);
        ToUpperTriangular(matrix);
    }
    public double[,] Task_4(double[,] matrix, int index)
    {
        // code here
        MatrixConverter[] mc = new MatrixConverter[]{ToUpperTriangular, ToLowerTriangular, ToLeftDiagonal, ToRightDiagonal};

        mc[index](matrix);
        // end

        return matrix;
    }
    #endregion
}
