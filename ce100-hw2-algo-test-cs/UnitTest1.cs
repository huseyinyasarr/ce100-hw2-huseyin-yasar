using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using ce100_hw2_algo_lib_cs;

namespace ce100_hw2_algo_test_cs

/**
* @file CE100-hw2-algo-test-cs
* @author huseyin yasar
* @date 31 March 2022
*
* @brief <b> HW-2 Unit Tests </b>
*
* HW-2 Sample Unit Test
*
* @see http://bilgisayar.mmf.erdogan.edu.tr/en/
*
*/


{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void HeapSortTestMethod1()
        {
            int[] array = new int[] { 1, 5, 3, 2, 0 };
            int n = array.Length;
            int[] expected = new int[] { 0,1,2,3,5};

            Class1.heapSort(array,n);
            CollectionAssert.AreEqual(expected,array);
        }

        [TestMethod]
        public void HeapSortTestMethod2()
        {
            int[] array = new int[] { 55 ,25 ,89 ,34 ,12 ,19 ,78 ,95 ,1 ,100 };
            int n = array.Length;
            int[] expected = new int[] { 1 ,12 ,19 ,25 ,34 ,55 ,78 ,89 ,95 ,100 };

            Class1.heapSort(array, n);
            CollectionAssert.AreEqual(expected, array);
        }
        [TestMethod]
        public void HeapSortTestMethod3()
        {
            int[] array = new int[] { -486, 145, -23, 62, 0, 2445, 354};
            int n = array.Length;
            int[] expected = new int[] { -486, -23, 0, 62, 145, 354, 2445};

            Class1.heapSort(array, n);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void RadixSortTestMethod1()
        {
            int[] array = new int[] { 4,8,1,3,7,6,2,9,5,0};
            int n = array.Length;
            int[] expected = new int[] { 0,1,2,3,4,5,6,7,8,9};

            Class1.radixsort(array, n);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void RadixSortTestMethod2()
        {
            int[] array = new int[] { 170, 45, 75, 90, 802, 24, 2, 66 };
            int n = array.Length;
            int[] expected = new int[] { 2, 24, 45, 66, 75, 90, 170, 802 };

            Class1.radixsort(array, n);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void RadixSortTestMethod3()
        {
            int[] array = new int[] { 4501,123,865,423,1532,901,4,56};
            int n = array.Length;
            int[] expected = new int[] { 4,56,123,423,865,901,1532,4501};

            Class1.radixsort(array, n);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void CountingSortTestMethod1()
        {
            int[] array = new int[] { 4501, 123, 865, 423, 1532, 901, 4, 56 };
            int[] expected = new int[] { 4, 56, 123, 423, 865, 901, 1532, 4501 };

            Class1.countingsort(array);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void CountingSortTestMethod2()
        {
            int[] array = new int[] { 55, 25, 89, 34, 12, 19, 78, 95, 1, 100 };
            int[] expected = new int[] { 1, 12, 19, 25, 34, 55, 78, 89, 95, 100 };

            Class1.countingsort(array);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void CountingSortTestMethod3()
        {
            int[] array = new int[] { 486,123,458,2,145,44,33,852,612,952};
            int[] expected = new int[] { 2,33,44,123,145,458,486,612,852,952};

            Class1.countingsort(array);
            CollectionAssert.AreEqual(expected, array);
        }



        [TestMethod]
        public void QuickSortLomutoTestMethod1()
        {
            int[] array = new int[] { 10, 7, 8, 9, 1, 5 };
            int[] expected = new int[] { 1,5,7,8,9,10 };
            int n = array.Length;
            Class1.quickSortLomuto(array,0,n-1);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void QuickSortLomutoTestMethod2()
        {
            int[] array = new int[] { 42,15,73,28,95,39,1,61};
            int[] expected = new int[] { 1,15,28,39,42,61,73,95};
            int n = array.Length;
            Class1.quickSortLomuto(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void QuickSortLomutoTestMethod3()
        {
            int[] array = new int[] { 170, 45, 75, 90, 802, 24, 2, 66 };
            int[] expected = new int[] { 2, 24, 45, 66, 75, 90, 170, 802 };
            int n = array.Length;
            Class1.quickSortLomuto(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }



        [TestMethod]
        public void QuickSortHoaresTestMethod1()
        {
            int[] array = new int[] { 10, 7, 8, 9, 1, 5 };
            int[] expected = new int[] { 1, 5, 7, 8, 9, 10 };
            int n = array.Length;
            Class1.quickSortHoares(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }


        [TestMethod]
        public void QuickSortHoaresTestMethod2()
        {
            int[] array = new int[] { 42, 15, 73, 28, 95, 39, 1, 61 };
            int[] expected = new int[] { 1, 15, 28, 39, 42, 61, 73, 95 };
            int n = array.Length;
            Class1.quickSortHoares(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }


        [TestMethod]
        public void QuickSortHoaresTestMethod3()
        {
            int[] array = new int[] { 170, 45, 75, 90, 802, 24, 2, 66 };
            int[] expected = new int[] { 2, 24, 45, 66, 75, 90, 170, 802 };
            int n = array.Length;
            Class1.quickSortHoares(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }


        [TestMethod]
        public void multiplyMatrixRecTestMethod1()
        {
            int row1 = 3, col1 = 3,
            row2 = 3, col2 = 3;
            int[,] A = {{1, 2, 3},
                    {4, 5, 6},
                    {7, 8, 9}};

            int[,] B = {{1, 2, 3},
                    {4, 5, 6},
                    {7, 8, 9}};

            int[,] C =new int[row1,col2];

            int[,] expected = { { 30, 36, 42 }, { 66, 81, 96 }, { 102, 126, 150 } };            
            Class1.multiplyMatrixRec(row1,col1,A,row2,col2,B,C);
            CollectionAssert.AreEqual(expected,C);
        }



        [TestMethod]
        public void multiplyMatrixRecTestMethod2()
        {
            int row1 = 3, col1 = 3,
            row2 = 3, col2 = 3;
            int[,] A = {{10, 11, 12},
                       {13, 14, 15},
                       {16, 17, 18}};
            Class1.i = 0; Class1.j = 0; Class1.k = 0;
            int[,] B = {{10, 11, 12},
                       {13, 14, 15},
                       {16, 17, 18}};

            int[,] C = new int[row1, col2];

            int[,] expected = { { 435, 468, 501 }, { 552, 594, 636 }, { 669, 720, 771 } };
            Class1.multiplyMatrixRec(row1, col1, A, row2, col2, B, C);
            CollectionAssert.AreEqual(expected, C);
        }


        [TestMethod]
        public void multiplyMatrixRecTestMethod3()
        {
            int row1 = 2, col1 = 2,
            row2 = 2, col2 = 2;
            int[,] A = {{12,32},
                    {41, 3}};
            Class1.i = 0; Class1.j = 0; Class1.k = 0;

            int[,] B = {{23,41},
                    {9,12}};

            int[,] C = new int[row1,col2];

            int[,] expected = {{ 564, 876}, { 970, 1717} };
            Class1.multiplyMatrixRec(row1, col1, A, row2, col2, B, C);
            CollectionAssert.AreEqual(expected, C);
        }





        [TestMethod]
        public void multiplyMatrixIterativeTestMethod1()
        {
            
            int[,] A = {{ 1, 1, 1, 1 },
                         { 2, 2, 2, 2 },
                         { 3, 3, 3, 3 },
                         { 4, 4, 4, 4 }};
            Class1.i = 0; Class1.j = 0; Class1.k = 0;
            int N = 4;
            int[,] B = { {1, 1, 1, 1 },
                         { 2, 2, 2, 2 },
                         { 3, 3, 3, 3 },
                         { 4, 4, 4, 4 }};

            int[,] C = new int[N, N];

            int[,] expected = { { 10,10,10,10 },
                                { 20,20,20,20},
                                { 30,30,30,30},
                                { 40,40,40,40} };
            Class1.multiply(A,B,C);
            CollectionAssert.AreEqual(expected, C);
        }



        [TestMethod]
        public void multiplyMatrixIterativeTestMethod2()
        {

            int[,] A = {{ 1, 2, 3, 4 },
                         { 5, 6, 7, 8 },
                         { 1, 2, 3, 4 },
                         { 5, 6, 7, 8 }};
            Class1.i = 0; Class1.j = 0; Class1.k = 0;
            int N = 4;
            int[,] B = { { 1, 2, 3, 4 },
                         { 5, 6, 7, 8 },
                         { 1, 2, 3, 4 },
                         { 5, 6, 7, 8 }};

            int[,] C = new int[N, N];

            int[,] expected = { { 34,44,54,64 },
                                { 82,108,134,160},
                                { 34,44,54,64 },
                                { 82,108,134,160} };
            Class1.multiply(A, B, C);
            CollectionAssert.AreEqual(expected, C);
        }



        [TestMethod]
        public void multiplyMatrixIterativeTestMethod3()
        {

            int[,] A = {{ 10,11,12,13},
                         { 13,14,15,16},
                         { 16,17,18,19},
                         { 20,21,22,23} };
            Class1.i = 0; Class1.j = 0; Class1.k = 0;
            int N = 4;
            int[,] B = { { 10,11,12,13},
                         { 13,14,15,16},
                         { 16,17,18,19},
                         {20,21,22,23 } };

            int[,] C = new int[N, N];

            int[,] expected = { { 695,741,787,833 },
                                { 872,930,988,1046},
                                { 1049,1119,1189,1259 },
                                { 1285,1371,1457,1543} };
            Class1.multiply(A, B, C);
            CollectionAssert.AreEqual(expected, C);
        }




        [TestMethod]
        public void MatrixMulStrTest1()
        {
            float[,] array = { { 1, 1 },
                           { 2, 2 } };
            

            float[,] array2 = { { 3, 3 },
                            { 4, 4 } };

            float[,] result = { { 7, 7 },
                            { 14, 14 } };

            float[,] exp = Class1.strassen(array, array2, 2);
            CollectionAssert.AreEqual(exp, result);
        }



        [TestMethod]
        public void MatrixMulStrTest2()
        {
            float[,] array = { { 4, 3 },
                           { 1, 5 } };
            

            float[,] array2 = { { 6, 2 },
                            { 7, 0 } };

            float[,] result = { { 45, 8 },
                            { 41, 2 } };

            float[,] exp = Class1.strassen(array, array2, 2);
            CollectionAssert.AreEqual(exp, result);
        }

        [TestMethod]
        public void MatrixMulStrTest3()
        {
            float[,] array = { { 5, 3 },
                           { 9, 6 } };
            

            float[,] array2 = { { 1, 0 },
                            { 5, 7 } };

            float[,] result = { { 20, 21 },
                            { 39, 42 } };

            float[,] exp = Class1.strassen(array, array2, 2);
            CollectionAssert.AreEqual(exp, result);
        }


        [TestMethod]
        public void RandomizedQuickSortLomutoTestMethod1()
        {
            int[] array = new int[] { 170, 45, 75, 90, 802, 24, 2, 66 };
            int[] expected = new int[] { 2, 24, 45, 66, 75, 90, 170, 802 };
            int n = array.Length;
            Class1.randomQuickSortLomuto(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }

        [TestMethod]
        public void RandomizedQuickSortLomutoTestMethod2()
        {
            int[] array = new int[] { 42, 15, 73, 28, 95, 39, 1, 61 };
            int[] expected = new int[] { 1, 15, 28, 39, 42, 61, 73, 95 };
            int n = array.Length;
            Class1.randomQuickSortLomuto(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }


        [TestMethod]
        public void RandomizedQuickSortLomutoTestMethod3()
        {
            int[] array = new int[] { 10, 7, 8, 9, 1, 5 };
            int[] expected = new int[] { 1,5,7,8,9,10 };
            int n = array.Length;
            Class1.randomQuickSortLomuto(array, 0, n - 1);
            CollectionAssert.AreEqual(expected, array);
        }


        [TestMethod]
        public void RandomizedQuickSortHoaresTestMethod1()
        {
            int[] array = new int[] { 10, 7, 8, 9, 1, 5 };
            int[] expected = new int[] { 1, 5, 7, 8, 9, 10 };
            int n = array.Length;
            Class1.randomQuicksortHoare(array,0,n-1);
            CollectionAssert.AreEqual(array, expected);
        }


        [TestMethod]
        public void RandomizedQuickSortHoaresTestMethod2()
        {
            int[] array = new int[] { 42, 15, 73, 28, 95, 39, 1, 61 };
            int[] expected = new int[] { 1, 15, 28, 39, 42, 61, 73, 95 };
            int n = array.Length;
            Class1.randomQuicksortHoare(array, 0, n - 1);
            CollectionAssert.AreEqual(array, expected);
        }


        [TestMethod]
        public void RandomizedQuickSortHoaresTestMethod3()
        {
            int[] array = new int[] { 486, 123, 458, 2, 145, 44, 33, 852, 612, 952 };
            int[] expected = new int[] { 2, 33, 44, 123, 145, 458, 486, 612, 852, 952 };
            int n = array.Length;
            Class1.randomQuicksortHoare(array, 0, n - 1);
            CollectionAssert.AreEqual(array, expected);
        }


        [TestMethod]
        public void PrioritryHeapSortTestMethod1()
        {
            
            int[] arr = new int[] { 45,12,78,27,63};
            Class1.insert(78);
            Class1.insert(63);
            Class1.insert(45);
            Class1.insert(27);
            Class1.insert(12);
            int expected = Class1.extractMax();
            Assert.AreEqual(78, expected);
        }


        [TestMethod]
        public void PrioritryHeapSortTestMethod2()
        {

            int[] arr = new int[] { 75,14,51,9,86,67,44,32,23,99};
            Class1.insert(99);
            Class1.insert(86);
            Class1.insert(75);
            Class1.insert(67);
            Class1.insert(51);
            Class1.insert(44);
            Class1.insert(32);
            Class1.insert(23);
            Class1.insert(14);
            Class1.insert(9);
            int expected = Class1.extractMax();
            Assert.AreEqual(99, expected);
        }


        [TestMethod]
        public void PrioritryHeapSortTestMethod3()
        {

            int[] arr = new int[] { 345,178,546,83,186,943,782,623};
            Class1.insert(943);
            Class1.insert(782);
            Class1.insert(623);
            Class1.insert(546);
            Class1.insert(345);
            Class1.insert(186);
            Class1.insert(178);
            Class1.insert(83);
            int expected = Class1.extractMax();
            Assert.AreEqual(943, expected);
        }



        [TestMethod]
        public void SelectionProblemTest1()
        {
            int[] arr = { 12, 3, 5, 7, 4, 19, 26 };
            int n = arr.Length, k = 3;
            int result = Class1.SelectionProblem(arr, 0, n - 1, k);
            Assert.AreEqual(5, result);
        }


        [TestMethod]
        public void SelectionProblemTest2()
        {
            int[] arr = { 41,32,84,12,25,65};
            int n = arr.Length, k = 2;
            int result = Class1.SelectionProblem(arr, 0, n - 1, k);
            Assert.AreEqual(25, result);
        }


        [TestMethod]
        public void SelectionProblemTest()
        {
            int[] arr = { 254,433,842,742,965,821,539};
            int n = arr.Length, k = 7;
            int result = Class1.SelectionProblem(arr, 0, n - 1, k);
            Assert.AreEqual(965, result);
        }
    }
}
