﻿using System;

namespace Utils
{
    public class Matrix
    {
        public static void Print(int[,] matrix)
        {
            for (int row = 0; row < matrix.GetLength(0); row++)
            {
                for (int col = 0; col < matrix.GetLength(1); col++)
                {
                    Console.Write(matrix[row, col].ToString().PadRight(3));
                }
                Console.WriteLine();
            }
        }
    }
}
