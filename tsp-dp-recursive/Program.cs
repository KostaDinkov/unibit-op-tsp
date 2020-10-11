using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Text;
using Utils;

namespace GasStationsDemo
{
    /*
     * The code below is mostly the same as the code in the GasStations project
     * There are some additional information being printed on the console like
     * the memo table (if n <= 5),
     * the tour itself, e.g. 0-3-4-1-5,
     * and some benchmarking timings
     * The input is not read from the console, but being imported
     * from a class with some aditional test data.
     * 
     * author Kosta Dinkov
     * 
     */
    public class Program
    {
        public static void Main()
        {
            // The TestData class in the Utils project contains additional
            // test cases besides the one from the project specification
            // See Utils.TestData for additional info
            var distanceMatrix = TestData.Test1;
            var N = distanceMatrix.GetLength(0);
            var stopwatch = new Stopwatch();
            stopwatch.Start();
            var startNode = 0;
            var solver = new TspDp(startNode, distanceMatrix);

            if (N <= 5)
            {
                solver.PrintMemoTbl();
            }

            stopwatch.Stop();
            Console.WriteLine("Tour: " + solver.ReconstructTour());
            Console.WriteLine("Tour cost: " + solver.TourCost);
            Console.WriteLine($"Execution time for N={N}: {stopwatch.ElapsedMilliseconds / 1000d} seconds");
            var cycles = (N - 1) * (N - 2) * BigInteger.Pow(2, N - 3) + (N - 1);
            var cyclePerSecond = cycles / stopwatch.ElapsedMilliseconds * 1000;
            Console.WriteLine($"Total Tour Computations = {cycles}. Computations per second: {cyclePerSecond}");
        }
    }
    
    public class TspDp
    {
        private readonly int[,] distance;
        private readonly int[,] memo;
        private readonly int nodeCount;
        private readonly int startNode;

        public TspDp(int[,] distance) : this(0, distance)
        {
        }

        public TspDp(int startNode, int[,] distance)
        {
            nodeCount = distance.GetLength(0);
            TourCost = ushort.MaxValue;
            this.startNode = startNode;
            this.distance = distance;
            memo = new int[nodeCount, 1 << nodeCount];
            ComputeMinTourCost();
        }


        public int TourCost { get; private set; }

        public string ReconstructTour()
        {
            var tour = new List<int>();
            var lastNode = startNode;
            var state = (1 << nodeCount) - 1;
            tour.Add(startNode);

            while (state > 1)
            {
                var node = -1;
                for (var j = 0; j < nodeCount; j++)
                {
                    if (j == startNode || NotIn(j, state)) continue;
                    if (node == -1) node = j;
                    var prevDist = memo[node, state] + distance[node, lastNode];
                    var newDist = memo[j, state] + distance[j, lastNode];
                    if (newDist < prevDist)
                    {
                        node = j;
                    }
                }

                tour.Add(node);
                state = state ^ (1 << node);
                lastNode = node;
            }

            tour.Add(startNode);
            tour.Reverse();
            return string.Join(" => ", tour);
        }


        private void ComputeMinTourCost()
        {
            if (nodeCount == 2)
            {
                TourCost = distance[1, 0] * 2;
                return;
            }

            var fullSet = (1 << nodeCount) - 1;

            for (var nextNode = 0; nextNode < nodeCount; nextNode++)
            {
                if (nextNode == startNode) continue;
                memo[nextNode, (1 << startNode) | (1 << nextNode)] = distance[startNode, nextNode];
            }

            for (var s = 3; s <= nodeCount; s++)
            {
                foreach (var subset in Combinations(s, nodeCount))
                {
                    if (NotIn(startNode, subset)) continue;
                    for (var nextEnd = 0; nextEnd < nodeCount; nextEnd++)
                    {
                        if (nextEnd == startNode || NotIn(nextEnd, subset)) continue;
                        var subMinusNextEnd = subset ^ (1 << nextEnd);
                        var minDist = int.MaxValue;
                        for (var lastEnd = 0; lastEnd < nodeCount; lastEnd++)
                        {
                            if (lastEnd == startNode || lastEnd == nextEnd || NotIn(lastEnd, subset)) continue;
                            var newDistance = memo[lastEnd, subMinusNextEnd] + distance[lastEnd, nextEnd];
                            if (newDistance < minDist)
                            {
                                minDist = newDistance;
                            }
                        }

                        memo[nextEnd, subset] = minDist;

                        // This works too but is almost twice slower
                        //
                        //memo[next, subset] = Enumerable.Range(0, nodeCount)
                        //    .Where(node => node != startNode && node != next && !NotIn(node, subset))
                        //    .Select(node => memo[node, subsetWithoutNext] + distance[node, next])
                        //    .Min();
                    }
                }
            }

            for (var node = 0; node < nodeCount; node++)
            {
                if (node == startNode) continue;
                var tourCost = memo[node, fullSet] + distance[node, startNode];
                if (tourCost < TourCost)
                {
                    TourCost = tourCost;
                }
            }
        }
        
        public static bool NotIn(int node, int subset)
        {
            return ((1 << node) & subset) == 0;
        }

        public static List<int> Combinations(int r, int n)
        {
            var subsets = new List<int>();
            Combinations(0, 0, r, n, subsets);
            return subsets;
        }

        private static void Combinations(int set, int at, int r, int n, List<int> subsets)
        {
            var elementsLeftToPick = n - at;
            if (elementsLeftToPick < r) return;

            if (r == 0)
            {
                subsets.Add(set);
            }
            else
            {
                for (var i = at; i < n; i++)
                {
                    set ^= 1 << i;

                    Combinations(set, i + 1, r - 1, n, subsets);

                    set ^= 1 << i;
                }
            }
        }

        public void PrintMemoTbl()
        {
            Console.OutputEncoding = Encoding.UTF8;
            var result = new StringBuilder();
            Console.Write(" ".PadRight(3));

            for (var h = 3; h < memo.GetLength(1); h++)
            {
                if (h % 2 == 0) continue;
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.Write(Convert.ToString(h).PadLeft(4));
                Console.ResetColor();
                Console.Write("\u2502".PadLeft(3));
            }

            Console.Write("\n" + new string('\u2501', memo.GetLength(1) / 2 * 7 - 4) + "\n");

            for (var row = 1; row < memo.GetLength(0); row++)
            {
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.Write(row.ToString().PadRight(3));
                Console.ResetColor();
                for (var col = 3; col < memo.GetLength(1); col++)
                {
                    if (col % 2 == 0) continue;
                    if (memo[row, col] != 0) Console.ForegroundColor = ConsoleColor.Cyan;
                    Console.Write(memo[row, col].ToString().PadLeft(4));
                    Console.ResetColor();
                    Console.Write("\u2502".PadLeft(3));
                }

                Console.WriteLine();
            }
        }
    }
}