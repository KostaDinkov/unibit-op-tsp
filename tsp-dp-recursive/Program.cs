using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Text;


namespace GasStationsDemo
{

    public class Program
    {
        public static void Main()
        {
            int[,] distanceMatrix = Utils.TestData.Test3;
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
            BigInteger cycles = (N - 1) * (N - 2) * BigInteger.Pow(2, N - 3) + (N - 1);
            BigInteger cyclePerSecond = (cycles / (stopwatch.ElapsedMilliseconds)) * 1000;
            Console.WriteLine($"Total Tour Computations = {cycles}. Computations per second: {cyclePerSecond}");
        }
    }


    public class TspDp
    {
        private readonly int[,] distance;
        private readonly int[,] memo;
        private readonly int nodeCount;
        private readonly int start;

        public TspDp(int[,] distance) : this(0, distance)
        {
        }

        public TspDp(int start, int[,] distance)
        {
            nodeCount = distance.GetLength(0);
            TourCost = ushort.MaxValue;
            this.start = start;
            this.distance = distance;
            memo = new int[nodeCount, 1 << nodeCount];
            ComputeMinTourCost();
        }


        public int TourCost { get; private set; }

        public string ReconstructTour()
        {
            var tour = new List<int>();
            var lastNode = start;
            var state = (1 << nodeCount) - 1;
            tour.Add(start);


            while (state > 1)
            {
                var node = -1;
                for (var j = 0; j < nodeCount; j++)
                {
                    if (j == start || NotIn(j, state)) continue;
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

            tour.Add(start);

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
            var finalState = (1 << nodeCount) - 1;



            for (var end = 0; end < nodeCount; end++)
            {
                if (end == start) continue;
                memo[end, (1 << start) | (1 << end)] = distance[start, end];
            }


            for (var r = 3; r <= nodeCount; r++)
            {
                foreach (var subset in Combinations(r, nodeCount))
                {
                    if (NotIn(start, subset)) continue;
                    for (var next = 0; next < nodeCount; next++)
                    {
                        if (next == start || NotIn(next, subset)) continue;
                        var subsetWithoutNext = subset ^ (1 << next);
                        var minDist = int.MaxValue;
                        for (var end = 0; end < nodeCount; end++)
                        {
                            if (end == start || end == next || NotIn(end, subset)) continue;
                            var newDistance = memo[end, subsetWithoutNext] + distance[end, next];
                            if (newDistance < minDist)
                            {
                                minDist = newDistance;
                            }
                        }
                        memo[next, subset] = minDist;


                        // This works too but is almost twice slower
                        //
                        //memo[next, subset] = Enumerable.Range(0, nodeCount)
                        //    .Where(node => node != start && node != next && !NotIn(node, subset))
                        //    .Select(node => memo[node, subsetWithoutNext] + distance[node, next])
                        //    .Min();

                    }

                }
            }


            for (var i = 0; i < nodeCount; i++)
            {
                if (i == start) continue;
                var tourCost = memo[i, finalState] + distance[i, start];
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