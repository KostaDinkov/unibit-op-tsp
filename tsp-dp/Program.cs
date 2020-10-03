using System;
using System.Collections.Generic;
using System.Text;

namespace tsp_dp
{
    /**
     * An implementation of the traveling salesman problem in Java using dynamic programming to improve
     * the time complexity from O(n!) to O(n^2 * 2^n).
     * <p>
     *     Time Complexity: O(n^2 * 2^n) Space Complexity: O(n * 2^n)
     *     @author William Fiset, william.alexandre.fiset@gmail.com
     */
    public class Program
    {
        public static void Main()
        {
            double[,] distanceMatrix = TestData.Test5;
            

            var startNode = 0;
            var solver = new TspDp(startNode, distanceMatrix);
            //solver.PrintMemoTbl();
            Console.WriteLine("Tour: " + solver.ReconstructTour());
            Console.WriteLine("Tour cost: " + solver.TourCost);
        }
    }


    public class TspDp
    {
        private readonly double[,] distance;
        private readonly double[,] memo;
        private readonly int NodeCount;
        private readonly int start;

        public TspDp(double[,] distance) : this(0, distance)
        {
        }

        public TspDp(int start, double[,] distance)
        {
            NodeCount = distance.GetLength(0);
            TourCost = double.PositiveInfinity;
            this.start = start;
            this.distance = distance;
            memo = new double[NodeCount, 1 << NodeCount];
            ComputeMinTourCost();
        }

        // Returns the minimal tour cost.
        public double TourCost { get; private set; }

        // Returns the optimal tour for the traveling salesman problem.
        public string ReconstructTour()
        {
            var tour = new List<int>();
            var lastNode = start;
            var state = (1 << NodeCount) - 1;
            tour.Add(start);

            // Reconstruct TSP path from memo table.
            while (state > 1)
            {
                // the current number of the gas station to be added to the tour
                var node = -1;
                for (var j = 0; j < NodeCount; j++)
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

        // Solves the traveling salesman problem and caches solution.
        private void ComputeMinTourCost()
        {
            if (NodeCount == 2)
            {
                TourCost = distance[1, 0] * 2;
                return;
            }

            ;

            var finalState = (1 << NodeCount) - 1;


            // Add all outgoing edges from the starting node to memo table.
            for (var end = 0; end < NodeCount; end++)
            {
                if (end == start) continue;
                memo[end, (1 << start) | (1 << end)] = distance[start, end];
            }


            for (var r = 3; r <= NodeCount; r++)
            {
                foreach (var subset in Combinations(r, NodeCount))
                {
                    if (NotIn(start, subset)) continue;
                    for (var next = 0; next < NodeCount; next++)
                    {
                        if (next == start || NotIn(next, subset)) continue;
                        var subsetWithoutNext = subset ^ (1 << next);
                        var minDist = double.PositiveInfinity;
                        for (var end = 0; end < NodeCount; end++)
                        {
                            if (end == start || end == next || NotIn(end, subset)) continue;
                            var newDistance = memo[end, subsetWithoutNext] + distance[end, next];
                            if (newDistance < minDist)
                            {
                                minDist = newDistance;
                            }
                        }

                        memo[next, subset] = minDist;
                    }
                }
            }

            // Connect tour back to starting node and minimize cost.
            for (var i = 0; i < NodeCount; i++)
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

        // This method generates all bit sets of size n where r bits
        // are set to one. The result is returned as a list of integer masks.
        public static List<int> Combinations(int r, int n)
        {
            var subsets = new List<int>();
            Combinations(0, 0, r, n, subsets);
            return subsets;
        }

        // To find all the combinations of size r we need to recurse until we have
        // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
        // an element which is found after the position of our last selected element
        private static void Combinations(int set, int at, int r, int n, List<int> subsets)
        {
            // Return early if there are more elements left to select than what is available.
            var elementsLeftToPick = n - at;
            if (elementsLeftToPick < r) return;

            // We selected 'r' elements so we found a valid subset!
            if (r == 0)
            {
                subsets.Add(set);
            }
            else
            {
                for (var i = at; i < n; i++)
                {
                    // Try including this element
                    set ^= 1 << i;

                    Combinations(set, i + 1, r - 1, n, subsets);

                    // Backtrack and try the instance where we did not include this element
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