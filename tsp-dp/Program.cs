using System;
using System.Collections.Generic;
using Console = Colorful.Console;
using System.Drawing;
using System.Text;

namespace tsp_dp
{
    /**
 * An implementation of the traveling salesman problem in Java using dynamic programming to improve
 * the time complexity from O(n!) to O(n^2 * 2^n).
 *
 * <p>Time Complexity: O(n^2 * 2^n) Space Complexity: O(n * 2^n)
 *
 * @author William Fiset, william.alexandre.fiset@gmail.com
 */
    public class Program
    {
        public static void Main()
        {
            // Create adjacency matrix
            //int n = 6;
            //double[,] distanceMatrix = new double[n, n];

            //for (int i = 0; i < n; i++)
            //{
            //    for (int j = 0; j < n; j++)
            //    {
            //        distanceMatrix[i, j] = 10000;
            //    }
            //}
            //distanceMatrix[5, 0] = 10;
            //distanceMatrix[1, 5] = 12;
            //distanceMatrix[4, 1] = 2;
            //distanceMatrix[2, 4] = 4;
            //distanceMatrix[3, 2] = 6;
            //distanceMatrix[0, 3] = 8;

            double[,] distanceMatrix = new double[,]
            {
                { 0, 1, 4, 5, 10},
                { 1, 0, 2, 6, 3 },
                { 4, 2, 0, 7, 8 },
                { 5, 6, 7, 0, 9 },
                {10, 3, 8, 9, 0 }
            };

            int startNode = 0;
            TspDp solver =
                new TspDp(startNode, distanceMatrix);

            Console.WriteLine("Tour: " + solver.getTour());
            Console.WriteLine("Tour cost: " + solver.getTourCost());
        }
    }


    public class TspDp
    {

        private readonly int N;
        private readonly int start;
        private readonly double[,] distance;
        private List<int> tour = new List<int>();
        private double minTourCost = Double.PositiveInfinity;
        private bool ranSolver = false;


        public TspDp(double[,] distance) : this(0, distance)
        {

        }

        public TspDp(int start, double[,] distance)
        {
            //Броят на бензиностанциите
            N = distance.GetLength(0);

            // Ако бензиностанциите са 2, минималното разстояние 
            // е равно на разстоянието от първата до втората бензиностанция, 
            // умножено по 2 (заради връщането обратно към първата) 
            if (N == 2)
            {
                minTourCost = distance[1, 0] * 2;
                ranSolver = true;
            };

            this.start = start;
            this.distance = distance;
        }

        // Returns the optimal tour for the traveling salesman problem.
        public string getTour()
        {
            if (!ranSolver) solve();
            tour.Reverse();
            return String.Join("=>", tour);
        }

        // Returns the minimal tour cost.
        public double getTourCost()
        {
            if (!ranSolver) solve();
            return minTourCost;
        }

        // Solves the traveling salesman problem and caches solution.
        public void solve()
        {

            if (ranSolver) return;

            int END_STATE = (1 << N) - 1;

            // We construct a N by 2^N size matrix
            double[,] memo = new double[N, (1 << N)];

            // Add all outgoing edges from the starting node to memo table.
            for (int end = 0; end < N; end++)
            {
                if (end == start) continue;
                memo[end, (1 << start) | (1 << end)] = distance[start, end];
            }

            for (int r = 3; r <= N; r++)
            {
                foreach (int subset in combinations(r, N))
                {
                    if (notIn(start, subset)) continue;
                    for (int next = 0; next < N; next++)
                    {
                        if (next == start || notIn(next, subset)) continue;
                        int subsetWithoutNext = subset ^ (1 << next);
                        double minDist = double.PositiveInfinity;
                        for (int end = 0; end < N; end++)
                        {
                            if (end == start || end == next || notIn(end, subset)) continue;
                            double newDistance = memo[end, subsetWithoutNext] + distance[end, next];
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
            for (int i = 0; i < N; i++)
            {
                if (i == start) continue;
                double tourCost = memo[i, END_STATE] + distance[i, start];
                if (tourCost < minTourCost)
                {
                    minTourCost = tourCost;
                }
            }
            PrintMemoTbl(memo);
            int lastNode = start;
            int state = END_STATE;
            tour.Add(start);

            // Reconstruct TSP path from memo table.
            for (int i = 1; i < N; i++)
            {

                // the current number of the gass station to be added to the tour
                int node = -1;
                for (int j = 0; j < N; j++)
                {
                    if (j == start || notIn(j, state)) continue;
                    if (node == -1) node = j;
                    double prevDist = memo[node, state] + distance[node, lastNode];
                    double newDist = memo[j, state] + distance[j, lastNode];
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
            ranSolver = true;
        }

        private static bool notIn(int elem, int subset)
        {
            return ((1 << elem) & subset) == 0;
        }

        // This method generates all bit sets of size n where r bits
        // are set to one. The result is returned as a list of integer masks.
        public static List<int> combinations(int r, int n)
        {
            List<int> subsets = new List<int>();
            combinations(0, 0, r, n, subsets);
            return subsets;
        }

        // To find all the combinations of size r we need to recurse until we have
        // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
        // an element which is found after the position of our last selected element
        private static void combinations(int set, int at, int r, int n, List<int> subsets)
        {

            // Return early if there are more elements left to select than what is available.
            int elementsLeftToPick = n - at;
            if (elementsLeftToPick < r) return;

            // We selected 'r' elements so we found a valid subset!
            if (r == 0)
            {
                subsets.Add(set);
            }
            else
            {
                for (int i = at; i < n; i++)
                {
                    // Try including this element
                    set ^= (1 << i);

                    combinations(set, i + 1, r - 1, n, subsets);

                    // Backtrack and try the instance where we did not include this element
                    set ^= (1 << i);
                }
            }
        }

        public static void PrintMemoTbl(double[,] matrix)
        {
            var result = new StringBuilder();
            Console.Write(" ".PadRight(3));
            for (int h = 3; h < matrix.GetLength(1); h++)
            {
                if (h % 2 == 0) continue;
                Console.Write(Convert.ToString(h).PadLeft(2), Color.YellowGreen);
                
                Console.Write("|".PadLeft(2));
            }
            Console.WriteLine();
            for (int row = 1; row < matrix.GetLength(0); row++)
            {
                Console.Write(row.ToString().PadRight(3), Color.YellowGreen);
                for (int col = 3; col < matrix.GetLength(1); col++)
                {
                    if (col % 2 == 0) continue;
                    if(matrix[row,col]!= 0)
                    {
                        Console.Write(matrix[row, col].ToString().PadLeft(2), Color.Pink);
                    }
                    else
                    {
                        Console.Write(matrix[row, col].ToString().PadLeft(2), Color.DarkGray);
                    }
                    
                    
                    Console.ResetColor();
                    Console.Write("|".PadLeft(2));
                }
                Console.WriteLine();
            }

            
        }


    }


}
