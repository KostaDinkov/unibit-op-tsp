using System;
using System.Collections.Generic;
using System.Linq;

namespace GasStations
{
    /*
     * Solution to the "Gas Stations" problem,
     * presented in the Programming Basics course, semester one, Unibit university, Sofia, Bulgaria
     * The implementation is based on the algorithm by Bellman-Held-Karp for the traveling salesman problem.
     * Note that the code below is heavily commended due to the educational context of the exercise
     * 
     * The algorithm used here achieves time complexity of O(n^2 * 2^n) and space complexity: O(n * 2^n)
     * @author Kosta Dinkov
     */
    public class Program
    {
        public static void Main()
        {
            const int startNode = 0;

            // there will be T number of tests that will check our solution
            // where 1 <= T <= 10 
            int testCount = int.Parse(Console.ReadLine());

            // for each of the tests we read the test data, compute and print the result to the console
            for (int test = 0; test < testCount; test++)
            {
                int[,] distanceMatrix = GetTestData();
                var solver = new TspDp(startNode, distanceMatrix);
                Console.WriteLine(solver.TourCost);
            }
        }

        /// <summary>
        /// Reads test data from the console and returns distance matrix based on the input. 
        /// The console input must be in the format specified in the problem specification
        /// </summary>
        private static int[,] GetTestData()
        {
            int nodeCount = int.Parse(Console.ReadLine());
            int[,] distanceMatrix = new int[nodeCount + 1, nodeCount + 1];
            for (int node = 0; node < nodeCount; node++)
            {
                int[] distances = Console.ReadLine()
                    .Split(' ', StringSplitOptions.RemoveEmptyEntries)
                    .Select(int.Parse).ToArray();

                for (int nextNode = node + 1, i = 0; i < distances.Length; nextNode++, i++)
                {
                    distanceMatrix[node, nextNode] = distances[i];
                    distanceMatrix[nextNode, node] = distances[i];
                }
            }
            
            return distanceMatrix;
        }
    }
    
    public class TspDp
    {
        private readonly int[,] distance;
        private readonly int[,] memo;
        private readonly int nodeCount;
        private readonly int startNode;

        public TspDp(int startNode, int[,] distance)
        {
            nodeCount = distance.GetLength(0);
            TourCost = ushort.MaxValue;
            this.startNode = startNode;
            this.distance = distance;
            memo = new int[nodeCount, 1 << nodeCount];
            ComputeMinTourCost();
        }

        // Returns the minimal tour cost.
        public int TourCost { get; private set; }
        
        /// <summary>
        /// Computes the minimum cost cycle from a starting
        /// node that visits all other nodes and 
        /// </summary>
        private void ComputeMinTourCost()
        {
            // in the case of only 2 nodes we return the doubled distance between them
            if (nodeCount == 2)
            {
                TourCost = distance[1, 0] * 2;
                return;
            }

            // the fullSet contains all nodes,
            // therefore its bit string for 5 nodes will look like 11111 which is 31 in decimal
            var fullSet = (1 << nodeCount) - 1;

            // We start filling our memo table with the distances from the starting node to all other nodes,
            // which are known from the distance matrix
            for (var nextNode = 0; nextNode < nodeCount; nextNode++)
            {
                if (nextNode == startNode) continue;

                // We use bitwise operations to generate the column index of the memo table
                // for example if nextNode is 3 (of total 4) then
                // (1 << 0) will generate 00001, (1 << nextNode) will generate 01000,
                // and the bitwise or between them (00001 | 01000) will evaluate to 01001
                // which will encode the distance between node 0 and node 3, based on the indexes of 
                // the resulting bit string that have 1 in them ( we will read bit strings from right to left,
                // so in this case we have 1s at the 0 and 3 index of the bit string)
                // the resulting cell of the memo table will be memo[3, 0b01001] or memo[3,9]
                memo[nextNode, (1 << startNode) | (1 << nextNode)] = distance[startNode, nextNode];
            }

            // Let s be the size of the subset of nodes that we will solve for.
            // We will gradually solve for bigger s, and we will reuse the previous results for smaller s
            for (var s = 3; s <= nodeCount; s++)
            {
                // We generate all subsets of size s of the full set of size nodeCount
                foreach (var subset in Combinations(s, nodeCount))
                {
                    // foreach of the generated subsets, we check if the starting node is in the subset
                    // and if not, then we skip and continue with the next subset
                    if (NotIn(startNode, subset)) continue;

                    // we can think of nextEnd as the last node of the path of the current subset
                    // for example if we have a subset 11101 and nextEnd is 3
                    // that means that starting from 0 we visit nodes 2 and 4 in some order and end on 3
                    for (var nextEnd = 0; nextEnd < nodeCount; nextEnd++)
                    {
                        // Naturally nextEnd cannot be the starting node, 
                        // and must be in the current subset. As from the example above
                        // if subset = 11101 and nextEnd is 1, there cannot be a path in this subset that ends on 1
                        // since 1 is not in the subset (we have 0 at index 1)
                        if (nextEnd == startNode || NotIn(nextEnd, subset)) continue;

                        // This is the part that we reuse calculations from a previous iteration of s
                        // If s is 3 than we have already calculated the minimum paths for all subsets of size 2
                        // and we can reuse the results.
                        // So if the current subset is 11101, then the previous subset, that does not include nextEnd=3
                        // is subMinusNextEnd = 10101
                        var subMinusNextEnd = subset ^ (1 << nextEnd);
                        var minDist = int.MaxValue;

                        // The lastEnd is the last node in the minimum cost path from the previous subset
                        // To continue the example from above, the last subset that does not include nextNode
                        // was 10101, so the lastEnd can be either 2 or 4 for the 0-2-4 path or 0-4-2 path.
                        for (var lastEnd = 0; lastEnd < nodeCount; lastEnd++)
                        {
                            // Again, lastEnd cannot be the starting node, or the nextNode that we are computing for
                            // Also lastEnd must be in the subMinusNextEnd subset
                            // Otherwise we skip and continue to check for another possible lastEnd
                            if (lastEnd == startNode || lastEnd == nextEnd || NotIn(lastEnd, subMinusNextEnd)) continue;

                            // Below we can see that we take the already computed in a previous iteration 
                            // minimum value for the distance of a path that starts at 0, visits all nodes
                            // in subMinusNextEnd and ends at lastEnd
                            // Lets suppose that the current subset is 11111, so we are at the final
                            // possible value for s with n=5. Let the nextEnd be 3. The previous subset that
                            // we have a memoized result in the memo table and does not include node 3 is 10111
                            // So the possible values for lastEnd are 1,2 and 4
                            // For each of this possible values we have a min cost path that starts at 0
                            // visits all other nodes and ends at lastEnd. For example, we have the min cost path
                            // for 0-{1,2}-4, no mater what the order of 1 and 2 is.
                            // (We have found the exact order of 1 and 2 in the previous iteration for s)
                            // This result is at memo[4,0b10111] or memo[4, 23]
                            // Now for each possible value for lastEnd, we sum the min path from 0 to lastEnd
                            // with the distance from lastEnd to nextEnd

                            var newDistance = memo[lastEnd, subMinusNextEnd] + distance[lastEnd, nextEnd];

                            // Of all the possible results we find the one with the minimum value ...
                            if (newDistance < minDist)
                            {
                                minDist = newDistance;
                            }
                        }

                        // ... and add it to the memo table. 
                        // As per the example above, if  subset = 11111, nextEnd=3, subMinusNextEnd = 10111
                        // memo[3,31] = minimum ot the (memo[lastEnd,23] + distance[lastEnd,nextEnd]), for all possible values of lastEnd
                        memo[nextEnd, subset] = minDist;
                    }
                }
            }

            // When we complete the iterations over subsets of all sizes up to s = n
            // we end up with a total of n-1 minimum cost paths, that start at 0
            // and end at each of the nodes {1,2,...,n}
            // So to find the optimal tour, we need to connect each of those paths back to 0
            // while adding the cost of that last edge that connects to 0, to the cost of the path.
            // For each path, this will result in a tour from 0 to 0 that visits all other nodes.
            // Finally we choose the minimal cost tour as our optimal.
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

        /// <summary>
        /// Checks if a node is not in a subset of nodes
        /// </summary>
        public static bool NotIn(int node, int subset)
        {
            return ((1 << node) & subset) == 0;
        }

        /// <summary>
        /// Generates all subsets with size r of the full set of nodes.
        /// The subsets are represented by bit strings of size n where only r positions are set to 1.
        /// Returns list with decimal representations of the bit strings.
        /// </summary>
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
    }
}