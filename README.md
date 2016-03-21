# SensorPlacements

Implements methods from Krause, Singh, and Guestrin paper *Near-Optimal Sensor Placements in Gaussian Processes:
Theory, Efficient Algorithms and Empirical Studies* (2008) published in Journal of Machine Learning Research.
The paper deals with placing sensors optimally from a finite set of points. This problem is NP-hard, but they provide some heuristic algorithms with some guarantees. The optimality criteria they use is mutual information which measures the information gained at all points we are interested in.

The greedy algorithm simply sequentially picks the point that provides the largest increase in mutual information until the desired number of sensors has been placed. Code for this is included in the file `SensorPlacement.R`.
