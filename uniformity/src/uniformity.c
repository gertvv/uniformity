#include <R.h>
#include <math.h>
#include <stdlib.h>

typedef struct Matrix {
	double * const data;
	int const nRow;
	int const nCol;
} Matrix;

/**
 * @param i Row index.
 * @param j Column index.
 */
inline double *get(Matrix *m, int i, int j) {
	return m->data + j * (m->nRow) + i;
}

inline void swap(int *arr, int i, int j) {
	int tmp = arr[i];
	arr[i] = arr[j];
	arr[j] = tmp;
}

/**
 * Let the rows of the matrix represent points in Euclidians space.
 * Calculate the (squared) distance between the p0-th and p1-st row.
 * Avoids calculating the sqrt() since the squared distance preserves ordering.
 */
double uniformity_sqdist(Matrix *m, int p0, int p1) { 
	double dist = 0.0;
	for (int j = 0; j < m->nCol; ++j) {
		double const d = *get(m, p0, j) - *get(m, p1, j);
		dist += d * d;
	}
	return dist;
} 

/**
 * Calculate the minimum spanning tree of the given set of points.
 * Points are given as coordinate vectors in Euclidian space. The MST is
 * calculated using Prim's algorithm and returned as a set of edges.
 * @param points A matrix in which each row represents a coordinate vector.
 * @param nDim The dimension of the Euclidian space.
 * @param nPoint The number of points.
 * @param edges The return value: predecessor list (starting with 2nd point).
 */
void uniformity_mst(double *points, int *nPoint, int *nDim, int *edges) {
	Matrix m = { points, *nPoint, *nDim };

	if (m.nRow == 0) return;

	int *q = malloc(sizeof(int) * m.nRow); // Index list of visited / non-visited vertices
	int *n = malloc(sizeof(int) * m.nRow); // Nearest neighbour list
	double *d = malloc(sizeof(double) * m.nRow); // Shortest-distance list

	for (int i = 0; i < m.nRow; ++i) {
		q[i] = i;
		n[i] = 0;
		d[i] = uniformity_sqdist(&m, 0, i);
	}
	
	for (int i = 1; i < m.nRow; ++i) {
		// Find the currently closest point
		int closest = i;
		double minDist = d[q[i]];
		for (int p = i + 1; p < m.nRow; ++p) {
			if (d[q[p]] < minDist) {
				closest = p;
				minDist = d[q[p]];
			}
		}
		edges[q[closest] - 1] = n[q[closest]] + 1; // Add edge (R indices)
		swap(q, i, closest); // Move closest point to list of visited

		// Calculate distances from newly visited point
		for (int p = i + 1; p < m.nRow; ++p) {
			double dist = uniformity_sqdist(&m, q[i], q[p]);
			if (dist < d[q[p]]) {
				d[q[p]] = dist;
				n[q[p]] = q[i];
			}
		}
	}

	free(q);
	free(n);
	free(d);
}

/**
 * For each point in the given set of points, find its nearest neighbour and
 * record the distance to the nearest neighbour.
 * @param points A matrix in which each row represents a coordinate vector.
 * @param nDim The dimension of the Euclidian space.
 * @param nPoint The number of points.
 * @param nearest For each vertex, the nearest vertex.
 * @param dist For each vertex, the (squared) distance to the nearest vertex.
 */
void uniformity_nn(double *points, int *nPoint, int *nDim, int *nearest, double *dist) {
	Matrix m = { points, *nPoint, *nDim };

	for (int i = 0; i < m.nRow; ++i) {
		nearest[i] = NA_INTEGER;
		dist[i] = NA_REAL;
	}

	for (int i = 0; i < m.nRow; ++i) {
		for (int j = i + 1; j < m.nRow; ++j) {
			double d = uniformity_sqdist(&m, i, j);
			if (ISNA(dist[i]) || d < dist[i]) {
				dist[i] = d;
				nearest[i] = j;
			}
			if (ISNA(dist[j]) || d < dist[j]) {
				dist[j] = d;
				nearest[j] = i;
			}
		}
	}
}
