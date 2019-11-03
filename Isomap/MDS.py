import numpy as np


class MDS():
    def __init__(self, data):
        self.data = data
        m = self.data.shape[1]
        dist = np.zeros((m, m))
        self.B = np.zeros((m, m))
        for i in range(m):
            for j in range(i, m):
                vector1 = self.data[:, i]
                vector2 = self.data[:, j]
                distance = np.sum(np.square(vector1 - vector2))
                dist[i][j] = distance

        for i in range(m):
            for j in range(i):
                dist[i][j] = dist[j][i]

        self.dist = dist

    def mds_transform(self, dims=2):
        dist = self.dist
        m = dist.shape[0]
        if m < dims:
            print("Too few samples")
            return None

        dist_i = []
        dist_j = []

        dist_ij = np.sum(dist) / (m * m)

        for i in range(m):
            dist_i.append(np.sum(dist[i, :]) / m)

        for j in range(m):
            dist_j.append(np.sum(dist[:, j]) / m)

        for i in range(m):
            for j in range(m):
                self.B[i][j] = -0.5 * (dist[i][j] - dist_i[i] - dist_j[j] + dist_ij)
        vals, vectors = np.linalg.eig(self.B)
        sorted_vals = np.argsort(vals)
        val_matrix = np.diag(vals[sorted_vals[-1: -dims-1: -1]])
        vec_matrix = vectors[:, sorted_vals[-1: -dims-1: -1]]
        Z = np.dot(np.sqrt(val_matrix), vec_matrix.T)
        return Z.T


if __name__ == "__main__":
    data = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [2, 2, 0, 0]])
    mds = MDS(data.T)
    print(mds.mds_transform(2))



