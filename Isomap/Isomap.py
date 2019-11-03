import numpy as np
from MDS import MDS


class Isomap(MDS):
    def __init__(self, data):
        MDS.__init__(self, data)

    def get_geodist(self, k_near):                         # Floyd Algorithm for distance computatioon
        dist = np.sqrt(self.dist)
        m = dist.shape[0]

        # Keep k nearest neighbor of a vertex, and set the remaining vertices to infinity
        for i in range(m):
            sorted_distance = np.argsort(dist[i, :])
            for j in range(k_near + 1, m):
                dist[i][sorted_distance[j]] = float('inf')

        for k in range(m):
            for i in range(m):
                for j in range(m):
                    if dist[i][k] + dist[k][j] < dist[i][j]:
                        dist[i][j] = dist[i][k] + dist[k][j]

        self.dist = np.square(dist)

    def iso_transform(self, knear, dim=2):
        if knear > self.dist.shape[0] - 1:
            print("Wrong k parameter")
            return None
        self.get_geodist(knear)
        if float('inf') in self.dist:
            print("Graph is not connected, select a larger k parameter")
            return None
        return self.mds_transform(dims=dim)


if __name__ == "__main__":
    data1 = np.arange(0, 5, 0.1)
    data2 = np.sin(data1)
    data3 = np.sin(data1)
    data = np.vstack([data1, data2, data3])
    iso = Isomap(data)
    print(iso.iso_transform(3, 1))



