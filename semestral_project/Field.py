import copy
import numpy as np
import matplotlib.pyplot as plt

c = 299792458  # m/s

c_us = 299.792458  # m/us


class Anchor:
    def __init__(self, name: int, number: int, coordinates: np.array, coordinates_unit: tuple, is_main: bool):
        self.is_main = is_main
        self.name = name
        self.number = number
        self.coordinates_m = coordinates
        self.coordinates_unit = coordinates_unit

    def __str__(self):
        return "{} {}".format("A" if self.is_main else "a", self.name)


class AnchorsField:
    # Assumed that there are only 9 anchors
    def __init__(self, main_anc: Anchor, other_ancs: list[Anchor], coors: np.array):
        self.anchor_names = {0x44: 0, 0x45: 1, 0x46: 2, 0x47: 3, 0x48: 4, 0x49: 5, 0x4a: 6, 0x4b: 7, 0x4c: 8}
        self.main_anc = copy.deepcopy(main_anc)
        self.ancs = copy.deepcopy(other_ancs)
        self.coors = coors  # for plotting only!

        # Time is corrected with respect to the main anchor
        self.TOF_table = None
        self.time_corrections = [None for _ in range(len(other_ancs))]
        self.time_corrections[main_anc.number] = 0

        self.init_TOF2main()
        self.TOF_table: np.array

    def init_TOF2main(self):
        self.main_anc: Anchor
        self.ancs: list[Anchor]
        res = np.zeros((len(self.ancs), len(self.ancs)))

        for i in range(len(self.ancs)):
            for j in range(len(self.ancs)):
                res[i, j] = np.linalg.norm(self.ancs[i].coordinates_m - self.ancs[j].coordinates_m) / c_us
        self.TOF_table = res

    def __str__(self):
        res_arr = [[None for _ in range(3)] for __ in range(3)]
        cds = self.main_anc.coordinates_unit
        res_arr[cds[0]][cds[1]] = self.main_anc.__str__()
        for anc in self.ancs:
            anc: Anchor
            cds = anc.coordinates_unit
            res_arr[cds[0]][cds[1]] = anc.__str__()
        res = ""
        for r in res_arr:
            r: list
            for e in r:
                e: Anchor
                res += "{}\t\t".format(e.__str__())
            res += "\n"
        return res

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(self.coors[0, :], self.coors[1, :], self.coors[2, :])
        ax.set_xlabel('x, [m]')
        ax.set_ylabel('y, [m]')
        ax.set_zlabel('z, [m]')
        plt.show()

    def update_corrections(self, n_from, n_to, t1, t2):
        n_from = self.anchor_names[n_from]
        n_to = self.anchor_names[n_to]

        # To correct any time, it is needed to ADD the corespondent t_error to received time
        if n_from == self.main_anc.number:
            t_ideal = t1 + self.TOF_table[n_from, n_to]
            correction = t_ideal - t2
            self.time_corrections[n_to] = correction
        elif n_to == self.main_anc.number:
            t_ideal = t2 - self.TOF_table[n_from, n_to]
            correction = t_ideal - t1
            self.time_corrections[n_from] = correction
        elif self.time_corrections[n_from] is not None:
            t_ideal = t1 + self.time_corrections[n_from] + self.TOF_table[n_from, n_to]
            correction = t_ideal - t2
            self.time_corrections[n_to] = correction
        elif self.time_corrections[n_to] is not None:
            t_ideal = t2 + self.time_corrections[n_to] - self.TOF_table[n_from, n_to]
            correction = t_ideal - t1
            self.time_corrections[n_from] = correction


    def get_tof(self, n_from, n_to):
        if n_from == self.main_anc.number:
            pass
        elif n_to == self.main_anc.number:
            pass
        else:
            raise Exception("Sorry, not yet implemented ger_tof for {} and {}".format(n_from, n_to))
