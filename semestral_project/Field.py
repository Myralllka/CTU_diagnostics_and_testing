import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

c = 299792458  # m/s
c_us = 299.792458  # m/us
c_ns = 0.299792458  # m/ns


def u2ms(inp):
    """
    Transform the DWM1000 time units to microseconds.
    To have all time in microseconds (as for )
    :param inp: DWM time units
    :return: microseconds
    """
    return inp / (499.2 * 128.0)


def u2us(inp):
    """
    Transform the DWM1000 time units to microseconds.
    To have all time in microseconds (as for )
    :param inp: DWM time units
    :return: microseconds
    """
    return inp / (499.2 * 128.0 * 1e3)


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

        # Time is corrected with respect to the main anchor. To find any time in a main time domain - add the el from
        # this file
        self.TOF_table = None
        self.l_time_corr = np.array([None for _ in range(len(other_ancs))])
        #
        self.time_sync = -np.inf
        self.l_time_corr[main_anc.number] = 0

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
        # plt.show()
        return ax

    def update_corrections(self, n_from, n_to, t1, t2):
        """

        :param n_from: node number of a tx
        :param n_to: node number of a rx
        :param t1: RTC time in ms
        :param t2: RTC time in ms
        :return:
        """
        n_from = self.anchor_names[n_from]
        n_to = self.anchor_names[n_to]

        # To correct any time, it is needed to ADD the corespondent t_error to received time
        if n_from == self.main_anc.number:
            t_ideal = t1 + self.TOF_table[n_from, n_to]
            correction = t_ideal - t2
            self.l_time_corr[n_to] = correction
        elif n_to == self.main_anc.number:
            t_ideal = t2 - self.TOF_table[n_from, n_to]
            correction = t_ideal - t1
            self.l_time_corr[n_from] = correction
        elif self.l_time_corr[n_from] is not None:
            t_ideal = t1 + self.l_time_corr[n_from] + self.TOF_table[n_from, n_to]
            correction = t_ideal - t2
            self.l_time_corr[n_to] = correction
        elif self.l_time_corr[n_to] is not None:
            t_ideal = t2 + self.l_time_corr[n_to] - self.TOF_table[n_from, n_to]
            correction = t_ideal - t1
            self.l_time_corr[n_from] = correction

    def locate_tag(self, tag_n: int, idxs: np.array, times: np.array):
        """

        :param tag_n: number of a tag to locate
        :param idxs: sorted idxs of anchors
        :param times: time when each anchor received the message from the moving tag
        :return:
        """
        # Look for x, y, z coordinates of a moving tag.
        idxs = [self.anchor_names[arg] for arg in idxs]
        # if not (self.main_anc.number in idxs):
        #     print("There is no main anchor in dataset. Feature is not supported yet.")
        #     return

        # extract anchors coordinates
        x_s = self.coors[0, idxs]
        y_s = self.coors[1, idxs]
        z_s = self.coors[2, idxs]

        # d = c * delta_t
        t = times.copy() + self.l_time_corr[idxs].copy()
        t = t.reshape(t.shape[0], 1)
        t_m = t @ np.ones(t.shape).T
        d_m = t_m - t_m.T
        d_m *= c_us

        bounds = [(-14, -10, 0), (0, 3, 6)]
        # initial guess - mean of inp
        x_in, y_in, z_in = np.mean(x_s), np.mean(y_s), np.mean(z_s)
        F = functions(x_s, y_s, z_s, d_m)
        J = jacobian(x_s, y_s, z_s, d_m)

        x = opt.least_squares(F, x0=np.array([x_in, y_in, z_in]), jac=J, bounds=bounds, method='dogbox')

        return x.x


def functions(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
    """ Given anchors at Xs[i], Ys[i], Zs[i] and TDoA between observers d01, d02, d12, this closure
        returns a function that evaluates the system of N hyperbolas for given x, y, z.
    """

    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            for j in range(i + 1, Xs.shape[0]):
                res.append(np.sqrt((x - Xs[j]) ** 2 + (y - Ys[j]) ** 2 + (z - Zs[j]) ** 2) -
                           np.sqrt((x - Xs[i]) ** 2 + (y - Ys[i]) ** 2 + (z - Zs[i]) ** 2) -
                           Ds[i, j])

        return res

    return fn


def jacobian(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
    """
    Collection of the partial derivatives of each function with respect to each independent variable.
    """

    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            for j in range(i + 1, Xs.shape[0]):
                nom = np.sqrt((x - Xs[j]) ** 2 + (y - Ys[j]) ** 2 + (z - Zs[j]) ** 2)
                den = np.sqrt((x - Xs[i]) ** 2 + (y - Ys[i]) ** 2 + (z - Zs[i]) ** 2)

                adx = (x - Xs[j]) / nom - (x - Xs[i]) / den
                ady = (y - Ys[j]) / nom - (y - Ys[i]) / den
                adz = (z - Zs[j]) / nom - (z - Zs[i]) / den

                res.append([adx, ady, adz])
        return res

    return fn
