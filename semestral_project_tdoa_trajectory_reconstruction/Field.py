import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

c = 299792458  # m/s
c_us = 299.792458  # m/us

# c_ns = 0.299792458  # m/ns
dw = 1 / (499.2 * 128.0 * 1e6)  # s
c_dw = c * dw  # m/dw = 0.00469176397 m/dw


def dw2us(t):
    return t / (499.2 * 128.0)


def dw2s(t):
    return t / (499.2 * 128.0 * 1e6)


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
    def __init__(self, main_anc: Anchor, other_ancs: list[Anchor], coors: np.array, names):
        self.anchor_names = names
        self.main_anc = copy.deepcopy(main_anc)
        self.ancs = copy.deepcopy(other_ancs)
        self.coors = coors  # for plotting only!

        # Time is corrected with respect to the main anchor. To find any time in a main time domain - add the el from
        # this file
        self.TOF_table = None
        self.TOF_dist = None
        self.l_time_corr = np.array([None for _ in range(len(other_ancs))])
        #
        self.prev_update = [-np.inf for _ in range(9)]
        self.time_sync = -np.inf
        self.l_time_corr[main_anc.number] = 0

        # TODO: specify speed!
        self.init_TOF2main(c_dw)
        self.TOF_table: np.array
        self.TOF_dist: np.array
        self.time_correction_functions = [None for _ in range(9)]

        self.time_corr_prev_local_t = [0 for _ in range(9)]

    def init_TOF2main(self, speed):
        self.main_anc: Anchor
        self.ancs: list[Anchor]
        res = np.zeros((len(self.ancs), len(self.ancs)))

        for i in range(len(self.ancs)):
            for j in range(len(self.ancs)):
                res[i, j] = np.linalg.norm(self.ancs[i].coordinates_m - self.ancs[j].coordinates_m)
        self.TOF_dist = res
        self.TOF_table = np.round(res / speed)
        # self.TOF_table = res / speed

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(self.coors[0, :], self.coors[1, :], self.coors[2, :])
        ax.set_xlabel('x, [m]')
        ax.set_ylabel('y, [m]')
        ax.set_zlabel('z, [m]')
        return ax

    def update_corrections(self, n_from, n_to, t_from, t_to, abs_t):
        """

        :param abs_t:
        :param n_from: node number of a tx
        :param n_to: node number of a rx
        :param t_from: RTC time in ms
        :param t_to: RTC time in ms
        :return:
        """
        n_from = self.anchor_names[n_from]
        n_to = self.anchor_names[n_to]

        # To correct any time, it is needed to ADD the corespondent t_error to received time
        if n_from == self.main_anc.number:
            t_ideal = t_from + self.TOF_table[n_from, n_to]
            correction = t_ideal - t_to
            self.upd_time_corr(n_to, t_to, correction, abs_t)
            self.prev_update[n_to] = abs_t
        elif n_to == self.main_anc.number:
            t_ideal = t_to - self.TOF_table[n_from, n_to]
            correction = t_ideal - t_from
            self.upd_time_corr(n_from, t_from, correction, abs_t)
            self.prev_update[n_from] = abs_t
        # elif ((self.l_time_corr[n_from] is not None) and
        #       (self.l_time_corr[n_to]) is None):
        #     t_ideal = t_from + self.get_time_corrections(n_from, abs_t, t_from) + self.TOF_table[n_from, n_to]
        #     correction = t_ideal - t_to
        #     self.upd_time_corr(n_to, t_to, correction, abs_t)
        # elif ((self.l_time_corr[n_to] is not None) and
        #       (self.l_time_corr[n_from] is None)):
        #     t_ideal = t_to + self.get_time_corrections(n_to, abs_t, t_to) - self.TOF_table[n_from, n_to]
        #     correction = t_ideal - t_from
        #     self.upd_time_corr(n_from, t_from, correction, abs_t)
        # elif ((self.l_time_corr[n_from] is not None) and
        #       (self.l_time_corr[n_to] is not None)):
        #     if self.prev_update[n_from] < self.prev_update[n_to]:
        #         t_ideal_from = t_to + self.get_time_corrections(n_to, abs_t, t_to) - self.TOF_table[n_from, n_to]
        #         correction_from = t_ideal_from - t_from
        #         self.upd_time_corr(n_from, t_from, correction_from, abs_t)
        #     else:
        #         t_ideal_to = t_from + self.get_time_corrections(n_from, abs_t, t_from) + self.TOF_table[n_from, n_to]
        #         correction_to = t_ideal_to - t_to
        #         self.upd_time_corr(n_to, t_to, correction_to, abs_t)

    def upd_time_corr(self, idx: int, sync_time, correction_time, abs_time):
        if self.l_time_corr[idx] is None:
            self.l_time_corr[idx] = correction_time
            return

        def getlinear(x, y):
            def inner(x1):
                return m * x1 + b

            top = (len(x) * np.sum(x * y) - np.sum(x) * np.sum(y))
            btm = (len(x) * np.sum(x * x) - np.sum(x) * np.sum(x))
            m = top / btm
            b = (np.sum(y) - m * np.sum(x)) / len(x)
            return inner

        # # x - local time
        # # y - error
        predict = getlinear(np.array([self.time_corr_prev_local_t[idx], sync_time]),
                            np.array([self.l_time_corr[idx], correction_time]))

        self.time_correction_functions[idx] = predict
        self.time_corr_prev_local_t[idx] = sync_time
        self.l_time_corr[idx] = correction_time
        # self.time_corr_prev_error[idx] = correction_time

    def get_time_corrections(self, idx: int, abs_t: int, local_t):
        if self.time_correction_functions[idx] is None:
            return self.l_time_corr[idx]
        else:
            return self.time_correction_functions[idx](local_t)

    def locate_tag(self, idxs: np.array, times: np.array, x_prev: np.array, speed, t_abs: int):
        """

        :param t_abs:
        :param x_prev: previous location of X
        :param idxs: sorted idxs of anchors
        :param times: time when each anchor received the message from the moving tag
        :param speed: speed of traveling
        :return:
        """
        # Look for x, y, z coordinates of a moving tag.
        idxs = [self.anchor_names[arg] for arg in idxs]

        # extract anchors coordinates
        x_s = self.coors[0, idxs]
        y_s = self.coors[1, idxs]
        z_s = self.coors[2, idxs]

        # d = c * delta_t
        corrections = np.array([self.get_time_corrections(i, t_abs, times[i]) for i in idxs])
        t = times.copy() + corrections
        # dsts = self.TOF_dist[4, idxs]
        t = t.reshape(t.shape[0], 1)
        t_m = t @ np.ones(t.shape).T
        d_m = np.abs(t_m - t_m.T)
        d_m = d_m * speed

        bounds = [(-18, -13, -2),
                  (6, 5, 5)]

        F = functions(x_s, y_s, z_s, d_m)
        J = jacobian(x_s, y_s, z_s, d_m)

        x = opt.least_squares(F,
                              x0=x_prev,
                              jac=J,
                              # bounds=bounds,
                              # loss="soft_l1",
                              # method='dogbox',
                              # tr_solver="exact"
                              )

        res_coors = x.x
        return res_coors


def functions(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
    """
    Given anchors at Xs[i], Ys[i], Zs[i] and TDoA between observers, this
    returns a function that evaluates the system of N hyperbolas for given x, y, z.
    """

    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            for j in range(i + 1, Xs.shape[0]):
                a1 = np.linalg.norm(np.array([x, y, z]) - np.array([Xs[j], Ys[j], Zs[j]]))
                a2 = np.linalg.norm(np.array([x, y, z]) - np.array([Xs[i], Ys[i], Zs[i]]))
                res.append(np.abs(a1 - a2) - Ds[i, j])

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
                sq_i = np.sqrt(np.power(x - Xs[i], 2.) + np.power(y - Ys[i], 2.) + np.power(z - Zs[i], 2.))
                sq_j = np.sqrt(np.power(x - Xs[j], 2.) + np.power(y - Ys[j], 2.) + np.power(z - Zs[j], 2.))

                common = (sq_j - sq_i) / np.abs(sq_j - sq_i)

                adx = ((x - Xs[j]) / sq_j - (x - Xs[i]) / sq_i) * common
                ady = ((y - Ys[j]) / sq_j - (y - Ys[i]) / sq_i) * common
                adz = ((z - Zs[j]) / sq_j - (z - Zs[i]) / sq_i) * common

                res.append([adx, ady, adz])
        return res

    return fn
