import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

c = 299792458  # m/s
c_us = 299.792458  # m/us

# c_ns = 0.299792458  # m/ns
dw = 1 / (499.2 * 128.0 * 1e6)  # s
# c_dw = c * dw  # m/dw = 0.00469176397 m/dw
c_dw = c / (499.2 * 128.0 * 1e6)  # s


def dw2us(t):
    return t / (499.2 * 128.0)


def dw2s(t):
    return t / (499.2 * 128.0 * 1e6)


class quu:
    def __init__(self, capacity):
        self.capacity = capacity
        self.storage = np.zeros(capacity)

    def append(self, item):
        self.storage = np.delete(np.append([item], self.storage), self.capacity)

    def get_data(self):
        return copy.deepcopy(self.storage)

    def last(self):
        return self.storage[0]


class Anchor:
    def __init__(self, name: int, number: int, coordinates: np.array, is_main: bool):
        self.is_main = is_main
        self.name = name
        self.number = number
        self.coordinates_m = coordinates


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
        self.init_TOF2main(c_dw)
        self.TOF_table: np.array
        self.l_time_corrections = [None for _ in range(9)]
        self.l_time_corrections[main_anc.number] = 0
        self.prev_update = [-np.inf for _ in range(9)]

        # for the time synchronisation
        # to fit a curve - at least 3 pts needed
        hist_size = 3
        self.hist_time_corrections = [quu(hist_size) for _ in range(9)]
        self.hist_local_times = [quu(hist_size) for _ in range(9)]

    def init_TOF2main(self, speed):
        res = np.zeros((len(self.ancs), len(self.ancs)))
        for i in range(len(self.ancs)):
            for j in range(len(self.ancs)):
                res[i, j] = np.linalg.norm(self.ancs[i].coordinates_m - self.ancs[j].coordinates_m)
        self.TOF_table = np.round(res / speed)

    def plot_2d(self):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.scatter(self.coors[0, :], self.coors[1, :])
        b = np.array(list(self.anchor_names.values())) + 68
        for t in range(len(b)):
            ax.text(self.coors[0, t] + 0.1, self.coors[1, t] + 0.1, str(b[t]))
        ax.set_xlabel('x, [m]')
        ax.set_ylabel('y, [m]')
        return ax

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(self.coors[0, :], self.coors[1, :], self.coors[2, :])
        ax.set_xlabel('x, [m]')
        ax.set_ylabel('y, [m]')
        ax.set_zlabel('z, [m]')
        return ax

    def get_t_correction(self, idx, t_local):
        # curve predict
        if idx == self.main_anc.number:
            return 0
        x = copy.deepcopy(self.hist_local_times[idx].get_data())
        y = copy.deepcopy(self.hist_time_corrections[idx].get_data())
        if np.any(x == 0):
            print("waiting for history to be full...")
            return self.l_time_corrections[idx]
        if np.any(y == 0):
            print("waiting for history to be full...")
            return self.l_time_corrections[idx]

        predict = np.poly1d(np.polyfit(x, y, 2))

        res = predict(t_local)
        return res

    def get_t_correction_1d(self, idx, t_local):
        """
        implementing the improvements for the time synchronisation from the report
        """

        # line predict
        def getlinear(x, y):
            def inner(x1):
                return m * x1 + b

            # sum_y =
            n = len(x)
            top = (n * sum(x * y) - sum(x) * sum(y))
            btn = (n * sum(x * x) - sum(x) * sum(x))

            m = top / btn
            b = (sum(y) - m * sum(x)) / n
            return inner

        x, y = self.hist_local_times[idx].get_data(), self.hist_time_corrections[idx].get_data()
        if np.any(x == 0):
            print("fuck")
            return self.l_time_corrections[idx]
        if np.any(y == 0):
            print("shit")
            return self.l_time_corrections[idx]

        predict = getlinear(x, y)
        a = predict(self.hist_local_times[idx].last())
        b = self.hist_time_corrections[idx].last()
        res = predict(t_local)
        return res

    def update_t_corrections_regressive(self, idx, t_abs, t_local, t_correction):
        if self.hist_local_times[idx].last() != t_local:
            self.hist_local_times[idx].append(t_local)
            self.l_time_corrections[idx] = t_correction
            # if self.hist_time_corrections[idx].last() != t_correction:
            self.hist_time_corrections[idx].append(t_correction)

    def update_t_corrections(self, n_from, n_to, t_from, t_to, t_abs=0):
        """
        implementing all formulas for time synchronisation from the report.
        :param t_abs:
        :param n_from: node number of a tx
        :param n_to: node number of a rx
        :param t_from: RTC time in ms
        :param t_to: RTC time in ms
        :return:
        """

        n_from = self.anchor_names[n_from]
        n_to = self.anchor_names[n_to]

        if n_from == self.main_anc.number:
            t_ideal = t_from + self.TOF_table[n_from, n_to]
            correction = t_ideal - t_to
            self.update_t_corrections_regressive(n_to, t_abs, t_to, correction)
            self.prev_update[n_to] = t_abs
            deb = t_to + correction - self.TOF_table[n_from, n_to]
            assert abs(deb - t_from) < 1, "wrong sync main n_from"
        elif n_to == self.main_anc.number:
            t_ideal = t_to - self.TOF_table[n_from, n_to]
            correction = t_ideal - t_from
            self.update_t_corrections_regressive(n_from, t_abs, t_from, correction)
            deb = t_from + correction + self.TOF_table[n_from, n_to]
            assert abs(deb - t_to) < 1, "wrong sync main n_to"
            self.prev_update[n_from] = t_abs
        elif ((self.l_time_corrections[n_from] is not None) and
              (self.l_time_corrections[n_to]) is None):
            t_ideal = t_from + self.get_t_correction(n_from, t_from) + self.TOF_table[n_from, n_to]
            correction = t_ideal - t_to
            self.update_t_corrections_regressive(n_to, t_abs, t_to, correction)
            self.prev_update[n_to] = self.prev_update[n_from]
            deb = t_to + correction - self.TOF_table[n_from, n_to]
            assert abs((t_from + self.get_t_correction(n_from, t_from)) - deb) < 1, "wrong sync from2None"
        elif ((self.l_time_corrections[n_to] is not None) and
              (self.l_time_corrections[n_from] is None)):
            t_ideal = t_to + self.get_t_correction(n_to, t_to) - self.TOF_table[n_from, n_to]
            correction = t_ideal - t_from
            self.update_t_corrections_regressive(n_from, t_abs, t_from, correction)
            self.prev_update[n_from] = self.prev_update[n_to]
            deb = t_from + correction + self.TOF_table[n_from, n_to]
            sec = (t_to + self.get_t_correction(n_to, t_to))
            assert abs(sec - deb) < 1, "wrong sync None2to"
        elif ((self.l_time_corrections[n_from] is not None) and
              (self.l_time_corrections[n_to] is not None)):
            if self.prev_update[n_from] < self.prev_update[n_to]:
                t_ideal_from = t_to + self.get_t_correction(n_to, t_to) - self.TOF_table[n_from, n_to]
                correction_from = t_ideal_from - t_from
                self.update_t_corrections_regressive(n_from, t_abs, t_from, correction_from)
            else:
                t_ideal_to = t_from + self.get_t_correction(n_from, t_from) + self.TOF_table[n_from, n_to]
                correction_to = t_ideal_to - t_to
                self.update_t_corrections_regressive(n_to, t_abs, t_to, correction_to)

    def locate_tag_tdoa_noneabs(self, idxs: np.array, t_local_noncorrected: np.array, x_prev: np.array, speed, t_abs=0):
        """

        :param idxs: sorted idxs of anchors
        :param t_local_noncorrected:
        :param x_prev: previous location of X
        :param speed: speed of traveling
        :param t_abs:
        :return:
        """
        # Look for x, y, z coordinates of a moving tag.
        idxs = np.array([self.anchor_names[arg] for arg in idxs])

        # extract anchors coordinates
        x_s = self.coors[0, idxs]
        y_s = self.coors[1, idxs]
        z_s = self.coors[2, idxs]

        corrections = [self.get_t_correction(idxs[i], t_local_noncorrected[i]) for i in range(len(idxs))]
        t = copy.deepcopy(t_local_noncorrected) + copy.deepcopy(corrections)
        t = t.reshape(t.shape[0], 1)
        t_m = t @ np.ones(t.shape).T
        d_m = - t_m + t_m.T
        d_m = d_m * speed

        F = tdoa_noneabs(x_s, y_s, z_s, d_m)
        J = tdoa_jacob_noneabs(x_s, y_s, z_s, d_m)

        bounds = [
            [-16, -10, 0],
            [4, 4, 5]
            ]

        x = opt.least_squares(F,
                              x0=copy.deepcopy(x_prev),
                              jac=J,
                              # bounds=bounds,
                              loss="soft_l1",
                              method='dogbox',
                              tr_solver="exact"
                              )

        res_coors = x.x
        bounds = np.array(bounds) * 2
        if np.any(res_coors < bounds[0]) or np.any(res_coors > bounds[1]):
            # print("wrong.")
            return None
        return res_coors

    def locate_tag_tdoa(self, idxs: np.array, t_local_noncorrected: np.array, x_prev: np.array, speed, t_abs=0):
        """

        :param idxs: sorted idxs of anchors
        :param t_local_noncorrected:
        :param x_prev: previous location of X
        :param speed: speed of traveling
        :param t_abs:
        :return:
        """
        # Look for x, y, z coordinates of a moving tag.
        idxs = np.array([self.anchor_names[arg] for arg in idxs])

        # extract anchors coordinates
        x_s = self.coors[0, idxs]
        y_s = self.coors[1, idxs]
        z_s = self.coors[2, idxs]

        corrections = [self.get_t_correction(idxs[i], t_local_noncorrected[i]) for i in range(len(idxs))]
        t = copy.deepcopy(t_local_noncorrected) + copy.deepcopy(corrections)
        t = t.reshape(t.shape[0], 1)
        t_m = t @ np.ones(t.shape).T
        d_m = np.abs(t_m - t_m.T)
        d_m = d_m * speed

        F = tdoa(x_s, y_s, z_s, d_m)
        J = tdoa_jacob(x_s, y_s, z_s, d_m)

        bounds = [
            [-16, -10, 0],
            [4, 4, 5]
            ]

        x = opt.least_squares(F,
                              x0=copy.deepcopy(x_prev),
                              jac=J,
                              # bounds=bounds,
                              loss="soft_l1",
                              method='dogbox',
                              tr_solver="exact"
                              )

        res_coors = x.x
        bounds = np.array(bounds)
        if np.any(res_coors < bounds[0]) or np.any(res_coors > bounds[1]):
            return None
        return res_coors

    def locate_tag_lq(self, idxs: np.array, t_local_noncorrected: np.array, x_prev: np.array, speed, t_abs=0):
        """

        :param idxs: sorted idxs of anchors
        :param t_local_noncorrected:
        :param x_prev: previous location of X
        :param speed: speed of traveling
        :param t_abs:
        :return:
        """
        idxs = np.array([self.anchor_names[arg] for arg in idxs])

        # extract anchors coordinates
        x_s = self.coors[0, idxs]
        y_s = self.coors[1, idxs]
        z_s = self.coors[2, idxs]

        corrections = [self.get_t_correction(idxs[i], t_local_noncorrected[i]) for i in range(len(idxs))]
        t = copy.deepcopy(t_local_noncorrected) + copy.deepcopy(corrections)

        t = t.reshape(t.shape[0], 1)
        t_m = t @ np.ones(t.shape).T
        d_m = t_m - t_m.T
        d_m = d_m * speed

        res = least_squares_fit(x_s, y_s, z_s, d_m)

        r = res[np.argmin(res[:, 3])]
        r = r[:-1]
        bounds = [
            [-16, -10, -4],
            [4, 4, 5]
            ]
        bounds = np.array(bounds) * 100
        if np.any(r < bounds[0]) or np.any(r > bounds[1]):
            # print("wrong.")
            return None
        if r is not None:
            return r
        print("error!!!!")
        return None


def least_squares_fit(Xs, Ys, Zs, Rs):
    """

    :param Xs:
    :param Ys:
    :param Zs:
    :param Rs: Rs[i, j] time of flight from node i to j
    :return:
    """

    Xs = copy.deepcopy(Xs)
    Ys = copy.deepcopy(Ys)
    Zs = copy.deepcopy(Zs)
    Rs = copy.deepcopy(Rs)

    def xp_n(i, j):
        return Xs[i] - Xs[j]

    def yp_n(i, j):
        return Ys[i] - Ys[j]

    def zp_n(i, j):
        return Zs[i] - Zs[j]

    def kp(i, j):
        return Xs[i] ** 2 - Xs[j] ** 2 + Ys[i] ** 2 - Ys[j] ** 2 + Zs[i] ** 2 - Zs[j] ** 2

    res = []
    for n in range(Xs.shape[0]):
        A = np.array([[xp_n(n, _), yp_n(n, _), zp_n(n, _), Rs[n, _]] for _ in range(Xs.shape[0]) if _ != n])
        k = np.array([[kp(n, _) + Rs[n, _] ** 2] for _ in range(Xs.shape[0]) if _ != n]) / 2
        try:
            inv = np.linalg.inv(A.T @ A)
        except:
            continue
        p = (inv @ A.T @ k)
        res.append(p.ravel())

    res = np.array(res)

    return res


def fit_point_between_circles(Xs, Ys, Zs, Rs):
    """
    Given spheres at Xs[i], Ys[i], Zs[i] with radius R[i], this
    """

    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            dst = np.linalg.norm(np.array([x, y, z]) - np.array([Xs[i], Ys[i], Zs[i]]))
            res.append(abs(dst - Rs[i]))
        return res

    return fn


def fit_point_between_circles_jacob(Xs, Ys, Zs, Rs):
    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            n = np.linalg.norm(np.array([x, y, z]) - np.array(Xs[i], Ys[i], Zs[i]))
            common = (n - Rs[i]) / (n * np.abs(n))
            dx = (x - Xs[i]) * common
            dy = (y - Ys[i]) * common
            dz = (z - Ys[i]) * common
            res.append([dx, dy, dz])
        return res

    return fn


def tdoa(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
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


def tdoa_jacob(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
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


def tdoa_noneabs(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
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
                res.append(a1 - a2 - Ds[i, j])

        return res

    return fn


def tdoa_jacob_noneabs(Xs: np.array, Ys: np.array, Zs: np.array, Ds: np.array):
    """
    Collection of the partial derivatives of each function with respect to each independent variable.
    """

    def fn(args):
        x, y, z = args
        res = []
        for i in range(Xs.shape[0]):
            for j in range(i + 1, Xs.shape[0]):
                den1 = np.linalg.norm(np.array([x, y, z]) - np.array([Xs[j], Ys[j], Zs[j]]))
                den2 = np.linalg.norm(np.array([x, y, z]) - np.array([Xs[i], Ys[i], Zs[i]]))

                adx = (x - Xs[j]) / den1 - (x - Xs[i]) / den2
                ady = (y - Ys[j]) / den1 - (y - Ys[i]) / den2
                adz = (z - Zs[j]) / den1 - (z - Zs[i]) / den2

                res.append([adx, ady, adz])
        return res

    return fn
