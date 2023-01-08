#!/bin/env python3

# The principle of this measurement consists in periodically sending messages (via a wireless channel) by a moving
# object ("mark") and a set of fixed measuring points ("anchor"). Messages are timestamped when both
# sent and received. The problem is that the individual elements of the network do not have synchronous time - each
# has its own time base, therefore even the anchors exchange synchronization messages with each other in order to
# recalculate all recorded time data into a common domain.
#
# In the attachment, I am sending a map of the location of the anchors and 2 different datasets. Finally,
# compare them with each other with regard to the resulting "quality" (I'll leave the choice of indicator up to you)
# of determining the position.
#
# The syncs.csv file contains sync messages between anchors in the format
# addr_tx,addr_rx,ts_tx,ts_rx,timestamp
# where addr represents the addresses of individual anchors (_tx - transmitter, _rx - receiver) and ts are the
# respective time stamps in time units of the DWM1000 module, which corresponds to a time of 1/(128*499.2*10^6) s, i.e.
# approx. 15.65 ps (verify this information in the datasheet for DWM1000). Timestamp is the timestamp assigned when the
# message is received by the parent system (Unix time in ms).
#
# The blinks.csv file contains location messages between markers and anchors in the format
# addr_tx,addr_rx,ts_rx,ts_tx,id
# where addr_tx represents the address of the tag, addr_rx the address of the anchor that captured the given message,
# ts are the respective timestamps, and id represents the sequence number of the sent message.
#
# The goal is to design and implement a method that would be able to reconstruct the trajectory of a moving object
# from the attached datasets. It can be used in practice, for example, when detecting the deviation of a robotic
# transporter from a predefined path.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
from scipy import stats
from Field2 import *
from scipy.optimize import curve_fit


mpl.use('Qt5Agg')

A_names = {0x44: 0, 0x45: 1, 0x46: 2, 0x47: 3, 0x48: 4, 0x49: 5, 0x4a: 6, 0x4b: 7, 0x4c: 8}


def init_field(main_N: int):
    # Assume only 9 specific anchors
    anchor_names = [0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4a, 0x4b, 0x4c]

    anchors_coors = np.array([[-1.97, -12.75, -12.77, -1.81, -6.86, -1.92, -6.87, -12.27, -6.77],
                              [-8.05, -8.05, 2.75, 2.75, -2.67, -2.67, -8.05, -2.67, 2.75],
                              [2.6, 2.6, 3.13, 3.13, 2.86, 2.86, 2.6, 2.86, 3.13]])
    ancs = []

    for i in range(anchors_coors.shape[1]):
        ancs.append(Anchor(anchor_names[i],
                           i,
                           anchors_coors[:, i],
                           False))

    main_anc = Anchor(anchor_names[main_N],
                      main_N,
                      anchors_coors[:, main_N],
                      True)
    res = AnchorsField(main_anc, ancs, anchors_coors,
                       {0x44: 0, 0x45: 1, 0x46: 2, 0x47: 3, 0x48: 4, 0x49: 5, 0x4a: 6, 0x4b: 7, 0x4c: 8})
    res.prev_update[main_anc.number] = np.inf
    return res


if __name__ == "__main__":

    fname_syncs = "syncs_walls5.csv"
    fname_blinks = "blinks_walls5.csv"
    # fname_syncs = "syncs_cold_start.csv"
    # fname_blinks = "blinks_cold_start.csv"

    # Data preprocessing: read files
    syncs = np.genfromtxt(fname_syncs, delimiter=',', dtype=np.float64)
    blinks = np.genfromtxt(fname_blinks, delimiter=',', dtype=np.float64)

    # remove repeating lines, sort with respect to the UNIX time
    new_array = [tuple(row) for row in blinks]
    blinks = np.unique(new_array, axis=0)
    blinks = blinks[blinks[:, 3].argsort()]

    new_array = [tuple(row) for row in syncs]
    syncs = np.unique(new_array, axis=0)
    syncs = syncs[syncs[:, 4].argsort()]

    # make a starting time unit in sync file as time 0
    # transfer all time to us.
    syncs[:, 4] *= 1000
    blinks[:, 3] *= 1000

    initial_time = syncs[0, 4]
    syncs[:, 4] -= initial_time
    blinks[:, 3] -= initial_time

    # t6874 = syncs[syncs[:, 0] == 68]
    # t6874 = t6874[t6874[:, 1] == 74]
    # x = t6874[800:807, 3]
    # y = t6874[800:807, 2] - t6874[800:807, 3]
    #
    # # plt.plot(x, y, 'r*')
    # # plt.plot(x, predict(x), 'g.')
    # # plt.show()
    # predict = np.poly1d(np.polyfit(x[:3], y[:3], 1))
    # plt.plot(y - predict(x), 'b.', label='1d polynom')
    # predict = np.poly1d(np.polyfit(x[:3], y[:3], 2))
    # plt.plot(y - predict(x), 'g.', label='2d polynom')
    # plt.ylabel("The difference between real and predicted values, [dw]")
    # plt.xlabel("68th anchor synchronisation time, [dw]")
    # plt.title("The difference between real and predicted values, trained on 3 elements")
    # plt.legend()
    # plt.show()

    # plt.ylabel("difference between anchors local time, own units")
    # plt.xlabel("absolute time, [s]")
    # plt.legend()
    # plt.show()

    # if __name__ == "a":

    # test_anc_n = 71
    # test_anc = syncs[syncs[:, 0] == test_anc_n]
    # syncs = syncs[syncs[:, 0] != test_anc_n]
    # syncs = syncs[syncs[:, 1] != test_anc_n]

    # # # # # # # # # # # # # # #
    # # # # # main part # # # # #
    # # # # # # # # # # # # # # #
    F = init_field(4)
    i = 0
    i_max = syncs.shape[0]
    i_blinks = 0
    x_prev = np.array([np.mean(F.coors[0]), np.mean(F.coors[1]), np.mean(F.coors[2])])
    res = []
    while True:
        n_active_ancs = [el is not None for el in F.l_time_corrections]
        # n_active_ancs[test_anc_n - 68] = True
        blinks_recorded = False
        if np.all(n_active_ancs):
            nrx_ttxs_trxs_id_ntx = []
            # while test_anc[i_blinks, 4] < syncs[i, 4]:
            #     nrx_ttxs_trxs_id_ntx.append([test_anc[i_blinks, 1],
            #                                  float(test_anc[i_blinks, 2]),
            #                                  float(test_anc[i_blinks, 3]),
            #                                  int(1),
            #                                  test_anc_n])
            while blinks[i_blinks, 3] < syncs[i, 4]:
                nrx_ttxs_trxs_id_ntx.append([blinks[i_blinks, 1],
                                             float(blinks[i_blinks, 3]),
                                             float(blinks[i_blinks, 2]),
                                             int(blinks[i_blinks, 4]),
                                             blinks[i_blinks, 0]])
                blinks_recorded = True
                i_blinks += 1
                if i_blinks >= blinks.shape[0]:
                    break
                # if i_blinks >= test_anc.shape[0]:
                #     break
            if blinks_recorded:
                blinks_recorded = False
                if len(nrx_ttxs_trxs_id_ntx) < 4:
                    print("Unable to locate the tag because at least 4 anchors are needed for that in 3D")
                else:
                    nrx_ttxs_trxs_id_ntx = np.array(nrx_ttxs_trxs_id_ntx)
                    a = np.all((nrx_ttxs_trxs_id_ntx[:, 3]) == nrx_ttxs_trxs_id_ntx[0, 3])
                    b = np.all((nrx_ttxs_trxs_id_ntx[:, 1]) == nrx_ttxs_trxs_id_ntx[0, 1])
                    if not a:
                        print("not all messages has the same id")
                        continue
                    if not b:
                        unique_tstpts = list(set(nrx_ttxs_trxs_id_ntx[:, 1]))
                        ns = []
                        for element in unique_tstpts:
                            n = np.count_nonzero(nrx_ttxs_trxs_id_ntx[:, 1] == element)
                            ns.append(n)
                        ns = np.array(ns)
                        ns = np.argmax(ns)
                        major_time = unique_tstpts[ns]
                        nrx_ttxs_trxs_id_ntx = nrx_ttxs_trxs_id_ntx[nrx_ttxs_trxs_id_ntx[:, 1] == major_time, :]
                    nrx_ttxs_trxs_id_ntx = (nrx_ttxs_trxs_id_ntx[nrx_ttxs_trxs_id_ntx[:, 0].argsort()])

                    l_x = F.locate_tag_lq(nrx_ttxs_trxs_id_ntx[:, 0],
                                          nrx_ttxs_trxs_id_ntx[:, 2],
                                          x_prev, c_dw,
                                          nrx_ttxs_trxs_id_ntx[:, 1])
                    if l_x is not None:
                        # print(i)
                        x_prev = l_x
                        res.append([l_x[0], l_x[1], l_x[2], nrx_ttxs_trxs_id_ntx[0][4]])
        # # # # # # # #
        # general part
        F.update_t_corrections(syncs[i, 0], syncs[i, 1], syncs[i, 2], syncs[i, 3], syncs[i, 4])
        i += 1
        # if i >= i_max or (i_blinks >= test_anc.shape[0]):
        if i >= i_max or (i_blinks >= blinks.shape[0]):
            break

    ax = F.plot_2d()
    res = np.array(res).T
    if res.size == 0:
        print("nothing happened: anchors were not synchronized")
    else:
        # res1 = res[0:3, res[-1] == test_anc_n]
        res1 = res[:, res[-1] == 30]
        print(res1.shape)
        # ax.scatter(res1[0, :], res1[1, :], res1[2, :], 'yo')

        # compute the absolute mean error
        s = 0
        # for each in res1.T:
        #     a = F.coors[:, test_anc_n - 68]
        #     b = each
        #     assert a.shape == b.shape, "wrong dims"
        #     s += np.linalg.norm(each - F.coors[:, test_anc_n - 68])
        # s = s / res1.shape[1]

        ax.scatter(res1[0, :], res1[1, :], color='y', alpha=0.3)
        plt.gca().set_aspect('equal')
        # plt.title("{} anchor localisation with {} anchor as main. MAE={:.2f}[m]".format(test_anc_n,
        # F.main_anc.number + 68, s))
        plt.title("tag number {} localisation with {} anchor as main.".format(30, F.main_anc.number + 68))
        plt.xlabel("x axis, [m]")
        plt.ylabel("y axis, [m]")
        plt.show()
