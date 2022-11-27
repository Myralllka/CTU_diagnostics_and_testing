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
from Field import AnchorsField, Anchor, u2ms

mpl.use('Qt5Agg')


def init_field(main_N: int):
    # Assume only 9 specific anchors
    anchor_names = [0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4a, 0x4b, 0x4c]

    anchors_coors = np.array([[-1.97, -12.75, -12.77, -1.81, -6.86, -1.92, -6.87, -12.27, -6.77],
                              [-8.05, -8.05, 2.75, 2.75, -2.67, -2.67, -8.05, -2.67, 2.75],
                              [2.6, 2.6, 3.13, 3.13, 2.86, 2.86, 2.6, 2.86, 3.13]])
    anchors_unit_coors = [[2, 2, 0, 0, 1, 1, 2, 1, 0],
                          [2, 0, 0, 2, 1, 2, 1, 0, 1]]
    ancs = []

    for i in range(9):
        ancs.append(Anchor(anchor_names[i],
                           i,
                           anchors_coors[:, i],
                           (anchors_unit_coors[0][i], anchors_unit_coors[1][i]),
                           False))

    main_anc = Anchor(anchor_names[main_N],
                      main_N,
                      anchors_coors[:, main_N],
                      (anchors_unit_coors[0][main_N], anchors_unit_coors[1][main_N]),
                      True)
    return AnchorsField(main_anc, ancs, anchors_coors)


if __name__ == "__main__":

    fname_syncs = "syncs_walls5.csv"
    fname_blinks = "blinks_walls5.csv"

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

    # plt.plot(syncs[:, 4], syncs[:, 2], 'r.')
    # plt.plot(syncs[:, 4], syncs[:, 3], 'y.')
    # plt.plot(blinks[:, 3], blinks[:, 2], 'g.')
    # plt.show()

    # make a starting time unit in sync file as time 0
    syncs[:, 4] *= 1000
    blinks[:, 3] *= 1000

    initial_time = syncs[0, 4]
    syncs[:, 4] -= initial_time
    blinks[:, 3] -= initial_time

    # Normalize: make all time units to be in us
    syncs[:, 2] = u2ms(syncs[:, 2])
    syncs[:, 3] = u2ms(syncs[:, 3])
    blinks[:, 2] = u2ms(blinks[:, 2])

    # # # # # # # # # # # # # # #
    # # # # # main part # # # # #
    # # # # # # # # # # # # # # #

    F = init_field(4)

    # For beginning let's ignore other nodes than 30
    node_name = 30

    last_i = syncs.shape[0]
    blinks_counter = 0
    for i in range(last_i):
        F.update_corrections(syncs[i, 0], syncs[i, 1], syncs[i, 2], syncs[i, 3], syncs[i, 4])
        anchors_active_idxs = [el is not None for el in F.l_time_corr]
        #
        while blinks[blinks_counter + 1, 3] <= syncs[i, 4]:
            if blinks[blinks_counter, 0] == node_name:
                pass
            blinks_counter += 1

        # if len([el is None for el in F.time_corrections]):
        #     continue
        print(F.l_time_corr)

# a = u2ms(syncs[:, 3])
# a = 1
# blinks_5 = blinks[blinks[:, 0] == 5]
# blinks_30 = blinks[blinks[:, 0] == 30]
# blinks_30 = blinks_30[blinks_30[:, 3].argsort()]
