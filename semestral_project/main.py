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

mpl.use('Qt5Agg')

anchor_names = [0x46, 0x4b, 0x45, 0x4c, 0x48, 0x4a, 0x47, 0x49, 0x44]
anchor_names = list(map(lambda x: int(x), anchor_names))
print(anchor_names)
anchors_coors = np.array([[-12.77, -12.77, -12.75, -6.77, -6.86, -6.87, -1.81, -1.92, -1.97],  # X coors
                          [2.75, -2.67, -8.05, 2.75, -2.67, -8.05, 2.75, -2.67, -8.05],  # Y coors
                          [3.13, 3.13, 3.13, 2.86, 2.86, 2.86, 2.6, 2.6, 2.6]])  # Z coors

c = 299792458  # m/s

fname_syncs = "syncs_walls5.csv"
fname_blinks = "blinks_walls5.csv"

syncs = np.genfromtxt(fname_syncs, delimiter=',')
blinks = np.genfromtxt(fname_blinks, delimiter=',')

new_array = [tuple(row) for row in blinks]
blinks = np.unique(new_array, axis=0)
blinks = blinks[blinks[:, 3].argsort()]

new_array = [tuple(row) for row in syncs]
syncs = np.unique(new_array, axis=0)
syncs = syncs[syncs[:, 4].argsort()]

if __name__ == "__main__":
    blinks_5 = blinks[blinks[:, 0] == 5]
    blinks_30 = blinks[blinks[:, 0] == 30]
    blinks_30 = blinks_30[blinks_30[:, 3].argsort()]

    N = 100
    blinks_30_70 = blinks_30[blinks_30[:, 1] == 72]
    a=0
    p = blinks[:, 3]
    print(np.min(p) - np.min(syncs[:, 4]))
    print(p[0] - syncs[0, 4])
    # print(np.max(p) - np.max(syncs[:, 4]))
    plt.plot(p, '.')
    plt.plot(syncs[:, 4], '.')
    plt.show()
    # for ()
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(blinks_30[:, 3], blinks_30[:, 2], blinks_30[:, 1])
# ax.set_xlabel('ts_rx')
# ax.set_ylabel('ts_tx')
# ax.set_zlabel('addr_rx')
# plt.show()
