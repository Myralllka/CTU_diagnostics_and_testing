class Node:
    def __init__(self, number: int):
        # Assume only 9 anchors!
        self.number = number
        self.time_to_anchors = []

    def __str__(self):
        return "({})".format(self.number)

    def add_distance(self, anc_name, d, idx):
        pass

    def compute_position(self):
        pass
