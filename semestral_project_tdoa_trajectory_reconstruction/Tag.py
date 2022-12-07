class Node:
    def __init__(self, number: int):
        # Assume only 9 anchors!
        # TODO: use for filtering
        self.number = number

    def __str__(self):
        return "({})".format(self.number)
