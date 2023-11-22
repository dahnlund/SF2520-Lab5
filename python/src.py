import numpy as np

class RelErrorTracker:
    def __init__(self, A, b):
        self.norms = np.array([])
        self.A = A
        self.b = b
    def callback(self,x):
        self.norms = np.append(self.norms, np.linalg.norm(self.A*x-self.b)/np.linalg.norm(self.b))
    """
    def callback2(self, norm):
        self.norms = np.append(self.norms, norm/np.linalg.norm(self.b))
    """