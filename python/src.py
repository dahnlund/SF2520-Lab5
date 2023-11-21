import numpy as np

class RelErrorTracker:
    def __init__(self, A, b):
        self.norms = np.array([])
        self.A = A
        self.b = b
    def callback1(self,x):
        self.norms = np.append(self.norms, np.linalg.norm(self.A*x-self.b))
    def callback2(self, residual_norm):
        self.norms = np.append(self.norms, residual_norm)