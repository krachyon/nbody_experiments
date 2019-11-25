import build_debug.nbody as nbody
import numpy as np

A=np.array([[1.,2],[1,2]])
B=np.array([1.,1.])

para = nbody.SimulationParameters()


nbody.EulerSimulation(para,A,A,B).step(1)
