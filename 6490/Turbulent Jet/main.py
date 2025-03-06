from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

def question1(u1,u2,x1,x2):
    u1_bar = np.mean(u1, axis=2)
    u2_bar = np.mean(u2, axis=2)
    u1_prime = u1 - u1_bar
    u2_prime = u2 - u2_bar





if(__name__ == "__main__"):
    file_path = "data_jet.mat"
    data = loadmat(file_path)
    u1 = np.array(data['u1'])
    u2 = np.array(data['u2'])
    x1 = np.array(data['x1'])[0,:]
    x2 = np.array(data['x2'])[:,0]
    print(u1.shape)
    print(x1.shape)
    print(x2.shape)

    question1(u1,u2,x1,x2)