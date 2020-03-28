import numpy as np
import sys

def generate_HMM_observation(num_obs, pi, T, E):
    def drawFrom(probs):
        return np.where(np.random.multinomial(1,probs) == 1)[0][0]

    obs = np.zeros(num_obs)
    states = np.zeros(num_obs)
    states[0] = drawFrom(pi)
    obs[0] = drawFrom(E[:, int(states[0])])
    for t in range(1,num_obs):
        states[t] = drawFrom(T[int(states[t-1]),:])
        obs[t] = drawFrom(E[:, int(states[t])])
    return obs, states


# number of states
N = int(sys.argv[1])

# number of observations
M = int(sys.argv[2])

true_pi = np.random.dirichlet(np.ones(M),size=1)
true_A = np.random.dirichlet(np.ones(M),size=M)
true_B = np.random.dirichlet(np.ones(N),size=M)

obs_seq, _ = generate_HMM_observation(50, true_pi, true_A, true_B)
with open('../input/observations.txt', 'w') as f:
    for item in obs_seq[:-1]:
        f.write("%i " % item)
    f.write("%i\n" % obs_seq[-1])

init_pi = np.random.dirichlet(np.ones(M),size=1)
init_A = np.random.dirichlet(np.ones(M),size=M)
init_B = np.random.dirichlet(np.ones(N),size=M)

with open('../input/sample.txt', 'w') as f:
    f.write(str(N) + " " + str(M) + "\n")
    f.write("\n")
    for item in init_pi[0][:-1]:
        f.write("%.4f " % item)
    f.write("%.4f" % init_pi[0][-1])
    f.write("\n")
    f.write("\n")

    for i in range(init_A.shape[0]):
        for item in init_A[i][:-1]:
            f.write("%.4f " % item)
        f.write("%.4f\n" % init_A[i][-1])
    f.write("\n")

    for i in range(init_B.shape[0]):
        for item in init_B[i][:-1]:
            f.write("%.4f " % item)
        f.write("%.4f\n" % init_B[i][-1])
