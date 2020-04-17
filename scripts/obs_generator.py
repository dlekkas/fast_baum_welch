import numpy as np
import sys

from HiddenMarkovModel import HiddenMarkovModel

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


# number of observations
N = int(sys.argv[1])

# number of states
M = int(sys.argv[2])
'''

true_pi = np.array([0.5, 0.5])

true_A = np.array([[0.85, 0.15],
                  [0.12, 0.88]])

true_B = np.array([[0.8, 0.1, 0.1],
                   [0.0, 0.0, 1.0]])
'''
num_obs = 50
iterations=20
diff_list = []
for it in range(iterations):

    print("Iteration: ", it)

    true_pi = np.random.dirichlet(np.ones(M),size=1)
    true_A = np.random.dirichlet(np.ones(M),size=M)
    true_B = np.random.dirichlet(np.ones(N),size=M)

    obs_seq, obs_states = generate_HMM_observation(num_obs, true_pi[0], true_A, np.transpose(true_B))

    init_pi = np.random.dirichlet(np.ones(M),size=1)
    init_A = np.random.dirichlet(np.ones(M),size=M)
    init_B = np.random.dirichlet(np.ones(N),size=M)

    model = HiddenMarkovModel(init_A, np.transpose(init_B), init_pi[0], epsilon=0.0001, maxStep=5)
    trans0, transition, emission, c = model.run_Baum_Welch_EM(obs_seq, summary=False, monitor_state_1=True)

    with open('../input/observations_' + str(it) + '.txt', 'w') as f:
        for item in obs_seq[:-1]:
            f.write("%i " % item)
        f.write("%i\n" % obs_seq[-1])

    with open('../input/sample_' + str(it) + '.txt', 'w') as f:
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

    with open('../output/output_' + str(it) + '.txt', 'a') as f:
        np.savetxt(f, transition, fmt="%.6f")
        f.write("\n")

        np.savetxt(f, np.transpose(emission), fmt="%.6f") 
        f.write("\n")

        np.savetxt(f, trans0, fmt="%.6f", newline = ' ')    
        f.write("\n")
        f.write("\n")

        pred = 1- (model.state_summary[-2] > 0.5)
        np.savetxt(f, pred, fmt="%d", newline=' ')    

    '''
    model = HiddenMarkovModel(true_A, np.transpose(true_B), true_pi[0])
    states_seq, state_prob = model.run_viterbi(obs_seq, summary=True)

    print("Viterbi states:")
    print(states_seq)

    print("True states:")
    print(obs_states)

    diff_array = np.array(states_seq) - np.array(obs_states)
    diff = np.sum(np.abs(diff_array))/num_obs
    diff_list.append(diff)

print(diff_list)
diff_array = np.array(diff_list)
print("Mean: ")
print(np.mean(diff_array))

print("Deviation: ")
print(np.std(diff_array))
'''
'''
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

init_pi = np.array([0.5, 0.5])
init_A = np.array([[0.5, 0.5], [0.5, 0.5]])
init_B = np.array([[0.11111111,0.33333333,0.55555556], [0.16666667, 0.33333333, 0.5]])
obs_seq = [0 ,1 ,2, 2, 2, 2, 2, 2, 2,2]

#obs_seq = [0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1, 0, 2, 1, 2, 0, 2, 0, 1, 2, 1, 2, 0, 2, 0, 2, 2, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 2, 2, 2, 0, 2, 2, 2, 1, 2, 0, 1, 0, 0, 2, 1, 2, 1, 1, 1, 0, 2, 0, 0, 1, 1, 2, 0, 1, 2, 0, 1, 0, 2, 1, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 0, 2, 2, 1, 0, 0, 1, 1, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 0, 2, 0, 1, 2, 0, 1, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 1, 2, 1, 0, 2, 2, 1, 2, 2, 2, 1, 0, 1, 2, 2, 2, 1, 0, 1, 0, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0, 2, 0, 1, 1, 2, 0, 0, 2, 2, 2, 1, 1, 0, 0, 1, 2, 1, 2, 1, 0, 2, 0, 2, 2, 0, 0, 0, 1, 0, 1, 1, 1, 2, 2, 0, 1, 2, 2, 2, 0, 1, 1, 2, 2, 0, 1, 2, 2, 2, 2, 2, 2, 0, 1, 2, 2, 0, 2, 0, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 2, 0, 0, 0, 2, 2, 1, 1, 2, 1, 0, 2, 1, 1, 1, 0, 1, 2, 1, 2, 1, 2, 2, 2, 0, 2, 0, 0, 2, 2, 2, 2, 2, 2, 1, 0, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 0, 1, 2, 0, 1, 2, 1, 2, 0, 2, 1, 0, 2, 2, 0, 2, 2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 1, 2, 0, 2, 1, 2, 2, 2, 1, 2, 2, 2, 0, 0, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 0, 2, 2, 1, 2, 2, 2, 2, 1, 2, 0, 2, 1, 2, 2, 0, 1, 0, 1, 2, 1, 0, 2, 2, 2, 1, 0, 1, 0, 2, 1, 2, 2, 2, 0, 2, 1, 2, 2, 0, 1, 2, 0, 0, 1, 0, 1, 1, 1, 2, 1, 0, 1, 2, 1, 2, 2, 0, 0, 0, 2, 1, 1, 2, 2, 1, 2 ]
model = HiddenMarkovModel(init_A, np.transpose(init_B), init_pi, epsilon=0.0001, maxStep=5)
trans0, transition, emission, c = model.run_Baum_Welch_EM(obs_seq, summary=False, monitor_state_1=True)

#pred = (1-model.state_summary[-2]) > 0.5

print("Transition Matrix: ")
print(transition)
print()
print("Emission Matrix: ")
print(emission)
print()
print("Reached Convergence: ")
print(c)

#print(model.state_summary)
pred = 1- (model.state_summary[-2] > 0.5)
print(pred)

states_seq, state_prob = model.run_viterbi(obs_seq, summary=True)
print("Viterbi results for the model (tf) : ")
print(states_seq)
#print(state_prob)
# our case:

obs_seq = [0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0]
pi = np.array([0.5, 0.5])
A = np.array([[0.85, 0.15], [0.12, 0.88]])
B = np.array([[0.8, 0.1, 0.1], [0.0, 0.0, 1.0]])
model = HiddenMarkovModel(A, np.transpose(B), pi, maxStep=20)
states_seq, state_prob = model.run_viterbi(obs_seq, summary=True)
print("Viterbi results for our model: ")
print(states_seq)
'''




