import numpy as np

def find_err(a, b):
    assert(a.shape == b.shape)
    diff = np.abs(a-b)
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)
    return mean_diff, max_diff

errors_mean = []
errors_max = []
for i in range(20):
    with open('../output/output_' + str(i) + '.txt', 'r') as f:
        a1 = np.loadtxt(f, max_rows=2)
        b1 = np.loadtxt(f, max_rows=2)
        pi1 = np.loadtxt(f, max_rows=1)
        g1 = np.loadtxt(f, max_rows=1)

    with open('../output/output_cpp_' + str(i) + '.txt', 'r') as f:
        a2 = np.loadtxt(f, max_rows=2)
        b2 = np.loadtxt(f, max_rows=2)
        pi2 = np.loadtxt(f, max_rows=1)
        g2 = np.loadtxt(f, max_rows=1)

    mean_a, max_a = find_err(a1, a2)    
    mean_b, max_b = find_err(b1, b2)    
    mean_pi, max_pi = find_err(pi1, pi2)    
    mean_g, max_g = find_err(g1, g2)    

    errors_mean.append([mean_a, mean_b, mean_pi, mean_g])
    errors_max.append([max_a, max_b, max_pi, max_g])

print('Mean Error: ')
print(errors_mean)
print('\n')

print('Max Error: ')
print(errors_max)
    
