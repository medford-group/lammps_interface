from lammps_interface.tools import make_fingerprint_matrix
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from ase.io import read

traj = read('interface_tip3p.traj', index='::10')

g2_etas=[0.005, 4.0, 20.0, 80.0]
g2_rs_s=[0] * 4
g4_etas=[0.005]
g4_zetas=[1.,4.]
g4_gammas=[1, -1]
cutoff=6.5

array_dict = make_fingerprint_matrix(traj,
        [g2_etas, g2_rs_s, g4_etas, cutoff, g4_zetas, g4_gammas],
        elements=['H','O'])


pca = PCA(n_components=2)
X_new = pca.fit_transform(array_dict['O'])


plt.scatter(X_new[:,0],X_new[:,1], alpha=0.05)

plt.title('PCA Decomposition of Fingerprints')
plt.ylabel('PC2')
plt.xlabel('PC1')
plt.savefig('PCA_O.png')
plt.show()
