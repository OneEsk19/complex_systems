# import packages
import csv
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from sklearn.model_selection import train_test_split
from irf import irf_utils, irf_jupyter_utils
from irf.ensemble import RandomForestClassifierWithWeights


# load data
with open('../../data/sample_d6_d11.csv', 'r') as f:
    idx = np.array([r for r in csv.reader(f, delimiter = ',')])
    # idx  = idx[(idx.REF == 'Control')&(idx.REF == "10x")]
with open('../../data/date_d6_d11.csv', 'r') as f:
    dmx = np.array([r for r in csv.reader(f, delimiter = ',')])

# prep
X = dmx[1:,1:].astype(float).T
# assigning a nurmerical variable to the string classifiers
# 0 for cntrol, 1 for compared
# biclass
y = np.array(itemgetter(*idx[1:,2])({'D11': 0, 'D06': 1}))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3)


# iRF
irfres = irf_utils.run_iRF(
    X_train = X_train, X_test = X_test,
    y_train = y_train, y_test = y_test,
    rf = RandomForestClassifierWithWeights(n_estimators = 100),
    K = 5, # number of iteration. This is recommended value by dev
    B = 100, # The number of bootstrap samples. Play around with this, see what changes
    M = 137, # number of trees (RIT) to build. Look into the effects of this parameter
    max_depth = 5,
)

rf_weights, K_iter_rf_data, rf_bootstrap_output, rit_bootstrap_output, stability_score = irfres


# feature importance
# i probably don't need to worry about the rest of this, i.e modifying etc
fids = dmx[1:,0]

iteration = 'rf_iter5'
impt = K_iter_rf_data[iteration]['feature_importances']
impt_std = K_iter_rf_data[iteration]['feature_importances_std']
impt_rank_idx = K_iter_rf_data[iteration]['feature_importances_rank_idx']

impt, impt_std = impt[impt_rank_idx], impt_std[impt_rank_idx]
impt, impt_std, impt_rank_idx = impt[impt > 0], impt_std[impt > 0], impt_rank_idx[impt > 0]

xids = np.arange(impt.shape[0])
plt.figure()
plt.bar(
    xids, impt, yerr = impt_std, align = 'center'
)
plt.xticks(xids, fids[impt_rank_idx], rotation = 'vertical')
plt.xlim([-1, xids.shape[0]])
plt.show()

for f,i in zip(fids[impt_rank_idx],impt): 
    print(f'feature [{f}] with importance {i:.3f}')

    
# interaction stability scores
for i,s in sorted(stability_score.items(), key = lambda x: -x[1]):
    print(f'interaction [{i}] with stability score {s:.3f}')

irf_jupyter_utils._get_histogram(
    {k: v for k,v in stability_score.items() if v > 0}, 
    sort = True,
)
