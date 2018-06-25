# Wilson Disease Machine Learning ####
# WD Blood Neurologic vs Hepatic
# Written by Yihui Zhu

# coding: utf-8

## Import python packages
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import LeaveOneOut
from sklearn import preprocessing
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

# training set
# load methylation data for each DMRs at each sample in training set
A = pd.read_csv("WDN_vs_WDH_Blood_Training_Sample_Methylation.csv")

# load training set diagnosis data
B = pd.read_csv("WDN_vs_WDH_Training_Samples.csv")

# testing set
# load methylation data for each DMRs at each sample in testing set
C = pd.read_csv("WDN_vs_WDH_Blood_Test_Sample_Methylation.csv")

# load testing set diagnosis data
D = pd.read_csv("WDN_vs_WDH_Test_Samples.csv")

# process data
# Input Neurologic as 0 and Hepatic as 1
n_tr = 40
n_te = 10
m = len(A)

X_tr = np.zeros((n_tr,m))
y_tr = np.zeros(n_tr)
dat_tr = np.array(A)

for k in range(n_tr):
    X_tr[k,:] = dat_tr[:,k+3]

for k in range(n_tr):
    y_tr[k] = (B.WD_Phenotype[k] == 'Hepatic')

X_train = X_tr
y_train = y_tr
    
X_test = np.zeros((n_te,m))
y_test = np.zeros(n_te)
dat_te = np.array(C)

for k in range(n_te):
    X_test[k,:] = dat_te[:,k+3]

for k in range(n_te):
    y_test[k] = (D.WD_Phenotype[k] == 'Hepatic')

# AdaBoost
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

adb = AdaBoostClassifier()
adb.fit(X_train,y_train)

adb_pred = adb.predict(X_test)
adb_prob = adb.predict_proba(X_test)

# precision
1 - np.linalg.norm(adb_pred - y_test, 0) / len(y_test)

# list the testing group diagnosis value
y_test

# list the predicted diagnosis value by AdaBoost
adb_pred

# list the AdaBoost predicted value probability
adb_prob[:,1]

# feature selection on AdaBoost model
feature = np.array(adb.feature_importances_)
plt.plot(adb.feature_importances_)
plt.show()

# plot the feature selection value
from matplotlib import pyplot as plt  
import matplotlib 
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
feature = np.array(adb.feature_importances_)
fig = plt.plot(adb.feature_importances_)
plt.xlabel('DMRs', fontsize=18)
plt.ylabel('Feature Importance', fontsize=18)
fig_size = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = fig_size
fig_size[0] = 12
fig_size[1] = 6
plt.text(1000, 0.05, "Precision = 0.9", fontsize=22,ha='left',  wrap=True)
plt.savefig('Figure 7A.png', figsize=(12, 6),dpi=800)

## select the none zero 
## use ind to select out those vriables 
ind = np.where(feature != 0)[0]
A.loc[ind,:]

## Save the selected 44 DMRs in new variables, sel
sel = A.loc[ind,:]
np.array(sel)

import pandas as pd 
df = pd.DataFrame(sel)
df.to_csv("Selected_DMR.csv")

## Use the new set of DMRs to see how to the new model's precision
# Redo processing data
n_tr = 40
n_te = 10
m = len(ind)

X_tr = np.zeros((n_tr,m))
y_tr = np.zeros(n_tr)
dat_tr = np.array(A)

for k in range(n_tr):
    X_tr[k,:] = dat_tr[ind,k+3]

for k in range(n_tr):
    y_tr[k] = (B.WD_Phenotype[k] == 'Hepatic')

X_retrain = X_tr
y_train = y_tr
    
X_retest = np.zeros((n_te,m))
y_test = np.zeros(n_te)
dat_te = np.array(C)

for k in range(n_te):
    X_retest[k,:] = dat_te[ind,k+3]

for k in range(n_te):
    y_test[k] = (D.WD_Phenotype[k] == 'Hepatic')

from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

adb = AdaBoostClassifier()
adb.fit(X_retrain,y_train)

# precision
1 - np.linalg.norm(adb_pred - y_test, 0) / len(y_test)

## The model before and after the selection gives the same precision. 
