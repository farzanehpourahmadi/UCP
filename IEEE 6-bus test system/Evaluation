import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import svm, metrics
from time import process_time
import time
from sklearn import preprocessing
from sklearn.metrics.pairwise import pairwise_kernels

from sklearn.utils import check_array

import seaborn as sns
from sklearn.model_selection import train_test_split

error = []
hl = []
hl_tr = []
hl_t = []
score_total = []
hinge_loss_m_total = []
error_tr = []
loss_test_total = []
# Our dataset and targets
t1_start = process_time()

df = pd.read_csv('C:\\Users\\dtn228\\PycharmProjects\\SVM\\wind_training_UC_different.csv')
df_test = pd.read_csv('C:\\Users\\dtn228\\PycharmProjects\\SVM\\wind_test_UC_different.csv')

X = df.values[:, [0, 1]]
X_test = df_test.values[:, [0, 1]]
lab = preprocessing.LabelEncoder()

for i in range(72):
    Y1 = df.values[:, i+2]
    Y_test1 = df_test.values[:, i + 2]
    Y = lab.fit_transform(Y1)
    Y_test = lab.fit_transform(Y_test1)
    #print(X_test)
    #print(Y_test)
    #print(i)
    #print(X)
    #print(Y)
    clf = svm.SVC(kernel="rbf", C=1000000, gamma='scale', decision_function_shape='ovr')
    # fit the classifier to your data
    clf.fit(X, Y)
    Y_pred_test = clf.predict(X_test)
    Y_pred_training = clf.predict(X)


    # get the predicted scores for your test data
    clf.decision_function_shape = "ovr"
    Y_pred_test_scores = clf.decision_function(X_test)
    Y_pred_test_scores = Y_pred_test_scores.reshape(1000)
    Y_pred_scores = clf.decision_function(X)

    #print(Y_pred_test_scores)
    #print(Y_test.shape)
    #print(Y_pred_test_scores.shape)
    loss_test = metrics.hinge_loss(Y_test, Y_pred_test_scores)
    loss_training = metrics.hinge_loss(Y, Y_pred_scores)
    score1 = clf.score(X, Y)
    print(i)
    print(score1)
    print(loss_test)
    print(loss_training)
    #print(Y_pred_test_scores)


    error.append(sum(abs(Y_pred_test - Y_test)))
    error_tr.append(sum(abs(Y_pred_training - Y)))

#    hl_tr.append(hinge_loss(Y_pred_training, Y))
    loss_test_total.append(loss_test)
    score_total.append(score1)
    #loss_test_total.append(loss_test)

    t1_stop = process_time()
    curr_time = time.time()


print(error)
print(error_tr)

print(score_total)

print(loss_test_total)


#print(hinge_loss_m_total)
