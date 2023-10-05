import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import svm
from time import process_time
import time
import seaborn as sns
from sklearn.model_selection import train_test_split


# Our dataset and targets
t1_start = process_time()

df = pd.read_csv('C:\\Users\\dtn228\\PycharmProjects\\SVM\\wind_training.csv')

X = df.drop('label', axis=1)
Y = df['label']


df_test = pd.read_csv('C:\\Users\\dtn228\\PycharmProjects\\SVM\\wind_test.csv')

X_test = df_test.drop('label', axis=1)
Y_test = df_test['label']


print(X)
print(Y)

# X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.20)
axs = plt.subplots(2, 1)

# figure number
fignum = 1

# fit the model
for kernel in ("linear", "poly", "rbf"):
    clf = svm.SVC(kernel="rbf", C=1000000, gamma='auto')
    clf.fit(X.values, Y.values)
    Y_pred_test = clf.predict(X_test.values)
    Y_pred_training = clf.predict(X.values)

    # plot the line, the points, and the nearest vectors to the plane
    plt.figure(fignum, figsize=(4, 3))
    plt.clf()

    #plt.scatter(
    #    clf.support_vectors_[:, 0],
    #    clf.support_vectors_[:, 1],
    #    s=80,
    #    facecolors="none",
    #    zorder=10,
    #    edgecolors="k",
    #)

    plt.scatter(df_test['w1'], df_test['w2'], c=Y_test, zorder=10, cmap=plt.cm.Paired, edgecolors="k")

    plt.axis("tight")
    x_min = 0
    x_max = 1
    y_min = 0
    y_max = 1

    XX, YY = np.mgrid[x_min:x_max:600j, y_min:y_max:600j]
    Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(XX.shape)
    plt.figure(fignum, figsize=(4, 3))
    plt.pcolormesh(XX, YY, Z > 0, cmap=plt.cm.Paired)
    plt.contour(
        XX,
        YY,
        Z,
        colors=["k", "k", "k"],
        linestyles=["--", "-", "--"],
        levels=[-0.5, 0, 0.5],
    )
    print(Y_pred_test)
    print(sum(abs(Y_pred_test-Y_test)))
    print(sum(abs(Y_pred_training - Y)))

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    plt.xticks(())
    plt.yticks(())

    csfont = {'fontname': 'Times New Roman', 'fontweight':'bold'}
    #plt.title("Kernelized SVM: test samples", **csfont, fontsize=22)
    plt.xlabel('Misclassification=14 out of 1000', **csfont, fontsize=20)
    #plt.xlabel(["first line \n second line"], **csfont, fontsize=15)
    #plt.ylabel("w2", **csfont, fontsize=20)

    fignum = fignum + 1
plt.show()

t1_stop = process_time()
curr_time = time.time()
print(curr_time)

print("Elapsed time:", t1_stop, t1_start)

print("Elapsed time during the whole program in seconds:",
      t1_stop - t1_start)