#!/usr/bin/env python3

from numpy import array
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler
import utils

def main():
    # define dataset
    print('____________________________________')
    print()
    print('Experiment 2')
    X = array(utils.parse_values("./interaction_data2.csv", [6]))
    y = array(utils.parse_is_interacted("./interaction_data2.csv"))
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    scores = []
    for train, test in kfold.split(X):
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]

        norm = MinMaxScaler().fit(X_train)
        X_train_norm = norm.transform(X_train)
        X_test_norm = norm.transform(X_test)

        model = LabelPropagation(max_iter = 1000, n_jobs = -1)
        # fit model on training dataset
        model.fit(X_train_norm, y_train)
        # make predictions on hold out test set
        yhat = model.predict(X_test_norm)
        # calculate score for test set
        score = accuracy_score(y_test, yhat)
        scores.append(score*100)
        # summarize score
        print('Accuracy: %.3f' % (score*100))
    print('Mean Accuracy: %.3f (SD: %.3f)' % (np.mean(scores), np.std(scores)))
    print('____________________________________')

if __name__ == "__main__":
    main()
