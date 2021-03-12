#!/usr/bin/env python3

from numpy import concatenate
from numpy import array
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
import pprint
import utils

__author__ = "Alexander Ayers"
__version__ = "Spring 2021"

def main():
  # define dataset
    X = array(utils.parse_interaction_score_and_similarity("./interaction_data.csv"))
    y = array(utils.parse_is_interacted("./interaction_data.csv"))
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    for train, test in kfold.split(X):
        
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]
        model = LabelPropagation(max_iter = 1000, n_jobs = -1)
    # fit model on training dataset
        model.fit(X_train, y_train)
    # make predictions on hold out test set
        yhat = model.predict(X_test)
    # calculate score for test set
        score = accuracy_score(y_test, yhat)
    # summarize score
        print('Accuracy: %.3f' % (score*100))

if __name__ == "__main__":
    main()
