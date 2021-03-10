#!/usr/bin/env python3

from numpy import concatenate
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
import pprint

__author__ = "Alexander Ayers"
__version__ = "Spring 2021"

def main():
  # define dataset
    X, y = make_classification(n_samples=1000, n_features=2, n_informative=2, n_redundant=0, random_state=1)
    print(y)
# split into train and test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.50, random_state=1, stratify=y)
# split train into labeled and unlabeled
    X_train_lab, X_test_unlab, y_train_lab, y_test_unlab = train_test_split(X_train, y_train, test_size=0.50, random_state=1, stratify=y_train)
# create the training dataset input
    X_train_mixed = concatenate((X_train_lab, X_test_unlab))
# create "no label" for unlabeled data
    nolabel = [-1 for _ in range(len(y_test_unlab))]
# recombine training dataset labels
    y_train_mixed = concatenate((y_train_lab, nolabel))
# define model
    model = LabelPropagation()
# fit model on training dataset
    model.fit(X_train_mixed, y_train_mixed)
# make predictions on hold out test set
    yhat = model.predict(X_test)
# calculate score for test set
    score = accuracy_score(y_test, yhat)
# summarize score
    print('Accuracy: %.3f' % (score*100))

if __name__ == "__main__":
    main()