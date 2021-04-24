#!/usr/bin/env python3

"""
Experiment 2 that utilizes the Jaccard and L3 score
and Label Propagation to compare the PPIs.
"""

from numpy import array
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler
import utils

def create_label_propagation_model(x_train_norm, y_train, x_test_norm):
    """
    creates a new label propagation model
    """
    model = LabelPropagation(max_iter=1000, n_jobs=-1)
    # fit model on training dataset
    model.fit(x_train_norm, y_train)
    # make predictions on hold out test set
    return model.predict(x_test_norm)

def main():
    """
    Execution point to the experiment.
    """
    # define dataset
    print('____________________________________')
    print()
    print('Experiment 2')
    x_value = array(utils.parse_values("./data/interaction_data2.csv", [5, 6]))
    y_value = array(utils.parse_is_interacted("./data/interaction_data2.csv"))
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    scores = []
    for train, test in kfold.split(x_value):
        x_train, x_test = x_value[train], x_value[test]
        y_train, y_test = y_value[train], y_value[test]

        norm = MinMaxScaler().fit(x_train)
        x_train_norm = norm.transform(x_train)
        x_test_norm = norm.transform(x_test)

        yhat = create_label_propagation_model(x_train_norm, y_train, x_test_norm)
        # calculate score for test set
        score = accuracy_score(y_test, yhat)
        scores.append(score*100)
        # summarize score
        print('Accuracy: %.3f' % (score*100))
    print('Mean Accuracy: %.3f (SD: %.3f)' % (np.mean(scores), np.std(scores)))
    print('____________________________________')


if __name__ == "__main__":
    main()
