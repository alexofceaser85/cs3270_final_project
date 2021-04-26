#!/usr/bin/env python3
'''
Base experiment that conducts the Label Propagation and prints results
'''
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler

__author__ = "Furichous Jones IV"
__version__ = "Spring 2021"

def calculate_scores(x_value, y_value):
    """
    calculates the accuracy scores
    """
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    scores = []
    for train, test in kfold.split(x_value):
        x_train, x_test = x_value[train], x_value[test]
        y_train, y_test = y_value[train], y_value[test]

        norm = MinMaxScaler().fit(x_train)
        x_train_norm = norm.transform(x_train)
        x_test_norm = norm.transform(x_test)

        yhat = create_label_propagation_model(
            x_train_norm, y_train, x_test_norm)
        # calculate score for test set
        score = accuracy_score(y_test, yhat)
        scores.append(score*100)
        # summarize score
    return scores


def create_label_propagation_model(x_train_norm, y_train, x_test_norm):
    """
    creates a new label propagation model
    """
    model = LabelPropagation(max_iter=1000, n_jobs=-1)
    # fit model on training dataset
    model.fit(x_train_norm, y_train)
    # make predictions on hold out test set
    return model.predict(x_test_norm)


def execute(x_value, y_value, experiment_name):
    """
    executes the label propagation and prints fold results and mean
    """
    print('____________________________________')
    print()
    print('Experiment ' + str(experiment_name))
    scores = calculate_scores(x_value, y_value)
    for score in scores:
        print('Accuracy: %.3f' % score)
    print('Mean Accuracy: %.3f (SD: %.3f)' % (np.mean(scores), np.std(scores)))
    print()
    print('____________________________________')
