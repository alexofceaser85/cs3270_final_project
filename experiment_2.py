#!/usr/bin/env python3

"""
Experiment 2 that utilizes the Jaccard and L3 score
and Label Propagation to compare the PPIs.
"""

from numpy import array
import utils
import base_experiment

def main():
    """
    Execution point to the experiment.
    """
    # define dataset
    x_values = array(utils.parse_values("./data/interaction_data2.csv", [5, 6]))
    y_values = array(utils.parse_is_interacted("./data/interaction_data2.csv"))

    base_experiment.execute(x_values, y_values, "2")


if __name__ == "__main__":
    main()
