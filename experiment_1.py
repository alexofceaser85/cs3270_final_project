#!/usr/bin/env python3
"""
Experiment 1 that utilizes the Sequence Similarity
and Label Propagation to compare the PPIs.
"""
from numpy import array
import utils
import base_experiment

__author__ = "Alexander Ayers"
__version__ = "Spring 2021"

def main():
    """
    Execution point to the experiment.
    """
    # define dataset
    x_values = array(utils.parse_values("./data/interaction_data2.csv", [3]))
    y_values = array(utils.parse_is_interacted("./data/interaction_data2.csv"))

    base_experiment.execute(x_values, y_values, "1")


if __name__ == "__main__":
    main()
