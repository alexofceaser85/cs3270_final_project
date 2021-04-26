#!/usr/bin/env python3
"""
Methods for computing sequence similarity of PPI.
"""
__author__ = "Furichous Jones IV"
__version__ = "Spring 2021"


def check_l_value(l_value, kappa, splice_site, length):
    """
    Returns true if valid l_value
    """
    return (l_value <= (round(kappa * abs(splice_site - l_value)) + 1)) and (l_value < length)

def acceptor_similarity(avec_first_half, bvec_first_half, a_length,
                        splice_site, length, degree, kappa):
    """
    Computes acceptor similarity.
    """
    sumk_acc = 0.0
    sum_acc = 0.0
    for k_value in range(1, degree+1):
        for i_value in range(0, a_length - k_value + 1):
            l_value = 1
            while check_l_value(l_value, kappa, splice_site, length):
                match = True
                j_value = i_value
                while (j_value < i_value + k_value) and match:
                    if j_value - l_value < 0:
                        match = False
                        break
                    match = (avec_first_half[j_value] ==
                             bvec_first_half[j_value - l_value])
                    j_value += 1
                if match:
                    delta = (1 / (2 * (l_value + 1)))
                    sumk_acc += (delta * (1 / a_length))
                l_value += 1

            l_value = 0
            while check_l_value(l_value, kappa, splice_site, length):
                match = True
                j_value = i_value
                while (j_value < i_value + k_value) and match:
                    if j_value + l_value >= a_length:
                        match = False
                        break
                    match = (avec_first_half[j_value] ==
                             bvec_first_half[j_value + l_value])
                    j_value += 1
                if match:
                    delta = (1 / (2 * (l_value + 1)))
                    sumk_acc += (delta * (1 / a_length))
                l_value += 1

        beta = (2 * (degree - k_value + 1)) / (degree * (degree + 1))
        sum_acc += sumk_acc * beta
    return sum_acc


def donor_similarity(avec_second_half, bvec_second_half, b_length, a_length,
                     splice_site, length, degree, kappa):
    """
    Computes donor similarity.
    """
    sumk_don = 0.0
    sum_don = 0.0
    for k_value in range(1, degree+1):
        for i_value in range(0, b_length - k_value + 1):
            l_value = 1
            while check_l_value(l_value, kappa, splice_site, length):
                match = True
                j_value = i_value
                while (j_value < i_value + k_value) and match:
                    if j_value - l_value < 0:
                        match = False
                        break
                    match = (avec_second_half[j_value]
                             == bvec_second_half[j_value - l_value])
                    j_value += 1
                if match:
                    delta = (1 / (2 * (l_value + 1)))
                    sumk_don += (delta * (1 / a_length))
                l_value += 1

            l_value = 0
            while check_l_value(l_value, kappa, splice_site, length):
                match = True
                j_value = i_value
                while (j_value < i_value + k_value) and match:
                    if j_value + l_value >= a_length:
                        match = False
                        break
                    match = (avec_second_half[j_value]
                             == bvec_second_half[j_value + l_value])
                    j_value += 1
                if match:
                    delta = (1 / (2 * (l_value + 1)))
                    sumk_don += (delta * (1 / a_length))
                l_value += 1

        beta = (2 * (degree - k_value + 1)) / (degree * (degree + 1))
        sum_don += sumk_don * beta
    return sum_don


def compute_similarity(a_string, b_string, degree, sigma, kappa):
    """
    Computes sequences between two protein sequences.
    """
    min_length = min(len(a_string), len(b_string))
    a_string = a_string[:min_length]
    b_string = b_string[:min_length]

    half_length = len(a_string) // 2

    length = 21

    avec_first_half = a_string[0: half_length]
    bvec_first_half = b_string[0: half_length]

    a_length = len(avec_first_half)
    b_length = len(bvec_first_half)
    splice_site = a_length // 2

    sum_acc = acceptor_similarity(
        avec_first_half, bvec_first_half, a_length, splice_site, length, degree, kappa)

    avec_second_half = a_string[half_length: len(a_string)]
    bvec_second_half = b_string[half_length: len(b_string)]

    a_length = len(avec_second_half)
    b_length = len(bvec_second_half)

    sum_don = donor_similarity(avec_second_half, bvec_second_half, b_length,
                               a_length, splice_site, length, degree, kappa)

    return sum_acc + sum_don


def main():
    """
    Main execution point for testing purposes.
    """
    avec = "MASEFKKKLFWRAVVAEFLATASIF"
    bvec = "TASIFKVKLFKKKLFWRASSAEFLA"
    degree = 33
    sigma = 1.0
    kappa = 0.37

    similarity = compute_similarity(avec, bvec, degree, sigma, kappa)
    print(similarity)


if __name__ == "__main__":
    main()
