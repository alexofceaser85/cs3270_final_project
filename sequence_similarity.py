import math

def compute_similarity(avec, bvec, degree, sigma, kappa):
	min_length = min(len(avec),len(bvec))
	avec = avec[:min_length]
	bvec = bvec[:min_length]
	
	half_length = len(avec) // 2

	L = 21
	'''EXTRACT ACCEPTORS'''

	avec_first_half = avec[0 : half_length]
	bvec_first_half = bvec[0 : half_length]

	alen = len(avec_first_half)
	blen = len(bvec_first_half)
	splice_site = alen // 2

	''' A C C E P T O R     S I M I L A R I T Y '''

	sumk_acc = 0.0
	sum_acc = 0.0

	for k in range(1, degree+1):
		for i in range(0, alen - k + 1):
			l = 1
			while (l <= (round(kappa * abs(splice_site - l)) + 1)) and (l < L):
				match = True
				j = i
				while (j < i + k) and match:
					if (j - l < 0):
						match = False
						break
					else:
						match = (avec_first_half[j] == bvec_first_half[j - l])
					j += 1
				if (match):
					delta = (1 / (2 * (l + 1)))
					sumk_acc += (delta * (1 / alen))
				l += 1

			l = 0
			while (l <= (round(kappa * abs(splice_site - l)) + 1)) and (l < L):
				match = True
				j = i
				while (j < i + k) and match:
					if (j + l >= alen):
						match = False
						break
					else:
						
						match = (avec_first_half[j] == bvec_first_half[j + l])
					j += 1
				if (match):
					delta = (1 / (2 * (l + 1)))
					sumk_acc += (delta * (1 / alen))
				l += 1
			
		beta = (2 * (degree - k + 1)) / (degree * (degree + 1))
		sum_acc += sumk_acc * beta

	''' EXTRACT DONORS '''
	avec_second_half = avec[half_length : len(avec)]
	bvec_second_half = bvec[half_length : len(bvec)]

	alen = len(avec_second_half)
	blen = len(bvec_second_half)

	''' D O N O R     S I M I L A R I T Y '''
	sumk_don = 0.0
	sum_don = 0.0

	for k in range(1, degree+1):
		for i in range(0, blen - k + 1):
			l = 1
			while (l <= (round(kappa * abs(splice_site - l)) + 1)) and (l < L):
				match = True
				j = i
				while (j < i + k) and match:
					if (j - l < 0):
						match = False
						break
					else:
						match = (avec_second_half[j] == bvec_second_half[j - l])
					j += 1
				if (match):
					delta = (1 / (2 * (l + 1)))
					sumk_don += (delta * (1 / alen))
				l += 1

			l = 0
			while (l <= (round(kappa * abs(splice_site - l)) + 1)) and (l < L):
				match = True
				j = i
				while (j < i + k) and match:
					if (j + l >= alen):
						match = False
						break
					else:
						match = (avec_second_half[j] == bvec_second_half[j + l])
					j += 1
				if (match):
					delta = (1 / (2 * (l + 1)))
					sumk_don += (delta * (1 / alen))
				l += 1

		beta = (2 * (degree - k + 1)) / (degree * (degree + 1))
		sum_don += sumk_don * beta

	return sum_acc + sum_don


def main():
	avec = "MASEFKKKLFWRAVVAEFLATASIF"
	bvec = "TASIFKVKLFKKKLFWRASSAEFLA"
	degree = 33
	sigma = 1.0
	kappa = 0.37

	similarity = compute_similarity(avec, bvec, degree, sigma, kappa)
	print(similarity)


if __name__ == "__main__":
    main()
