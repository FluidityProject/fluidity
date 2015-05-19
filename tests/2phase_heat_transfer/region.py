def val(x):
	tol=1.0e-05
	if (abs(x-0.0)<tol):
		region=1
	elif (abs(x-400.0)<tol):
		region=2
	else:
		region=3
	return region

