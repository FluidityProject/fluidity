def val(X):
  if abs(X)-1.0 < 100.0:
    return 1.0
  if abs(X)-1.0 < 300.0:
    return 5.0
  else:
    return 10.0
