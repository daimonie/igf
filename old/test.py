import numpy as np

gamma_left = np.array([[ 0.0102,  0.    ], [ 0.,      0.    ]])
ad_gf = np.array([[ 0.34114042-187.98045019j,  1.64280995 +38.98008398j], [ 1.64280995 +38.98008398j,  7.90760947  -8.44562354j]])
gamma_right = np.array([[ 0.,      0.    ], [ 0.,      0.0102]])
ret_gf = np.array([[ 0.34114042+187.98045019j,  1.64280995 -38.98008398j], [ 1.64280995 -38.98008398j,  7.90760947  +8.44562354j]])

print gamma_left
print ad_gf
print gamma_right
print ret_gf

left = np.dot(gamma_left, ad_gf)
right = np.dot(gamma_right, ret_gf)

tm = np.dot(left, right)

print "T(E0)=%2.3f" % np.real(np.trace(tm))