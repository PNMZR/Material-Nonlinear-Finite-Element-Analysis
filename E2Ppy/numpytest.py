# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:21:11 2020

@author: lenovo
"""
import numpy as np

test = np.array([[ 1,  2,  3,  4], 
                 [ 5,  6,  7,  8], 
                 [ 9, 10, 11, 12], 
                 [13, 14, 15, 16]])

print("<<<改变前>>>")
print(test)
#test[[1,3], :][:, [1,3]] = np.ones((2, 2))
test[np.ix_([0,1,3], [0, 1,3])] = np.ones((3, 3))
print("<<<改变后>>>")
print(test)

# print("<<<改变前>>>")
# print(test)
# test[[1,3], :] = np.array([[ 1,  1,  1,  1], 
#                  [ 1,  1,  1,  1]])
# print("<<<改变后>>>")
# print(test)