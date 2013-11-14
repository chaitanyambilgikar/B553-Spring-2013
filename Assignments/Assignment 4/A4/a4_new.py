from numpy import *
f = open("train.txt","r")
train_string = f.read()

z = train_string.split("\n")
cpts = {}
cpts[('d',0)] = 0.0
cpts[('d',1)] = 0.0
cpts[('p',0)] = 0.0
cpts[('p',1)] = 0.0
for string in z:
	if string[0] == '1':
		cpts[('d',0)] = cpts[('d',0)] + 1
	if string[2] == '1':
		cpts[('p',0)] = cpts[('p',0)] + 1
cpts[('d',0)]/=1000
cpts[('d',1)] = 1 - cpts[('d',0)]
cpts[('p',0)]/=1000
cpts[('p',1)] = 1 - cpts[('p',0)]
""" declaring a uniform distribution for all the cpt's of the variables """
cpts[('i',0)] = 0.5
cpts[('i',1)] = 0.5
cpts[('l',000)] = 1.0/8
cpts[('l',001)] = 1.0/8
print binary_repr(2)
""" cpts[('l',010)] = 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8)"""