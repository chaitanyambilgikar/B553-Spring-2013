
P_I = {}
P_G = {}
P_L = {}
P_S = {}
P_W = {}
P_A = {}



	


f = open("train.txt", "r")
""" read the training data in a String"""
trainString=f.read();

""" Z contains all the Z_i as a list """
Z = trainString.split("\n")

#D I P G L S W A

#Computing Probabilities for D and P
count=0.0;
countd =0.0
countp =0.0

for x in Z:
	y = x.split(" ");
	#print (x[0])
	countd = countd+int(x[0], base=10);
	countp = countp+int(x[0], base=10);
	count= count +1
	
#print(countd)
#print(countp)
#print(count)

#print("Probability of D =1",countd/count)


""" Making CPT for all the Nodes except D & P """

theta_d =[1-countd/count, countd/count]    							#{'0':countd/count, '1':1-countd/count}
theta_p =[1-countp/count, countp/count] 							#{'0':countp/count, '1':1-countp/count}
theta_i =[0.5,0.5]

#order L P G
theta_l =[1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8]
theta_g =[1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8]
theta_s =[1.0/4, 1.0/4, 1.0/4, 1.0/4]
theta_w =[1.0/4, 1.0/4, 1.0/4, 1.0/4]

#order G L S W A
theta_a =[1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32]









def OnlyImissing(Z_j,index ):
	
	"""
	Over here willll we do variable elimination indirectly
	i.e I wont use variable elimination algorithm 
	but construct simple formulas for variable elimination
	P(I) = Summation over all value(P_Algsw*P_Lgp*P_Gdi*P_Si*P_Wi*P_P*P_D)
	But over here all values except I are observed
	So we don't use variable elimination
	"""
	
	
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)] # observed probability for probvability of A giveb L,G,S,W
	
	Pl_pg=theta_l(int(Z_j[4]+Z_j[2]+Z_j[3],2)); # observed probability for Probability of L given P and G
	
	Pd=theta_d[int(Z_j[0],2)]; #observed probability of D
	
	Pp=theta_p[int(Z_j[2],2)]; #observed probability of P
	
	
	#Now the probabilities P_Gdi = theta_g *P_Si=theta_s *P_Wi = theta_w will will be calculated using variable elimintaion
	
	zero = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_w[int((Z_j[6]+str(0)),2)]*theta_s[int((Z_j[5]+str(0)),2)]
	one = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_w[int((Z_j[6]+str(1)),2)]*theta_s[int((Z_j[5]+str(1)),2)]
	
	P_I['index']=[zero,one];
	
	
	




def IandW(Z_j ):
	
	
	
	"""
	will have to marginalize this 
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)] # observed probability for probvability of A giveb L,G,S,W"""
	
	
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(0)+Z_j[7],2)] + theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(1)+Z_j[7],2)] 
	
	Pl_pg=theta_l(int(Z_j[4]+Z_j[2]+Z_j[3],2)); # observed probability for Probability of L given P and G
	
	Pd=theta_d[int(Z_j[0],2)]; #observed probability of D
	
	Pp=theta_p[int(Z_j[2],2)]; #observed probability of P
	
	
	#Now the probabilities P_Gdi = theta_g *P_Si=theta_s *P_Wi = theta_w will will be calculated using variable elimintaion
	#Calculating for I
	Izero = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*(theta_w[int((str(0)+str(0)),2)]+theta_w[int((str(1)+str(0)),2)])*theta_s[int((Z_j[5]+str(0)),2)]
	Ione = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*(theta_w[int((str(0)+str(1)),2)]+theta_w[int((str(1)+str(1)),2)])*theta_s[int((Z_j[5]+str(1)),2)]
	
	P_I['index']=[Izero,Ione];
	
	#Calculating for W
	Wzz= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(0)+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_s[int((Z_j[5]+str(0)),2)*theta_i[int(str(0))]
	Wzo= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(0)+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_s[int((Z_j[5]+str(1)),2)*theta_i[int(str(1))]
	Woz= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(1)+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_s[int((Z_j[5]+str(0)),2)*theta_i[int(str(0))]
	Woo= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+str(1)+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_s[int((Z_j[5]+str(1)),2)*theta_i[int(str(1))]
	
	P_W['index']=[Wzz,Wzo,Woz,Woo]



def IandS(Z_j ):
	
	
	
	"""
	will have to marginalize this 
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)] # observed probability for probvability of A giveb L,G,S,W"""
	
	
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+str(0)+Z_j[6]+Z_j[7],2)] + theta_a[int(Z_j[3]+Z_j[4]+str(1)+Z_j[6]+Z_j[7],2)] 
	
	Pl_pg=theta_l(int(Z_j[4]+Z_j[2]+Z_j[3],2)); # observed probability for Probability of L given P and G
	
	Pd=theta_d[int(Z_j[0],2)]; #observed probability of D
	
	Pp=theta_p[int(Z_j[2],2)]; #observed probability of P
	
	
	#Now the probabilities P_Gdi = theta_g *P_Si=theta_s *P_Wi = theta_w will will be calculated using variable elimintaion
	#Calculating for I
	Izero = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*(theta_s[int((str(0)+str(0)),2)]+theta_s[int((str(1)+str(0)),2)])*theta_w[int((Z_j[6]+str(0)),2)]
	Ione = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*(theta_s[int((str(0)+str(1)),2)]+theta_s[int((str(1)+str(1)),2)])*theta_w[int((Z_j[6]+str(1)),2)]
	
	P_I['index']=[Izero,Ione];
	
	#Calculating for S
	Szz= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+str(0)+Z_j[6]+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_w[int((Z_j[6]+str(0)),2)*theta_i[int(str(0))]
	Szo= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+str(0)+Z_j[6]+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_w[int((Z_j[6]+str(1)),2)*theta_i[int(str(1))]
	Soz= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+str(1)+Z_j[6]+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_w[int((Z_j[6]+str(0)),2)*theta_i[int(str(0))]
	Soo= Pp*Pd*Pl_pg*theta_a[int(Z_j[3]+Z_j[4]+str(1)+Z_j[6]+Z_j[7],2)]*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_w[int((Z_j[6]+str(1)),2)*theta_i[int(str(1))]
	
	P_S['index']=[Szz,Szo,Soz,Woo]


def IandL(Z_j ):
	
	#D I P G L S W A    1,3,4,5,6 can be x
	#0 1 2 3 4 5 6 7
	
	"""
	will have to marginalize this 
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)] # observed probability for probvability of A giveb L,G,S,W"""
	
	
	Pa_glsw= theta_a[int(Z_j[3]+str(0)+Z_j[5]+Z_j[6]+Z_j[7],2)] + theta_a[int(Z_j[3]+str(1)+Z_j[5]+Z_j[6]+Z_j[7],2)] 
	
	Pl_pg=theta_l(int(str(0)+Z_j[2]+Z_j[3],2))+theta_l(int(str(1)+Z_j[2]+Z_j[3],2)); # observed probability for Probability of L given P and G
	
	Pd=theta_d[int(Z_j[0],2)]; #observed probability of D
	
	Pp=theta_p[int(Z_j[2],2)]; #observed probability of P
	
	
	#Now the probabilities P_Gdi = theta_g *P_Si=theta_s *P_Wi = theta_w will will be calculated using variable elimintaion
	#Calculating for I
	Izero = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_w[int((Z_j[6]+str(0)),2)]*theta_s[int((Z_j[5]+str(0)),2)]
	Ione = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_w[int((Z_j[6]+str(1)),2)]*theta_s[int((Z_j[5]+str(1)),2)]
	
	P_I['index']=[Izero,Ione];
	
	#Calculating for Lpg
	
	Lzero=Pp*Pd*theta_a[int(Z_j[3]+str(0)+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(str(0)+Z_j[2]+Z_j[3],2))*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*(theta_w[int((Z_j[6]+str(0)),2)]+int((Z_j[6]+str(1)),2)])*(theta_s[int((Z_j[5]+str(0)),2)]+theta_s[int((Z_j[5]+str(1)),2)])*(theta_i[int(str(1))]+theta_i[int(str(0))])
	Lone=Pp*Pd*theta_a[int(Z_j[3]+str(1)+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(str(1)+Z_j[2]+Z_j[3],2))*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*(theta_w[int((Z_j[6]+str(0)),2)]+int((Z_j[6]+str(1)),2)])*(theta_s[int((Z_j[5]+str(0)),2)]+theta_s[int((Z_j[5]+str(1)),2)])*(theta_i[int(str(1))]+theta_i[int(str(0))]
	
	
	P_L['index']=[Lzero,Lone]


def IandG(Z_j ):
	
	#D I P G L S W A    1,3,4,5,6 can be x
	#0 1 2 3 4 5 6 7
	
	"""
	will have to marginalize this 
	Pa_glsw= theta_a[int(Z_j[3]+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)] # observed probability for probvability of A giveb L,G,S,W"""
	
	
	Pa_glsw= theta_a[int(str(0)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]+theta_a[int(str(1)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]
	
	Pl_pg=theta_l(int(Z_j[4]+Z_j[2]+str(0),2)) + theta_l(int(Z_j[4]+Z_j[2]+str(1),2)) # observed probability for Probability of L given P and G
	
	Pd=theta_d[int(Z_j[0],2)]; #observed probability of D
	
	Pp=theta_p[int(Z_j[2],2)]; #observed probability of P
	
	
	#Now the probabilities P_Gdi = theta_g *P_Si=theta_s *P_Wi = theta_w will will be calculated using variable elimintaion
	#Calculating for I
	Izero = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(0)),2)]*theta_w[int((Z_j[6]+str(0)),2)]*theta_s[int((Z_j[5]+str(0)),2)]
	Ione = Pp*Pd*Pa_glsw*Pl_pg*theta_g[int((Z_j[3]+Z_j[0]+str(1)),2)]*theta_w[int((Z_j[6]+str(1)),2)]*theta_s[int((Z_j[5]+str(1)),2)]
	
	P_I['index']=[Izero,Ione];
	
	#Calculating for Lpg
	
	
	
	Gzz=Pp*Pd*theta_a[int(str(0)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(Z_j[4]+Z_j[2]+str(0),2))*(theta_w[int((Z_j[6]+str(0)),2)])*(theta_s[int((Z_j[5]+str(0)),2)]),2)])*(theta_i[int(str(0))])
	Gzo=Pp*Pd*theta_a[int(str(0)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(Z_j[4]+Z_j[2]+str(0),2))*(theta_w[int((Z_j[6]+str(1)),2)])*(theta_s[int((Z_j[5]+str(1)),2)]),2)])*(theta_i[int(str(1))])
	Goz=Pp*Pd*theta_a[int(str(1)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(Z_j[4]+Z_j[2]+str(1),2))*(theta_w[int((Z_j[6]+str(0)),2)])*(theta_s[int((Z_j[5]+str(0)),2)]),2)])*(theta_i[int(str(0))])
	Goo=Pp*Pd*theta_a[int(str(1)+Z_j[4]+Z_j[5]+Z_j[6]+Z_j[7],2)]*theta_l(int(Z_j[4]+Z_j[2]+str(1),2))*(theta_w[int((Z_j[6]+str(1)),2)])*(theta_s[int((Z_j[5]+str(1)),2)]),2)])*(theta_i[int(str(1))])
	
	
	
	P_G['index']=[Gzz,Gzo,Goz,Goo]













""" CPT table complete"""


# E M Algorithm Starts

"""
E- Step

over here u initialized the CPT of the unobserved variables (eitehr by intution or either randomly) in theta_i at time t=0

Then you will calculate the Probability for each unobserved variable X_i for each Z_j i.e for each dataset
of the whole Data Set Z using marginal distribution over Theta_i at time t

"""


for i in range(len(Z)):
	
	
	y = Z[i].split(" ");
	
	#D I P G L S W A    1,3,4,5,6 can be x
	# as u know the values of atleast 6 variables you don't need variable elimination	
	# its like saying calculate probability of I=0/1 given rest values
	"""
	if y[3]== 'x':
		
		P_I.append([0.0,0.0])
		print(P_I[i])
	elif y[4]== 'x':
		P_I.append([0.0,0.0])
		print(P_I[i])
	elif y[5]== 'x':
		P_I.append([0.0,0.0])
		print(P_I[i])
	elif y[6]== 'x':
		P_I.append([0.0,0.0])
		print(P_I[i])
	else:
		OnlyImissing(Z[i].split(" "),i)
	"""