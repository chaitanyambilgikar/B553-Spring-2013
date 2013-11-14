""" ********************* Homework 2 *********************
Done By: Chaitanya Bilgikar (cbilgika) & Nasheed Moiz (nmoiz)
"""	

from __future__ import division
from collections import defaultdict
from itertools import *
from pprint import pprint
from optparse import OptionParser
import sys
import random
import os
import operator
def constant_factory(value):
	return repeat(value).next

#list of tags-global
taglist = ['ADJ','ADV','ADP','CONJ','DET','NOUN','NUM','PRON','PRT','VERB','X','.']

#get list of words from the file
def get_word_list(filename):
	main_list=list()
	f= open(filename,"r").read()
	sentence_list=f.split('. .\n\n')
	sent_number=1
	for sent in sentence_list:
		if sent:
			word_count=1
			temp_dict = {}
			for pair in sent.split('\n'):
				temp=pair.split()
				if temp:
					temp_dict[word_count]=(temp[0].lower(),temp[1])
					word_count+=1
			main_list.append(temp_dict)
	return main_list

"""********************************** STEP 1 - LEARNING ***************************** """
tag_prob={}
word_given_tag = defaultdict(constant_factory(0.00000001))
word_prob =defaultdict(constant_factory(0.00000001))
tag_priors= defaultdict(int)
s1_prob={}


def learning():
	#main list of senctences
	main_list = list()
	#filename to train
	filename='./data/bc.train'
	main_list=get_word_list(filename)
	"""start by initializing all the dictionaries
	tag_prob is the probability of all the tags (occuring at any position)
	s1_prob is the probability of each tag being the first tag in a sentence
	tag_priors is the probability of a tag given the previous tag
	word_given_tag is the probability of a word given the tag"""
	total_words=0
	
	tag_count={} #dictionary to keep count of no. of occurances of each tag-used in word_given_tag
	for i in taglist:
		tag_prob[i]=0
		tag_count[i]=0
		s1_prob[i]=0.0
  		#initialize every possible adjacent combination of tags
		for sv in taglist:
			tag_priors[(i,sv)]=0.0 

	#adding values to all the dictionaries defined above
	sentence_count=0
	for sent in main_list:
		s1_prob[sent[1][1]]+=1
		sentence_count+=1
		#defining s(i+1)|s(i), dictionary having keys as a tuple as (next tag,current tag)
	    #count up occurences of a word occuring as a tag, and occurences of the word as any tag
		for k,v in sent.iteritems():
			tag_prob[v[1]]+=1
			tag_count[v[1]]+=1
			total_words+=1
			if k+1 in sent.keys():
				tag_priors[(sent[k+1][1],sent[k][1])]+=1
			word_given_tag[(v[0],v[1])]+=1
			word_prob[v[0]]+=1



	for k in tag_prob.keys():
		tag_prob[k]/=total_words
	#defining the total number of words 
	#converting s1_prob to probability values
	for tag in s1_prob:
		s1_prob[tag]/=sentence_count
	#coverting the tag priors to probabilities
	for (tag1,tag2) in tag_priors.keys():
		if tag_count[tag2]!=0:
			tag_priors[(tag1,tag2)]/=tag_count[tag2]
	#coverting the count of each word into its probabilitites
	for eachword in word_prob.keys():
		word_prob[eachword]/=total_words
	#getting prob. of the words given the tag
	for key,value in word_given_tag.iteritems():
		if tag_count[key[1]]!=0:
			word_given_tag[key]/=tag_count[key[1]]

""" ***************************** STEP 2 - NAIVE BAYES ******************************"""
#helper function needed in naive bayes
def makedict():
	d={}
	for tag in taglist:
		d[tag]=0.0
	return d
def naive_bayes():
	file_to_test = './data/bc.test'
	test_data = get_word_list(file_to_test)
	correct_sent = 0 #keeping a count of the number of correctly predicted senteces
	sentence_count=0 #keep a count of number of sentences
	prob_tag_given_word = {}
	for sent in test_data:
		sentence_count+=1
		correct = 0
		n = len(sent.keys()) #n is the number of words in the sentence
		each_sent_dict={}
		original_word_list=list()
		original_tag_list = list()

		#predicting the tag for every word in the sentence
		for k,word in sent.iteritems():
			original_word_list.append(word[0])
			original_tag_list.append(word[1])
			temp_dict = {}
			temp_dict = makedict()
			if word_prob[word[0]]!=0:
				for eachtag in temp_dict.keys():
					temp_dict[eachtag] = (word_given_tag[(word[0],eachtag)]*tag_prob[eachtag])/word_prob[word[0]]
			else:
				for eachtag in temp_dict.keys():
					temp_dict[eachtag] = 0.000000001
			each_sent_dict[word[0]]=max(temp_dict.iteritems(),key=operator.itemgetter(1))[0]
			
			if word[1]==each_sent_dict[word[0]]:
				correct+=1

		#print the sentence
		print "\n******ORIGINAL SENTENCE*********** "
		for i in original_word_list:
			print i,

		#print the original tags for this sentence
		print "\nGROUND TRUTH: "
		for i in original_tag_list:
			print i,

		print "\n PREDICTED TAGS: "
		#print the tags we predict
		for i in original_word_list:
			print each_sent_dict[i],

		correct_word_per = (correct/n)*100
		if correct_word_per==100.0:
			correct_sent+=1
		print "\nPercentage of correct tags in sentence: ",correct_word_per
		print "\n*******************************************************************"
		prob_tag_given_word[k]=each_sent_dict
	print "\nPercentage of sentences predicted correctly in Naive Bayes: ",(correct_sent/sentence_count)*100
		
""" ****************************** STEP 3 & 4 - MAX MARGINAL AND SAMPLING **********************"""
#function to pick a tag based on the weight of the tag (weight is the probability here)
def sample(sampleProbability):
	for k,v in sampleProbability.iteritems():
		sampleProbability[k]=int((1/(1-sampleProbability[k]))%306259)
	a= random.choice([x for x in sampleProbability for y in range(sampleProbability[x])]);
	#print a
	return a;
#function to do sampling. Called from max_marginal()
def sampling(tau,sent):
	
	sample_from = {}
	sample_from[1]=tau[1]
	sampled = {}
	for index,word in sent.iteritems():
		temp_dict={}
		if index==1:
			fixed_s1 = (max(tau[1].iteritems(),key=operator.itemgetter(1))[0],max(tau[1].iteritems(),key=operator.itemgetter(1))[1])
			#fixed_s1 = sample(sample_from[index])
			#print "fixed: ",fixed_s1
			sampled[index] = fixed_s1[0]
		else:
			for tag in taglist:
				temp_dict[tag] = (word_given_tag[(word[0],tag)]*tag_priors[(tag,sampled[index-1])])/((word_prob[word[0]])*tag_prob[sampled[index-1]])

			sample_from[index] = temp_dict
			sampled[index] = sample(sample_from[index])
			#print "sampled[",index,"] = ",sampled[index]
	return sampled


#function that performs forward-backward for the network	
def max_marginal():

	file_to_test = './data/bc.test'
	test_data = get_word_list(file_to_test)
	correct_sent = 0
	sentence_count=0
	for sent in test_data:
		""" tau will be a dictionary of dictionaries. The first key is the index of the word.
		The second key (or the keys in the dictionary will be all the tags). We can get the 
		most probable tag by doing an arg max over each of the dictionaries within tau.
		tau will be defined recursively.
		in general, tau(Si) = Sum over all values of Si-1 (P(si|si-1)P(wi|si)tau(Si-1))
		This is for all tau from 2 to N
		For tau(S1) = P(S1|W1) = P(W1|S1)P(S1)
		Note that there is no normalization constant anywhere.
		"""
		sentence_count+=1
		tau = {} #main tau
		tau_f = {} #forward tau
		tau_r= {} #backward tau
		n = len(sent.keys()) #n is the number of words in the sentence
		original_word_list=list()
		original_tag_list = list()
		"""alpha = 0.0
		for index,word in sent.iteritems():
			alpha += word_prob[ word[0] ]
		alpha= 1/alpha"""


		#forward part
		for i in range(1,n+1): #going through all the words
			current_word = sent[i][0]
			original_word_list.append(sent[i][0])
			original_tag_list.append(sent[i][1])
			if i == 1:
				temp_tau = {}
				for tag in taglist:
					temp_tau[tag] = word_given_tag[(current_word,tag)]*s1_prob[tag]					
				tau_f[1] = temp_tau
			else:
				temp_tau = {}
				for tag in taglist:
					temp_sum=0.0
					for tag2 in taglist:
						temp_sum+=(tag_priors[(tag,tag2)]*tau_f[i-1][tag2]*word_given_tag[(current_word,tag)])

					temp_tau[tag] = temp_sum
				tau_f[i] = temp_tau
		#backward part
		for i in range(n,0,-1):
			current_word = sent[i][0]

			if i == n:
				temp_tau={}

				for tag in taglist:
					temp_tau[tag]=1.0
				tau_r[n]=temp_tau
			else:
				temp_tau={}
				next_word = sent[i+1][0]
				for tag in taglist:
					temp_sum=0.0
					for tag2 in taglist:
						temp_sum+= word_given_tag[(next_word,tag2)]*tag_priors[(tag2,tag)]*tau_r[i+1][tag2]

					temp_tau[tag]=temp_sum
				tau_r[i] = temp_tau
		for i in tau_f.keys():
			tau[i]={}
			for tag in tau_f[i].keys():
				tau[i][tag] = tau_f[i][tag]*tau_r[i][tag]
		#calculating the accuracy of forward backward
		correct = 0
		calculated={}
		for key,word in sent.iteritems():
			calculated[key] = max(tau[key].iteritems(),key=operator.itemgetter(1))[0]
			#print word[1],calculated
			if word[1] == calculated[key]:
				correct+=1
		correct_word_per = (correct/n)*100
		if correct_word_per == 100.0:
			correct_sent+=1
		#print the original sentence
		print "\n******ORIGINAL SENTENCE*********** "
		for i in original_word_list:
			print i,
		#print the original tags for this sentence
		print "\nGROUND TRUTH: "
		for i in original_tag_list:
			print i,
		#print the predicted tags
		print "\nPREDICTED TAGS: "
		for k in calculated.keys():
			print calculated[k],
		print "\nPercentage correct tags for this sentence = ",correct_word_per
		print "\n****************************************************************************"
		#call to sampling, done 5 times
		sent_corr = 0 # flag to check correctly predicted sentences in sampling
		for i in range(1,6):
			
			print "\nSample ",i,"\n"
			sampled = sampling(tau,sent)
			correct=0
			for k in sampled.keys():
				print sampled[k],
			for k,v in sent.iteritems():
				if sent[k][1]==sampled[k]:
					correct+=1
			print "\n Percentage Correct (sampling): ",(correct/n)*100
			#set if at least one sentence is predicted correctly
			if ((correct/n)*100) == 100.0:
				sent_corr = 1
		if sent_corr == 1:
			print "\n Percentage of correct sentences in sampling: 100.0 (Found atleast 1 correct match)"  

	print "\nPercentage of sentences Predicted correctly in Max Marignal: ",(correct_sent/sentence_count)*100
	print "\n************************************************************************************************"
		
if __name__ == '__main__':
	parser = OptionParser("usage: %prog [options] algorithm \n"\
			"eg \t %prog -c m \n"\
			"\t %prog -c n \n")

	parser.add_option("--c", "--C", "-C","-c",dest="algorithm",help="Specify which of 2 algorithms to evaluate testing data on \n")




	(options, args) = parser.parse_args()
	algorithm = options.algorithm

	learning()
	if options.algorithm == None:
		print "=================NAIVE BAYES================"
		naive_bayes()
	elif options.algorithm.startswith("m"):
		print "\n=================MAX MARGINAL================="
		max_marginal()
	elif options.algorithm.startswith("n"):
		print "=================NAIVE BAYES================"
		naive_bayes()
	
