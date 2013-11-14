from __future__ import division
import sys
import random
import os
taglist = ['ADJ','ADV','ADP','CONJ','DET','NOUN','NUM','PRON','PRT','VERB','X','.']
f= open('./data/bc.tiny.test',"r").read()
sentence_list=f.split('. .\n\n')
main_list = list()
sent_number=1
for sent in sentence_list:
	if sent:
		word_count=1
		temp_dict = {}
		for pair in sent.split('\n'):
			temp=pair.split()
			if temp:
				temp_dict[word_count]=(temp[0],temp[1])
				word_count+=1
		main_list.append(temp_dict)
#defining tag_dist as a dict to get P(S1)
tag_dist = {}
totalsentcount=0
for tag in taglist:
	tag_dist[tag]=0
#defining tag_priors as a dict to store P(Si+1|Si)
tag_priors = {}
for tag1 in taglist:
	for tag2 in taglist:
		tag_priors[(tag1,tag2)]=0
#print len(tag_priors.keys())

for sent in main_list:
	tag_dist[sent[1][1]]+=1
	totalsentcount+=1
	for index in sent.keys():
		if index != 1:
			tag_priors[(sent[index-1][1],sent[index][1])]+=1



for tag in tag_dist:
	tag_dist[tag]/=totalsentcount




