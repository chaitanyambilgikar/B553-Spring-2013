from __future__ import division
import sys
import random
import os
taglist = ['ADJ','ADV','ADP','CONJ','DET','NOUN','NUM','PRON','PRT','VERB','X','.']
f= open('./data/bc.test',"r").read()
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
				temp_dict[word_count]=(temp[0].lower(),temp[1])
				word_count+=1
		main_list.append(temp_dict)
#defining the total number of words 
sentence_count=0
#defining tag_dist as a dict to get P(S1)
tag_dist={}
tag_dist_copy ={}
for tag in taglist:
	tag_dist[tag]=0.0
for tag in taglist:
	tag_dist_copy[tag]=0.0
#adding values to tag_dist

for sent in main_list:
	tag_dist[sent[1][1]]+=1
	tag_dist_copy[sent[1][1]]+=1
	sentence_count+=1

#converting tag_dist to probability values
for tag in tag_dist:
	tag_dist[tag]/=sentence_count

#defining s(i+1)|s(i), dictionary having keys as a tuple as (next tag,current tag)
tag_priors={}
for si in taglist:
	for sv in taglist:
		tag_priors[(si,sv)]=0.0
for sent in main_list:
	for i in sent.keys():
		if i+1 in sent.keys():
			next_tag = sent[i+1][1]
			curr_tag = sent[i][1]
			tag_priors[(next_tag,curr_tag)]+=1
for (tag1,tag2) in tag_priors.keys():
	if tag_dist_copy[tag2]!=0:
		tag_priors[(tag1,tag2)]/=tag_dist_copy[tag2]

#dictionary having keys as tuples with (word,tag)
word_given_tag = {}
for sent in main_list:
	for k,v in sent.iteritems():
		if (v[0],v[1]) in word_given_tag.keys():
			word_given_tag[(v[0],v[1])]+=1
		else:
			word_given_tag[(v[0],v[1])]=1
for key,value in word_given_tag.iteritems():
	if tag_dist_copy[key[1]]!=0:
		word_given_tag[key]/=tag_dist_copy[key[1]]
print tag_dist