#!/usr/bin/python
import os
import random
from args import Config
import time
from train_eval import train
from transformers import BertForSequenceClassification
from utils import get_sent, get_demo_sent, get_data_loader, printm, format_time,get_input


# data_path = os.path.join(os.getcwd(), 'github_version\\data')
data_path='data'
bert_path = 'bert'

if __name__ =='__main__':

	config = Config(data_path, bert_path)
	# memoryUtil = printm(config)
  	# print(memoryUtil)	
	# if memoryUtil > 0.2:
	# 	try:
	# 		os._exit(0)
	# 	except:
	# 		print('MemoryUtil FULL')

	start_time = time.time()
	print("Loading data...")
	#train_sent, dev_sent, test_sent = get_sent(data_path)
	train_sent, dev_sent, test_sent = get_demo_sent(data_path)
	print(train_sent)
	
	train_dataloader = get_data_loader(config, train_sent)
	dev_dataloader = get_data_loader(config, dev_sent)
	test_dataloader = get_data_loader(config, test_sent)
	print("Time usage:", format_time(time.time() - start_time))

	# train
	model = BertForSequenceClassification.from_pretrained(config.bert_path)
	train(config, model, train_dataloader, dev_dataloader, test_dataloader)
	input_ids, segment_ids, input_masks, label_ids =get_input(config,[[['A 2-year-old boy is brought to the emergency department by his parents for 5 days of high fever and irritability. The physical exam reveals conjunctivitis, strawberry tongue, inflammation of the hands and feet, desquamation of the skin of the fingers and toes, and cervical lymphadenopathy with the smallest node at 1.5 cm.','The abdominal exam demonstrates tenderness and enlarged liver. Laboratory tests report elevated alanine aminotransferase, white blood cell count of 17,580/mm, albumin 2.1 g/dL, C-reactive protein 4.5 mg, erythrocyte sedimentation rate 60 mm/h, mild normochromic, normocytic anemia, and leukocytes in urine of 20/mL with no bacteria identified'],0]])
	outputs = model(input_ids=input_ids, attention_mask=input_masks, token_type_ids=segment_ids, labels=label_ids)
	print(input_ids)
	print(segment_ids)
	print(input_masks)
