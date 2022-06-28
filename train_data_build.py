import csv
from collections import Counter
import re

import pymongo
from lxml import etree
import pandas as pd
import string
from sklearn.model_selection import train_test_split

import nltk
import torch
from transformers import BertTokenizer, BertModel, BertForSequenceClassification

from nltk.tokenize import sent_tokenize
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import torch.nn.functional as F
from args import Config
from train_eval import test
from utils import get_input, get_examples, get_data_loader, get_examples_PMID, get_data_loader_PMID

def build_data(f_pm, f_qrels, mypapers, f_output):
    f1 = open(f_output, 'a', encoding='utf-8', newline='')
    # 2. 基于文件对象构建 csv写入对象
    csv_writer = csv.writer(f1)

    tree = etree.parse(f_qrels)
    topics_list = tree.xpath('//topic')

    f_pm_2017 = open(f_pm)
    topics_dict = {}
    for item in f_pm_2017.readlines():
        item_list = item.replace('\n', '').split(' ')
        if str(item_list[2]).isdigit():
            if item_list[0] in topics_dict.keys():
                temp_dict = topics_dict[item_list[0]]
                temp_dict[item_list[2]] = item_list[3]
                topics_dict[item_list[0]] = temp_dict
            else:
                temp_dict = {item_list[2]: item_list[3]}
                topics_dict[item_list[0]] = temp_dict

    i = 0
    for topics in topics_list:
        disease = topics.xpath('disease/text()')
        gene = topics.xpath('gene/text()')
        demographic = topics.xpath('demographic/text()')
        i = i + 1
        topics_ref = list(topics_dict[str(i)].keys())
        char = disease[0] + ' ' + gene[0] + ' ' + demographic[0]
        print(topics_ref)
        print(char)
        for paper in mypapers.find({'PMID': {'$in': topics_ref}}):
            if len(paper['Abstract']) > 0:
                if int(topics_dict[str(i)][str(paper['PMID'])]) == 2 or int(
                        topics_dict[str(i)][str(paper['PMID'])]) == 1:
                    csv_writer.writerow([paper['Abstract'], char, 1])
                else:
                    csv_writer.writerow([paper['Abstract'], char, 0])


def get_pmids(file_list, mypapers, mypapers_ref):
    PMID_list = []
    for f_pm in file_list:
        f = open(f_pm)
        for item in f.readlines():
            print(item)
            item_list = item.replace('\n', '').split(' ')
            PMID_list.append(item_list[2])
    result = list(Counter(PMID_list).keys())
    for paper in mypapers.find({'PMID': {'$in': result}}):
        ss = mypapers_ref.insert_one(paper)
        print(ss)


def read_data(data_file, train_size, test_size, random_state):
    dataset = pd.read_csv(data_file)
    dataset = dataset.dropna(inplace=False)
    train_x, test_x = train_test_split(dataset, train_size=train_size, test_size=test_size, random_state=random_state,
                                       stratify=dataset.loc[:, ['label']])
    return train_x, test_x


def get_age(text):
    reg = r'(?:9\d|8\d|7\d|6\d|5\d|4\d|3\d|2\d|1\d|1\d{2})(?:\s|-)?years?(?:\s|-)?old'
    age = re.findall(reg, text)
    return age


def get_disease_gene(abstract, disease_gene, stopwords, porter):
    abstract_tokens = nltk.word_tokenize(abstract.lower())
    abstract_clean = [clean_word(w) for w in abstract_tokens if w.lower() not in stopwords]
    abstract_stem = [porter.stem(t) for t in abstract_clean]

    disease_gene_tokens = nltk.word_tokenize(disease_gene.lower())
    disease_gene_clean = [clean_word(w) for w in disease_gene_tokens if w.lower() not in stopwords]
    disease_gene_stem = [porter.stem(t) for t in disease_gene_clean]

    result = list(set(abstract_stem).intersection(set(disease_gene_stem)))
    return result


def clean_word(word):
    return re.sub('[!"#$%&\'()*+,-./:;<=>?@，。?★、…【】《》？“”‘！[\\]^_`{|}~\s]+', "", word.lower())


def get_disease_gene_dict(disease_gene, stopwords, porter):
    disease_gene_dict = {}
    disease_gene_tokens = nltk.word_tokenize(disease_gene)
    disease_gene_tokens_lower = nltk.word_tokenize(disease_gene)
    disease_gene_clean = [clean_word(w) for w in disease_gene_tokens_lower if w.lower() not in stopwords]
    disease_gene_stem = [porter.stem(t) for t in disease_gene_clean]

    for count, dg in enumerate(disease_gene_tokens):
        disease_gene_dict[disease_gene_stem[count]] = dg.lower()
    return disease_gene_dict


def clean_text(text, stopwords):
    abstract_tokens = nltk.word_tokenize(text)
    abstract_clean = [clean_word(w) for w in abstract_tokens if w.lower() not in stopwords]
    return ' '.join(filter(str.isalnum, abstract_clean))


def clean_text_list(text, stopwords):
    abstract_tokens = nltk.word_tokenize(text)
    abstract_clean = [clean_word(w) for w in abstract_tokens if w.lower() not in stopwords]
    return list(filter(str.isalnum, abstract_clean))


def check_words(abstract, word, porter):
    result = False
    abstract_tokens = nltk.word_tokenize(abstract.lower())
    abstract_stem = [porter.stem(t) for t in abstract_tokens]

    word_stem = [porter.stem(str(word).lower())]
    if len((set(abstract_stem).intersection(set(word_stem)))) > 0:
        result = True

    return result


def get_age_class(demographic):
    result = ''
    isadult = False
    age = int("".join(list(filter(str.isdigit, demographic))))
    if age < 2:
        result = 'Infant'
    elif 2 <= age < 6:
        result = 'Preschool'
    elif 6 <= age < 13:
        result = 'Child'
    elif 13 <= age < 19:
        result = 'Adolescent'
    elif 19 <= age < 35:
        result = 'Young'
    elif 35 <= age < 60:
        result = 'Middle Aged'
    elif 60 <= age < 80:
        result = 'Aged'
    elif age >= 80:
        result = 'Aged, 80 and over'

    if age >= 18:
        isadult = True
    return result, isadult


def get_clean_abstract(text, stopwords):
    abstract_tokens = nltk.word_tokenize(text.lower())
    abstract_clean = [clean_word(w) for w in abstract_tokens if w.lower() not in stopwords]

    return ' '.join(abstract_clean), len(abstract_clean)


def process_data(f_pm, f_qrels, mypapers_ref, stopwords, porter, pad_size, f_output, tokenizer, model, pmid_dict,
                 query_dict,len_topics):
    '''
    Generate a training set with negative samples
    :param f_pm:topic_files
    :param f_qrels:qrels_files
    :param mypapers_ref:
    :param stopwords:
    :param porter:
    :param pad_size:bert_pad_size
    :param f_output:result_file
    :param tokenizer:bert_tokenizer
    :param model:
    :param pmid_dict:
    :param query_dict:
    :return:
    '''
    f1 = open(f_output, 'a', encoding='utf-8', newline='')
    # 2. 基于文件对象构建 csv写入对象
    csv_writer = csv.writer(f1)

    tree = etree.parse(f_qrels)
    topics_list = tree.xpath('//topic')

    f_p = open(f_pm)
    topics_dict = {}
    for item in f_p.readlines():
        item_list = item.replace('\n', '').split(' ')
        if str(item_list[2]).isdigit():
            if item_list[0] in topics_dict.keys():
                temp_dict = topics_dict[item_list[0]]
                temp_dict[item_list[2]] = item_list[3]
                topics_dict[item_list[0]] = temp_dict
            else:
                temp_dict = {item_list[2]: item_list[3]}
                topics_dict[item_list[0]] = temp_dict

    for count, topics in enumerate(topics_list):
        print(count)
        disease = topics.xpath('disease/text()')
        gene = topics.xpath('gene/text()')
        demographic = topics.xpath('demographic/text()')
        topics_ref = list(topics_dict[str(count + 1)].keys())
        dis_gene = disease[0] + ' ' + gene[0]
        age = str(demographic[0]).split(' ')[0]
        sex = str(demographic[0]).split(' ')[1]
        for paper in mypapers_ref.find({'PMID': {'$in': topics_ref}}):
            if len(paper['Abstract']) > 0:
                elements = ''
                dis_gene = ' '.join(clean_text_list(dis_gene, stopwords)).lower()
                disease_gene_dict = get_disease_gene_dict(dis_gene, stopwords, porter)
                result = get_disease_gene(str(paper['Abstract']), dis_gene, stopwords, porter)
                if len(result) > 0:
                    for item in disease_gene_dict:
                        if item in result:
                            elements = elements + disease_gene_dict[item] + ' '

                Mesh_list = []
                for Mesh in paper['MeshHeadingNameList']:
                    Mesh_list.append(Mesh['MeshHeadingName'])
                if 'Humans' in Mesh_list:
                    elements = elements + 'human' + ' '
                else:
                    if check_words(str(paper['Abstract']), 'humans', porter):
                        elements = elements + 'human' + ' '
                result, isadult = get_age_class(age)
                if result in Mesh_list:
                    elements = elements + result.lower() + ' '
                else:
                    age_list = get_age(str(paper['Abstract']))
                    if len(age_list) > 0:
                        result, isadult = get_age_class(age_list[0])
                        elements = elements + result.lower() + ' '
                if isadult and 'Adult' in Mesh_list:
                    elements = elements + 'adult' + ' '
                sex = string.capwords(sex)
                if sex in Mesh_list:
                    elements = elements + sex.lower() + ' '
                else:
                    if check_words(str(paper['Abstract']), sex, porter):
                        elements = elements + sex.lower() + ' '
                    else:
                        if check_words(str(paper['Abstract']), 'woman', porter):
                            elements = elements + 'female' + ' '
                        elif check_words(str(paper['Abstract']), 'man', porter):
                            elements = elements + 'male' + ' '
                        elif check_words(str(paper['Abstract']), 'boy', porter):
                            elements = elements + 'male' + ' '
                        elif check_words(str(paper['Abstract']), 'girl', porter):
                            elements = elements + 'female' + ' '
                elements = clean_text(elements, stopwords)
                ab_clean, ab_length = get_clean_abstract(str(paper['Abstract']), stopwords)

                if ab_length + len(elements) <= pad_size:
                    if len(elements) > 0:
                        abst = elements + '.' + ab_clean
                    else:
                        abst = ab_clean
                else:
                    abstract_summary = summary(str(paper['Abstract']), stopwords, len(elements), tokenizer, model)
                    if len(elements) > 0:
                        abst = elements + '.' + abstract_summary
                    else:
                        abst = abstract_summary

                Mesh_age, isadult = get_age_class(age)

                if int(topics_dict[str(count + 1)][str(paper['PMID'])]) == 2 or int(
                        topics_dict[str(count + 1)][str(paper['PMID'])]) == 1:
                    csv_writer.writerow(
                        [str(paper['PMID']), abst, str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower(),
                         1])
                    for item in pmid_dict:
                        if int(item) != int(count + len_topics+1):
                            if str(paper['PMID']) not in pmid_dict[item]:
                                csv_writer.writerow(
                                    [str(paper['PMID']), abst, query_dict[item], 0])
                else:
                    csv_writer.writerow(
                        [str(paper['PMID']), abst, str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower(),
                         0])



def summary(text, stopwords, number_elements, tokenizer, model):
    sentences = []
    # Tokenize sentences
    sentences_raw = sent_tokenize(text)
    text_length = 0
    for sentence in sentences_raw:
        tokens = nltk.word_tokenize(sentence)
        abstract_clean = [w.lower() for w in tokens if w.lower() not in stopwords]
        sentences.append(' '.join(abstract_clean))
        text_length = text_length + len(abstract_clean)
    results = []
    for sentence in sentences:
        if len(sentence) > 512:
            sentence = sentence[0:511]
            input_id = torch.tensor(tokenizer.encode(sentence)).unsqueeze(0)
            outputs = model(input_id)
            pooled_output = outputs[1]  # 句向量
            results.append(pooled_output.tolist()[0])
        else:
            input_id = torch.tensor(tokenizer.encode(sentence)).unsqueeze(0)
            outputs = model(input_id)
            pooled_output = outputs[1]  # 句向量
            results.append(pooled_output.tolist()[0])
    # Average the word emebddings in each sentence to create a single array for a sentence
    averaged = []
    for sent in results:
        averaged.append(np.mean(sent, axis=0, dtype=np.float64))

    # # Cluster the data
    n_clusters = int(
        np.ceil(len(averaged) * (350 / (
                text_length - number_elements))))  # Our output summary will be 30% of the size of the input corpus
    if n_clusters>len(averaged):
        n_clusters=len(averaged)
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans = kmeans.fit(np.array(np.array(averaged).reshape(-1, 1)))
    #
    # # Summarization
    avg = []
    for j in range(n_clusters):
        idx = np.where(kmeans.labels_ == j)[0]
        avg.append(np.mean(idx))
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, np.array(averaged).reshape(-1, 1))
    ordering = sorted(range(n_clusters), key=lambda k: avg[k])
    summary = ' '.join([sentences[closest[idx]] for idx in ordering])

    # Display the summary
    return summary


def insert_top_1000(mycount, top_1000):
    mydoc = mycount.find().sort("bm25_score", -1).limit(1000)
    for x in mydoc:
        mydict = {"PMID": x['PMID'], "related": x['related'], "bm25_score": x['bm25_score']}
        top_1000.insert_one(mydict)


def get_top_1000_PMIDS(top_1000):
    result = {}
    mydoc = top_1000.find()
    for x in mydoc:
        result[x['PMID']] = x['related']
    return result


def cal_com_score(top_1000, top_1000_s, result_1000):
    result_relate = {}
    result_bm25 = {}
    result_bert = {}
    mydoc = top_1000.find()
    for x in mydoc:
        result_relate[x['PMID']] = x['related']
        result_bm25[x['PMID']] = x['bm25_score']
    mydoc_s = top_1000_s.find()
    for x in mydoc_s:
        result_bert[x['PMID']] = x['score']
    for item in result_bm25:
        mydict = {"PMID": item, "related": result_relate[item],
                  "com_score": float(result_bm25[item]) * float(result_bert[item])}
        result_1000.insert_one(mydict)


def data_builder(disease, gene, demographic, mypapers_ref, topics_ref, pad_size, tokenizer, model):
    data_sent = []
    dis_gene = disease + ' ' + gene
    age = str(demographic).split(' ')[0]
    sex = str(demographic).split(' ')[1]
    print(len(topics_ref))
    for paper in mypapers_ref.find({'PMID': {'$in': topics_ref}}):
        if len(paper['Abstract']) > 0:
            elements = ''
            dis_gene = ' '.join(clean_text_list(dis_gene, stopwords)).lower()
            disease_gene_dict = get_disease_gene_dict(dis_gene, stopwords, porter)
            result = get_disease_gene(str(paper['Abstract']), dis_gene, stopwords, porter)
            if len(result) > 0:
                for item in disease_gene_dict:
                    if item in result:
                        elements = elements + disease_gene_dict[item] + ' '

            Mesh_list = []
            for Mesh in paper['MeshHeadingNameList']:
                Mesh_list.append(Mesh['MeshHeadingName'])
            if 'Humans' in Mesh_list:
                elements = elements + 'human' + ' '
            else:
                if check_words(str(paper['Abstract']), 'humans', porter):
                    elements = elements + 'human' + ' '
            result, isadult = get_age_class(age)
            if result in Mesh_list:
                elements = elements + result.lower() + ' '
            else:
                age_list = get_age(str(paper['Abstract']))
                if len(age_list) > 0:
                    result, isadult = get_age_class(age_list[0])
                    elements = elements + result.lower() + ' '
            if isadult and 'Adult' in Mesh_list:
                elements = elements + 'adult' + ' '
            sex = string.capwords(sex)
            if sex in Mesh_list:
                elements = elements + sex.lower() + ' '
            else:
                if check_words(str(paper['Abstract']), sex, porter):
                    elements = elements + sex.lower() + ' '
                else:
                    if check_words(str(paper['Abstract']), 'woman', porter):
                        elements = elements + 'female' + ' '
                    elif check_words(str(paper['Abstract']), 'man', porter):
                        elements = elements + 'male' + ' '
                    elif check_words(str(paper['Abstract']), 'boy', porter):
                        elements = elements + 'male' + ' '
                    elif check_words(str(paper['Abstract']), 'girl', porter):
                        elements = elements + 'female' + ' '
            elements = clean_text(elements, stopwords)
            ab_clean, ab_length = get_clean_abstract(str(paper['Abstract']), stopwords)

            if ab_length + len(elements) <= pad_size:
                if len(elements) > 0:
                    abst = elements + '.' + ab_clean
                else:
                    abst = ab_clean
            else:
                abstract_summary = summary(str(paper['Abstract']), stopwords, len(elements), tokenizer, model)
                if len(elements) > 0:
                    abst = elements + '.' + abstract_summary
                else:
                    abst = abstract_summary

            Mesh_age, isadult = get_age_class(age)
            sentence1, sentence2 = abst, str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower()
            data_sent.append([[sentence1, sentence2], 0, paper['PMID']])
        else:
            Mesh_age, isadult = get_age_class(age)
            data_sent.append(
                [['None', str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower()], 0, paper['PMID']])
    print(len(data_sent))
    return data_sent


def rerank(top_1000, top_1000_s, dis_gene, gene, demographic, mypapers, pad_size, tokenizer, model, config, sim_model):
    PMIDS = get_top_1000_PMIDS(top_1000)
    print(PMIDS)
    data = data_builder(dis_gene, gene, demographic, mypapers, list(PMIDS.keys()), pad_size, tokenizer,
                        model)
    data_summary = get_data_loader_PMID(config, data)
    result = test(config, sim_model, data_summary)
    for item in result:
        mydict = {"PMID": str(item), "related": PMIDS[str(item)], "score": result[item]}
        top_1000_s.insert_one(mydict)


def get_PMIDs_query_dict(file_list, topic_file, stopwords):
    pmid_list = []
    for file in file_list:
        f_p = open(file)
        topic = 1
        topics_dict = {}
        for item in f_p.readlines():
            item_list = item.replace('\n', '').split(' ')
            if topic != int(item_list[0]):
                topic += 1
            if str(item_list[2]).isdigit() and str(item_list[3]) != '0':
                if topic in topics_dict.keys():
                    temp_list = topics_dict[topic]
                    temp_list.append(item_list[2])
                    topics_dict[topic] = temp_list
                else:
                    temp_list = [item_list[2]]
                    topics_dict[topic] = temp_list
        pmid_list.append(topics_dict)
    t = 1
    pmid_dict, query_dict, disease_dict, gene_dict = {}, {}, {}, {}
    for item in pmid_list:
        for _topic in item:
            pmid_dict[t] = item[_topic]
            t += 1
    q = 1
    for f_qrels in topic_file:
        tree = etree.parse(f_qrels)
        topics_list = tree.xpath('//topic')
        for count, topics in enumerate(topics_list):
            disease = topics.xpath('disease/text()')
            gene = topics.xpath('gene/text()')
            demographic = topics.xpath('demographic/text()')
            dis_gene = disease[0] + ' ' + gene[0]
            age = str(demographic[0]).split(' ')[0]
            sex = str(demographic[0]).split(' ')[1]
            dis_gene = ' '.join(clean_text_list(dis_gene, stopwords)).lower()
            Mesh_age, isadult = get_age_class(age)
            query_dict[q] = str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower()
            if disease[0] not in disease_dict:
                disease_dict[disease[0]] = [q]
            else:
                temp_list = disease_dict[disease[0]]
                temp_list.append(q)
                disease_dict[disease[0]] = temp_list
            gene_list = str(gene[0]).split(', ')
            print(gene_list)
            for gene in gene_list:
                if gene not in gene_dict:
                    gene_dict[gene] = [q]
                else:
                    temp_list = gene_dict[gene]
                    temp_list.append(q)
                    gene_dict[gene] = temp_list
            q += 1

    return pmid_dict, query_dict, disease_dict, gene_dict


def process_data_classify(f_pm, f_qrels, mypapers_ref, stopwords, porter, pad_size, f_output, tokenizer, model,
                          pmid_dict,
                          query_dict, disease_dict, gene_dict):
    f1 = open(f_output, 'a', encoding='utf-8', newline='')
    # 2. 基于文件对象构建 csv写入对象
    csv_writer = csv.writer(f1)

    tree = etree.parse(f_qrels)
    topics_list = tree.xpath('//topic')

    f_p = open(f_pm)
    topics_dict = {}
    for item in f_p.readlines():
        item_list = item.replace('\n', '').split(' ')
        if str(item_list[2]).isdigit():
            if item_list[0] in topics_dict.keys():
                temp_dict = topics_dict[item_list[0]]
                temp_dict[item_list[2]] = item_list[3]
                topics_dict[item_list[0]] = temp_dict
            else:
                temp_dict = {item_list[2]: item_list[3]}
                topics_dict[item_list[0]] = temp_dict

    for count, topics in enumerate(topics_list):
        print(count)
        disease = topics.xpath('disease/text()')
        gene = topics.xpath('gene/text()')
        gene_list = str(gene[0]).split(', ')
        demographic = topics.xpath('demographic/text()')
        topics_ref = list(topics_dict[str(count + 1)].keys())
        dis_gene = disease[0] + ' ' + gene[0]
        age = str(demographic[0]).split(' ')[0]
        sex = str(demographic[0]).split(' ')[1]
        for paper in mypapers_ref.find({'PMID': {'$in': topics_ref}}):
            if len(paper['Abstract']) > 0:
                elements = ''
                dis_gene = ' '.join(clean_text_list(dis_gene, stopwords)).lower()
                disease_gene_dict = get_disease_gene_dict(dis_gene, stopwords, porter)
                result = get_disease_gene(str(paper['Abstract']), dis_gene, stopwords, porter)
                if len(result) > 0:
                    for item in disease_gene_dict:
                        if item in result:
                            elements = elements + disease_gene_dict[item] + ' '

                Mesh_list = []
                for Mesh in paper['MeshHeadingNameList']:
                    Mesh_list.append(Mesh['MeshHeadingName'])
                if 'Humans' in Mesh_list:
                    elements = elements + 'human' + ' '
                else:
                    if check_words(str(paper['Abstract']), 'humans', porter):
                        elements = elements + 'human' + ' '
                result, isadult = get_age_class(age)
                if result in Mesh_list:
                    elements = elements + result.lower() + ' '
                else:
                    age_list = get_age(str(paper['Abstract']))
                    if len(age_list) > 0:
                        result, isadult = get_age_class(age_list[0])
                        elements = elements + result.lower() + ' '
                if isadult and 'Adult' in Mesh_list:
                    elements = elements + 'adult' + ' '
                sex = string.capwords(sex)
                if sex in Mesh_list:
                    elements = elements + sex.lower() + ' '
                else:
                    if check_words(str(paper['Abstract']), sex, porter):
                        elements = elements + sex.lower() + ' '
                    else:
                        if check_words(str(paper['Abstract']), 'woman', porter):
                            elements = elements + 'female' + ' '
                        elif check_words(str(paper['Abstract']), 'man', porter):
                            elements = elements + 'male' + ' '
                        elif check_words(str(paper['Abstract']), 'boy', porter):
                            elements = elements + 'male' + ' '
                        elif check_words(str(paper['Abstract']), 'girl', porter):
                            elements = elements + 'female' + ' '
                elements = clean_text(elements, stopwords)
                ab_clean, ab_length = get_clean_abstract(str(paper['Abstract']), stopwords)

                if ab_length + len(elements) <= pad_size:
                    if len(elements) > 0:
                        abst = elements + '.' + ab_clean
                    else:
                        abst = ab_clean
                else:
                    abstract_summary = summary(str(paper['Abstract']), stopwords, len(elements), tokenizer, model)
                    if len(elements) > 0:
                        abst = elements + '.' + abstract_summary
                    else:
                        abst = abstract_summary

                Mesh_age, isadult = get_age_class(age)
                if int(topics_dict[str(count + 1)][str(paper['PMID'])]) == 2 or int(
                        topics_dict[str(count + 1)][str(paper['PMID'])]) == 1:
                    csv_writer.writerow(
                        [str(paper['PMID']), abst, str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower(),
                         1])
                    class_list = set()
                    for gene in gene_list:
                        class_list = class_list.union(set(gene_dict[gene]))
                    class_list = list(class_list.union(set(disease_dict[disease[0]]).union()))
                    print(class_list)
                    for same in class_list:
                        if int(same) != int(count + 31):
                            if str(paper['PMID']) not in pmid_dict[same]:
                                csv_writer.writerow(
                                    [str(paper['PMID']), abst, query_dict[same], 0])
                else:
                    csv_writer.writerow(
                        [str(paper['PMID']), abst, str(dis_gene).lower() + ' ' + Mesh_age.lower() + ' ' + sex.lower(),
                         0])


def process_rank_data(mypapers_ref, stopwords, pad_size, tokenizer, model,topic_papers,topk):
        topics_ref=[]
        pers= {}
        vector_list=[]
        related={}
        i=0
        for x in topic_papers.find().sort('score',-1).limit(topk):
            topics_ref.append(x['PMID'])
            pers[i]=float((x['score']))
            related[i]=x['related']
            i+=1

        for paper in mypapers_ref.find({'PMID': {'$in': topics_ref}}):
            abst=''
            if len(paper['Abstract']) > 0:
                if len(clean_text(str(paper['Abstract']),stopwords)) <= pad_size:
                        abst = clean_text(str(paper['Abstract']),stopwords)
                else:
                    abst = summary(str(paper['Abstract']), stopwords, 0, tokenizer, model)
                    abst_list=abst.split(' ')
                    if len(abst_list)>512:
                        abst_list=abst_list[0:511]
                        abst=' '.join(abst_list)
            input_id = torch.tensor(tokenizer.encode(abst)).unsqueeze(0)
            outputs = model(input_id)
            pooled_output = outputs[1]  # 句向量
            vector_list.append(pooled_output.tolist()[0])
        return vector_list,pers,related


if __name__ == "__main__":
    data_path = 'data'
    bert_path = 'biobert'
    myclient = pymongo.MongoClient("mongodb://localhost:27017/")
    mydb = myclient["pubmed"]
    mypapers_ref = mydb["papers"]
    stopwords = nltk.corpus.stopwords.words('english')
    porter = nltk.PorterStemmer()
    #
    config = Config(data_path, bert_path)
    sim_model = BertForSequenceClassification.from_pretrained(config.bert_path)

    tokenizer = BertTokenizer.from_pretrained(bert_path)  # 分词词
    model = BertModel.from_pretrained(bert_path)  # 模型

    mytopicdb = myclient["cs2019_bert_nomalize"]

