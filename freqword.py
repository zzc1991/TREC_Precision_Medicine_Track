# -*- coding: utf-8 -*-
import pymongo
import re
import nltk


myclient = pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords1"]
mypapers = mydb["papers"]

rstr = r"[\/\\\:\*\?\"\<\>\|.%$]"


def is_number(s):
    try:
        float(s)
        return False
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return False
    except (TypeError, ValueError):
        pass

    return True


def content_fraction(text):
    stopwords = nltk.corpus.stopwords.words('english')
    pattern = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9^-]")
    text1 = re.sub(pattern, ' ', text).lower()
    tokens = nltk.word_tokenize(text1)  # 抽词
    content = [w for w in tokens if w.lower() not in stopwords and is_number(w.lower()) and len(w) > 1]
    return content

s=0
for x in mypapers.find({},
                       {'PMID', 'ChemicalNameList', 'OtherAbstract', 'Abstract', 'KeywordsList', 'MeshHeadingNameList',
                        'ArticleTitle', 'JournalTitle'}, no_cursor_timeout=True):
    abstract = x['Abstract'] + ' ' + x['OtherAbstract'] + ' ' + x['ArticleTitle']
    words = content_fraction(abstract)
    dic = {}
    for word in words:
        if word not in dic:
            dic[word] = 1
        else:
            dic[word] += 1
    dic_freq = {k: v for k, v in dic.items() if v > 0}

    mydict = {"wordfreq": dic_freq, "abstract": abstract, "PMID": x['PMID'], "ChemicalNameList": x['ChemicalNameList'],
              "KeywordsList": x['KeywordsList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
              "JournalTitle": x['JournalTitle']}
    y = mywords.insert_one(mydict)
    print(str(s) + '---------' + str(y))
    s += 1
