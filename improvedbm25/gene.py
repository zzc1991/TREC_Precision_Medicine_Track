import pymongo
import numpy as np
import re

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords3"] #pubmed中所有的词频、化学词、关键词和主题词表

wordfreq_list=[]
ChemicalName_List=[]
MeshHeadingName_List=[]
Keywords_List =[]
i=0
cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
for x in mywords.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                      no_cursor_timeout=True):
    wordfreq = x['wordfreq']
    ChemicalNameList = x['ChemicalNameList']
    MeshHeadingNameList = x['MeshHeadingNameList']
    KeywordsList = x['KeywordsList']
    if len(wordfreq)>0:
        wordfreq_list.append(len(wordfreq))
    if len(ChemicalNameList)>0:
        ChemicalName_List.append(len(ChemicalNameList))
    if len(MeshHeadingNameList)>0:
        MeshHeadingName_List.append(len(MeshHeadingNameList))
    if len(KeywordsList)>0:
        Keywords_List.append(len(KeywordsList))
    i=i+1
    print(i)

print('wordfreq:'+str(np.mean(wordfreq_list)))
print('ChemicalName_List:'+str(np.mean(ChemicalName_List)))
print('MeshHeadingName_List:'+str(np.mean(MeshHeadingName_List)))
print('Keywords_List:'+str(np.mean(Keywords_List)))
