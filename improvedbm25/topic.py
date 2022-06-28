# -*- coding: utf-8 -*-
import sys
import pymongo
import os

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mycol = mydb["topics2019"]

data = []
for line in open("2019.txt","r"): #设置文件对象并读取每一行文件
    data.append(line)               #将每一行文件加入到list中
    list=line.strip().split(' ')
    mydict = {"topic": str(list[0]),"PMID":str(list[2]),"relate":str(list[3])}
    x = mycol.insert_one(mydict)
    print(str(list[0]) + ',' +str(x))