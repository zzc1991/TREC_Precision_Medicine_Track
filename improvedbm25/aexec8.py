import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_8"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_8_related"]#聚类后对应与主题相关联的文献





def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_bladder = log((29138919 - 176012 + 0.5) / (176012 + 0.5), 10)
    idf_fgfr3 = log((29138919 - 1977 + 0.5) / (1977 + 0.5), 10)
    idf_s249c = log((29138919 - 37 + 0.5) / (37 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 37086 + 0.5) / (37086 + 0.5), 10)
    idf_ele_4 = log((13670358 - 8455 + 0.5) / (8455 + 0.5), 10)
    idf_ele_5 = log((13670358 - 22004 + 0.5) / (22004 + 0.5), 10)
    idf_ele_6 = log((13670358 - 22 + 0.5) / (22 + 0.5), 10)
    idf_ele_7 = log((13670358 - 24866 + 0.5) / (24866 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 51048 + 0.5) / (51048 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 9290 + 0.5) / (9290 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 1444 + 0.5) / (1444 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 22004 + 0.5) / (22004 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 24866 + 0.5) / (24866 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 23278 + 0.5) / (23278 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        len_freq=0
        bladder_score=0
        fgfr3_score=0
        s249c_score=0

        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        bladder = [True for x in wordfreq.items() if 'bladder' in x]


        # ---------------摘要统计-------------------#


        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'bladder' in key1:
                bladder_score = bladder_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'fgfr3' in key1:
                fgfr3_score = fgfr3_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 's249c' in key1:
                s249c_score = s249c_score + wordfreq[key]

        bm25_bladder_score = (((k1+1)*bladder_score)/((k1*(b1+(1-b1)*(len_freq/85)))+bladder_score))
        bm25_fgfr3_score =(((k1+1)*fgfr3_score)/((k1*(b1+(1-b1)*(len_freq/85)))+fgfr3_score))
        bm25_s249c_score = (((k1+1)*s249c_score)/((k1*(b1+(1-b1)*(len_freq/85)))+s249c_score))

        bm25_ab_score =idf_bladder*bm25_bladder_score+idf_fgfr3*bm25_fgfr3_score+idf_s249c*bm25_s249c_score

        idf_para=[{str(bladder_score):idf_bladder},{str(fgfr3_score):idf_fgfr3},{str(s249c_score):idf_s249c}]

        # ---------------共现分析摘要-------------------#
        if len(bladder)!=0 and bladder[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'fgfr3' in key:
                    gx = idf_fgfr3
                if 's249c' in key:
                    gx1 = idf_s249c

    # ---------------共现分析化学-------------------#
        if len(bladder) != 0 and bladder[0]:
            for ele in ChemicalNameList:
                if 'FGFR3' in ele['NameOfSubstance']:
                    gx = idf_fgfr3
                if 'S249C' in ele['NameOfSubstance']:
                    gx1 = idf_s249c
                    break
    # ---------------共现分析关键字-------------------#
        if len(bladder) != 0 and bladder[0]:
            for eleK in KeywordsList:
                if 's249c' in str(eleK).lower():
                    gx1 = idf_s249c
                if 'fgfr3' in str(eleK).lower():
                    gx = idf_fgfr3

     # ---------------共现分析医学主题词-------------------#
        if len(bladder) != 0 and bladder[0]:
            for eleM in MeshHeadingNameList:
                if 'S249C' in eleM['MeshHeadingName']:
                    gx1 = idf_s249c
                if 'FGFR3' in eleM['MeshHeadingName']:
                    gx = idf_fgfr3
                    break

        for ele in ChemicalNameList:
            if 'Urinary Bladder Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Rare Diseases' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Receptor, Fibroblast Growth Factor, Type 3' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break


        for eleM in MeshHeadingNameList:
            if 'Urinary Bladder Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Rare Diseases' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Receptor, Fibroblast Growth Factor, Type 3' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            if 'Male' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            if 'Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break


        for eleK in KeywordsList:
            if 'bladder' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 's249c' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break

        total_gx=gx1+gx2+gx3+gx
        cmk_len=len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
        bm25_cmk_len=ss1 + ss2 + ss4
        bm25_cmk_score = (((k2 + 1) * bm25_cmk_len) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + bm25_cmk_len))
        bm25_score=bm25_ab_score+bm25_cmk_score+total_gx
        if(bm25_score>yuzhi):
            mydict = {"PMID": x['PMID'],"ab_score":bm25_ab_score,"idf_para":idf_para,
                      "cmk_len":cmk_len,"cmk_freq":bm25_cmk_len,"bm25_cmk_score":bm25_cmk_score,"gx":total_gx,"bm25_score":bm25_score,
                       "ChemicalNameList":x['ChemicalNameList'],"MeshHeadingNameList":x['MeshHeadingNameList'],"KeywordsList":x['KeywordsList']}
            y = mydata.insert_one(mydict)
            k=k+1
            print(str(y) + '---------' + str(k))

def count(mysort,mycount,topic):
    for x in mysort.find({}, {'PMID', 'ab_score','idf_para', 'cmk_len', 'cmk_freq', 'bm25_cmk_score','gx','bm25_score',
                              'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
        kk = 0
        for y in mytopic.find({"topic": topic}, {'PMID', 'relate'}):
            if x['PMID'] == y['PMID']:
                mydict = {"PMID": x['PMID'], "related": y['relate'], "ab_score":x["ab_score"],"idf_para":x['idf_para'],
                          "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],'bm25_cmk_score':x['bm25_cmk_score'],'gx':x['gx'],
                          "bm25_score": x['bm25_score'],
                          "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                          "KeywordsList": x['KeywordsList']}
                ss = mycount.insert_one(mydict)
                print(ss)
                kk = kk + 1
        if (kk == 0):
            mydict = {"PMID": x['PMID'], "related": -1, "ab_score": x["ab_score"], "idf_para": x['idf_para'],
                      "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'bm25_cmk_score': x['bm25_cmk_score'],
                      'gx': x['gx'],
                      "bm25_score": x['bm25_score'],
                      "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                      "KeywordsList": x['KeywordsList']}
            ss = mycount.insert_one(mydict)
            print(ss)

if __name__ == '__main__':
    sortsecond(mywords,mydata,6)
    count(mydata,mycount,"8")



