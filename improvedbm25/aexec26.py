import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_26"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_26_related"]#聚类后对应与主题相关联的文献

def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_squamous = log((29138919 - 122216 + 0.5) / (122216 + 0.5), 10)
    idf_cell = log((29138919 - 5400577 + 0.5) / (5400577 + 0.5), 10)
    idf_lung = log((29138919 - 494482 + 0.5) / (494482 + 0.5), 10)
    idf_fgfr1 = log((29138919 - 36925 + 0.5) / (36925 + 0.5), 10)
    idf_amplification = log((29138919 - 2 + 0.5) / (2 + 0.5), 10)


    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0+ 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_6 = log((13670358 - 2126 + 0.5) / (2126 + 0.5), 10)


    idf_eleM_1 = log((25389659 - 46148 + 0.5) / (46148 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 200349 + 0.5) / (200349+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 85840 + 0.5) / (85840 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 16241 + 0.5) / (16241 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 123459 + 0.5) / (123459 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 2126 + 0.5) / (2126 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 55 + 0.5) / (55 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 222 + 0.5) / (222 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        gx4=0
        len_freq=0
        squamous_score=0
        cell_score = 0
        lung_score =0
        fgfr1_score = 0
        amplification_score=0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        squamous = [True for x in wordfreq.items() if 'squamous' in x]
        cell = [True for x in wordfreq.items() if 'cell' in x]
        lung = [True for x in wordfreq.items() if 'lung' in x]
        # ---------------摘要统计-------------------#


        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            if 'squamous' in key:
                squamous_score = squamous_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'cell' in key1:
                cell_score = cell_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'lung' in key1:
                lung_score = lung_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'fgfr1' in key1:
                fgfr1_score = fgfr1_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'amplification' in key1:
                amplification_score = amplification_score + wordfreq[key]

        bm25_squamous_score = (((k1+1)*squamous_score)/((k1*(b1+(1-b1)*(len_freq/85)))+squamous_score))
        bm25_cell_score = (((k1 + 1) * cell_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + cell_score))
        bm25_lung_score = (((k1 + 1) * lung_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + lung_score))
        bm25_fgfr1_score = (((k1 + 1) * fgfr1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + fgfr1_score))
        bm25_amplification_score = (((k1 + 1) * amplification_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + amplification_score))

        bm25_ab_score =idf_squamous*bm25_squamous_score+idf_cell*bm25_cell_score+idf_lung*bm25_lung_score+idf_fgfr1*bm25_fgfr1_score+idf_amplification*bm25_amplification_score

        idf_para=[{str(squamous_score):idf_squamous},{str(cell_score):idf_cell},{str(lung_score):idf_lung},{str(fgfr1_score):idf_fgfr1},{str(amplification_score):idf_amplification}]

        # ---------------共现分析摘要-------------------#
        if len(squamous) != 0 and squamous[0] and len(cell) != 0 and cell[0] and len(lung) != 0 and lung[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'fgfr1' in key:
                    gx = idf_fgfr1
                if 'amplification' in key:
                    gx1 = idf_fgfr1

    # ---------------共现分析化学-------------------#
        if len(squamous) != 0 and squamous[0] and len(cell) != 0 and cell[0] and len(lung) != 0 and lung[0]:
            for ele in ChemicalNameList:
                if 'FGFR1' in ele['NameOfSubstance']:
                    gx = idf_fgfr1
                    break

    # ---------------共现分析关键字-------------------#
        if len(squamous) != 0 and squamous[0] and len(cell) != 0 and cell[0] and len(lung) != 0 and lung[0]:
            for eleK in KeywordsList:
                if 'fgfr1' in str(eleK).lower():
                    gx = idf_fgfr1
                    break

     # ---------------共现分析医学主题词-------------------#
        if len(squamous) != 0 and squamous[0] and len(cell) != 0 and cell[0] and len(lung) != 0 and lung[0]:
            for eleM in MeshHeadingNameList:
                if 'FGFR1' in eleM['MeshHeadingName']:
                    gx = idf_fgfr1
                    break

        for ele in ChemicalNameList:
            if 'Carcinoma, Non-Small-Cell Lung' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Lung Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Carcinoma, Squamous Cell' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'Epithelial Cells' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            if 'Gene Amplification' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break
        for ele in ChemicalNameList:
            if 'Receptor, Fibroblast Growth Factor, Type 1' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_6
                break


        for eleM in MeshHeadingNameList:
            if 'Carcinoma, Non-Small-Cell Lung' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Lung Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Epithelial Cells' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if 'Gene Amplification' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_4
                break

        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            if 'Female' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            if 'Middle Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if 'Carcinoma, Squamous Cell' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_8
                break


        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break
        for eleM in MeshHeadingNameList:
            if 'Receptor, Fibroblast Growth Factor, Type 1' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_10
                break

        for eleK in KeywordsList:
            if 'squamous cell lung cancer' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'fgfr1' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break


        total_gx=gx1+gx2+gx3+gx+gx4
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
    sortsecond(mywords,mydata,14)
    count(mydata,mycount,"26")



