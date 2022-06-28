from lxml import etree
import pymongo
import os

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mycol = mydb["papers"]
xml_list=[]
path = ''
path_list=sorted(os.listdir('baseline'))
for filename in path_list:
    if os.path.splitext(filename)[1] == '.xml':
        xml_list.append('baseline/'+filename)
i=1

for xml in xml_list:
    tree = etree.parse(xml)
    articleslist=tree.xpath('//PubmedArticle')
    for articles in tree.xpath('//PubmedArticle'):
        PMID=str(articles.xpath('MedlineCitation/PMID/text()')[0])
        DateCompleted,DateRevised,JournalISSN,JournalIssueVolume,JournalIssueIssue='null','null','null','null','null'
        PubDateYear, PubDateMonth, PubDateMedlineDate, PubDateSeason, PubDateSuffix = 'null', 'null', 'null', 'null', 'null'
        ChemicalNameList, MeshHeadingNameList, KeywordsList, ReferencesList = [], [], [], []
        try:
            if(len(articles.xpath('MedlineCitation/Article/Journal/ISSN/text()'))!=0):
                JournalISSN=str(articles.xpath('MedlineCitation/Article/Journal/ISSN/text()')[0])
        except AttributeError:
            JournalISSN = 'null'
        try:
            if(len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/Volume/text()'))!=0):
                JournalIssueVolume=str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/Volume/text()')[0])
        except AttributeError:
            JournalIssueVolume = 'null'
        try:
            if(len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/Issue/text()')) != 0):
                JournalIssueIssue=str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/Issue/text()')[0])
        except AttributeError:
            JournalIssueIssue = 'null'

        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Year/text()')) != 0):
                PubDateYear=str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Year/text()')[0])
        except AttributeError:
            PubDateYear = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Month/text()')) != 0):
                PubDateMonth = str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Month/text()')[0])
        except AttributeError:
            PubDateMonth = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate/text()')) != 0):
                PubDateMedlineDate = str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate/text()')[0])
        except AttributeError:
            PubDateMedlineDate = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Season/text()')) != 0):
                PubDateSeason = str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Season/text()')[0])
        except AttributeError:
            PubDateSeason = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Suffix/text()')) != 0):
                PubDateSuffix = str(articles.xpath('MedlineCitation/Article/Journal/JournalIssue/PubDate/Suffix/text()')[0])
        except AttributeError:
            PubDateSuffix = 'null'


        JournalTitle='null'
        JournalISOAbbreviation = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/Title/text()')) != 0):
                JournalTitle = str(articles.xpath('MedlineCitation/Article/Journal/Title/text()')[0])
        except AttributeError:
            JournalTitle = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Journal/ISOAbbreviation/text()')) != 0):
                JournalISOAbbreviation = str(articles.xpath('MedlineCitation/Article/Journal/ISOAbbreviation/text()')[0])
        except AttributeError:
            JournalISOAbbreviation = 'null'


        ArticleTitle='null'
        PaginationMedlinePgn='null'
        ELocationID='null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/ArticleTitle/text()')) != 0):
                ArticleTitle = str(articles.xpath('MedlineCitation/Article/ArticleTitle/text()')[0])
        except AttributeError:
            ArticleTitle = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Pagination/MedlinePgn/text()')) != 0):
                PaginationMedlinePgn = str(articles.xpath('MedlineCitation/Article/Pagination/MedlinePgn/text()')[0])
        except AttributeError:
            PaginationMedlinePgn = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/ELocationID/text()')) != 0):
                ELocationID = str(articles.xpath('MedlineCitation/Article/ELocationID/text()')[0])
        except AttributeError:
            ELocationID = 'null'


        Abstract=''
        CopyrightInformation='null'
        AbstractTexts = articles.xpath('MedlineCitation/Article/Abstract/AbstractText/text()')
        try:
            for AbstractText in AbstractTexts:
                Abstract=Abstract+str(AbstractText)
        except AttributeError:
            Abstract = 'null'
        try:
            if(len(articles.xpath('MedlineCitation/Article/Abstract/CopyrightInformation/text()')) != 0):
                CopyrightInformation=str(articles.xpath('MedlineCitation/Article/Abstract/CopyrightInformation/text()')[0])
        except AttributeError:
            CopyrightInformation = 'null'


        AuthorNameList=[]
        Authors = articles.xpath('MedlineCitation/Article/AuthorList/Author')

        try:
            for Author in Authors:
                   author_name =  Author.find('LastName').text + ' ' + Author.find('ForeName').text
                   AuthorNameList.append({"AuthorName":str(author_name),"Initials":str(Author.find('Initials').text)})
        except:
            AuthorNameList=['null']


        AffiliationList = []
        try:
            for Author in Authors:
                affiliations = [x.xpath('Affiliation/text()')[0] for x in Author.xpath('AffiliationInfo')]
                for affiliation in affiliations:
                    AffiliationList.append({"Affiliation":str(affiliation)})
        except Exception as e:
            pass
            continue
        CollectiveName = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/AuthorList/Author/CollectiveName/text()')) != 0):
                CollectiveName = str(articles.xpath('MedlineCitation/Article/AuthorList/Author/CollectiveName/text()')[0])
        except AttributeError:
            CollectiveName = 'null'

        Language='null'
        try:
            if (len(articles.xpath('MedlineCitation/Article/Language/text()')) != 0):
                Language = str(articles.xpath('MedlineCitation/Article/Language/text()')[0])
        except AttributeError:
            Language = 'null'

        GrantNameList=[]
        Grants = articles.xpath('MedlineCitation/Article/GrantList/Grant')
        try:
            for Grant in  Grants:
                GrantNameList.append({"GrantID":str(Grant.find('GrantID').text),"Agency":str(Grant.find('Agency').text),"Country":Grant.find('Country').text})
        except AttributeError:
            GrantName='null'
        #print(PMID+','+GrantName)
        PublicationTypeNameList=[]
        PublicationTypes = articles.xpath('MedlineCitation/Article/PublicationTypeList/PublicationType/text()')
        try:
            for PublicationType in PublicationTypes:
                PublicationTypeNameList.append({"PublicationType":PublicationType})
        except AttributeError:
            PublicationTypeName = 'null'

        MedlineJournalInfoCountry='null'
        MedlineJournalInfoMedlineTA='null'
        MedlineJournalInfoNlmUniqueID='null'
        MedlineJournalInfoISSNLinking='null'
        try:
            if (len(articles.xpath('MedlineCitation/MedlineJournalInfo/Country/text()')) != 0):
                MedlineJournalInfoCountry = str(articles.xpath('MedlineCitation/MedlineJournalInfo/Country/text()')[0])
        except AttributeError:
            MedlineJournalInfoCountry = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/MedlineJournalInfo/MedlineTA/text()')) != 0):
                MedlineJournalInfoMedlineTA = str(articles.xpath('MedlineCitation/MedlineJournalInfo/MedlineTA/text()')[0])
        except AttributeError:
            MedlineJournalInfoMedlineTA = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/MedlineJournalInfo/lmUniqueID/text()')) != 0):
                MedlineJournalInfoNlmUniqueID = str(articles.xpath('MedlineCitation/MedlineJournalInfo/lmUniqueID/text()')[0])
        except AttributeError:
            MedlineJournalInfoNlmUniqueID = 'null'
        try:
            if (len(articles.xpath('MedlineCitation/MedlineJournalInfo/ISSNLinking/text()')) != 0):
                MedlineJournalInfoISSNLinking = str(articles.xpath('MedlineCitation/MedlineJournalInfo/ISSNLinking/text()')[0])
        except AttributeError:
            MedlineJournalInfoISSNLinking = 'null'
        Chemicals=articles.xpath('MedlineCitation/ChemicalList/Chemical')
        try:
            for Chemical in  Chemicals:
                ChemicalNameList.append({"RegistryNumber":Chemical.find('RegistryNumber').text,"NameOfSubstance": Chemical.find('NameOfSubstance').text})
        except AttributeError:
            ChemicalName = 'null'
        MeshHeadings = articles.xpath('MedlineCitation/MeshHeadingList/MeshHeading')
        try:
            for MeshHeading in  MeshHeadings:
                meshheading=MeshHeading.find('DescriptorName').text
                qualifiername=''
                for QualifierName in MeshHeading.findall('.//QualifierName'):
                    qualifiername=qualifiername+'&&'+QualifierName.text
                MeshHeadingNameList.append({"MeshHeadingName":meshheading,"QualifierName":qualifiername})
        except AttributeError:
            MeshHeadingName = 'null'
        OtherAbstract=''
        OtherAbstracts=articles.xpath('MedlineCitation/OtherAbstract/AbstractText/text()')
        try:
            for AbstractText in OtherAbstracts:
                OtherAbstract=OtherAbstract+str(AbstractText)
        except AttributeError:
            OtherAbstract = 'null'
        KeywordLists=articles.xpath('MedlineCitation/KeywordList/Keyword/text()')
        for Keyword in KeywordLists:
            KeywordsList.append({"KeyWord":Keyword})
        ReferenceList=articles.xpath('PubmedData/ReferenceList/Reference')
        try:
            for ref in ReferenceList:
                references=ref.find('Citation').text
                refids=ref.xpath('ArticleIdList/ArticleId/text()')
                ids=''
                for refid in refids:
                    ids=ids+'&&'+str(refid)
                ReferencesList.append({"references":references,"referencesids":ids})
        except AttributeError:
            References = 'null'
        mydict = {"PMID": PMID,
                  # "DateCompleted":DateCompleted,
                  # "DateRevised":DateRevised,
                  "JournalISSN":JournalISSN,
                  "JournalIssueVolume": JournalIssueVolume,
                  "JournalIssueIssue":JournalIssueIssue,
                  "PubDateYear":PubDateYear,
                  "PubDateMonth":PubDateMonth,
                  "PubDateMedlineDate":PubDateMedlineDate,
                  "PubDateSeason":PubDateSeason,
                  "PubDateSuffix":PubDateSuffix,
                  "JournalTitle":JournalTitle,
                  "JournalISOAbbreviation":JournalISOAbbreviation,
                  "ArticleTitle":ArticleTitle,
                  "PaginationMedlinePgn":PaginationMedlinePgn,
                  "ELocationID":ELocationID,
                  "Abstract":Abstract,
                  "CopyrightInformation":CopyrightInformation,
                  "AuthorNameList":AuthorNameList,
                  "AffiliationList":AffiliationList,
                  "CollectiveName":CollectiveName,
                  "Language":Language,
                  "GrantNameList":GrantNameList,
                  "PublicationTypeNameList":PublicationTypeNameList,
                  "MedlineJournalInfoCountry":MedlineJournalInfoCountry,
                  "MedlineJournalInfoMedlineTA":MedlineJournalInfoMedlineTA,
                  "MedlineJournalInfoNlmUniqueID":MedlineJournalInfoNlmUniqueID,
                  "MedlineJournalInfoISSNLinking":MedlineJournalInfoISSNLinking,
                  "ChemicalNameList":ChemicalNameList,
                  "MeshHeadingNameList":MeshHeadingNameList,
                  "OtherAbstract":OtherAbstract,
                  "KeywordsList":KeywordsList,
                  "ReferencesList":ReferencesList}

        x = mycol.insert_one(mydict)
        print(str(i)+','+PMID+','+str(x))
    i=i+1

