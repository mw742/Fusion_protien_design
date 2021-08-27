#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 17:39:40 2021

@author: wmm
"""
#################################################################################
## For a protein dataset which can all be found on Uniprot,
## the program will extact all the loop areas information of these proteins,
## and save as an excel file
#################################################################################

import pandas as pd
import requests
import xlwt
from bs4 import BeautifulSoup
import urllib3
urllib3.disable_warnings()






#######################################################################
##输入蛋白文件，程序就会提取相应的uniprot库中的id，
##之后用于从中爬取相应id蛋白信息中的linker信息
#######################################################################

#打开数据库文件protein_class_Predicted.tsv，
#从中依次提取库中每一个蛋白对应的uniprot id，
#便于下一步在uniprot库中提取需要的信息
protein_data=pd.read_csv('protein_class_Predicted.tsv', sep='\t')
data_length=len(protein_data)
uniprot_id=[]
i=0
while i < data_length:
    #有一些蛋白没有储存在uniprot库里，对应的uniprot列没有信息
    #有记录的蛋白的uniprot值为字符串类型，因此可以从中提取出来
    if isinstance(protein_data["Uniprot"][i], str):
        uniprot_id.append(protein_data["Uniprot"][i])
    else:
        pass
    i+=1
#print(uniprot_id)






#######################################################################
##依次提取上文收集的id list中对应的蛋白的信息
#######################################################################

#利用上面提取的uniprot id，依次从uniprot网站提取每个蛋白对应的信息
#并把提取到的linker序列储存在linker_database列表里
linker_database=[]
#把表格中的linker提取到一个excel表里
#创建一个excel表格
#第一个sheet为sheet1
work_book=xlwt.Workbook(encoding='utf-8')
sheet=work_book.add_sheet('sheet1', cell_overwrite_ok=True)
n=0
#开始挨个查找
for ID in uniprot_id:
    print("uniprot_id", ID)
    #依次输入搜索的id，获取相应的网页源码
    html=requests.get("https://www.uniprot.org/uniprot/"+ID, headers={"Connection":"close"}, verify=False)
    res=html.text
    #再次封装，利用beautifulsoup包解析获取具体标签内的内容
    bs=BeautifulSoup(res,'lxml')
    #找出所有表格中，行的数据
    tr=bs.select('a.position.tooltipped')
    #挨个遍历每一行，利用关键字domain定位要找的数值
    #把domain序列位置的信息储存在domain_position里
    domain_positions=[]
    for item in tr:
        htmlitem_str=str(item)
        if htmlitem_str.find("key=Domain") > -1:
            domain_number=htmlitem_str.split(">")[1].split("<")[0]
            domain_positions.append(domain_number)
        else:
            pass
    #如果只有一个doomain或者没有domain，则认为结果中没有loops
    #否则，loops为连接两两个domains之间的部分
    if len(domain_positions) <= 1:
        pass
    else:
        #在domain序列数中，收集连接它们之间的序列的数
        anchors=[]
        total_order=len(domain_positions)
        current_order=0
        while current_order < total_order:
            left_anchor=int(domain_positions[current_order].split("–")[0])
            right_anchor=int(domain_positions[current_order].split("–")[1])
            anchors.append(left_anchor)
            anchors.append(right_anchor)
            current_order+=1
        #除掉首尾的数
        anchors=anchors[1:-1]
        #用linker_positions列表，把linker的首尾序列数储存起来，便于后期从整个蛋白序列中提取linker序列
        pairs=int(len(anchors)/2)
        current_pair=0
        left_pos=0
        linker_positions=[]
        while current_pair < pairs:
            linker_left_pos=int(anchors[left_pos]+1)
            linker_right_pos=int(anchors[left_pos+1]-1)
            if linker_right_pos <= linker_left_pos:
                pass
            else:
                linker_positions.append([linker_left_pos, linker_right_pos])
            left_pos+=2
            current_pair+=1
        #print(linker_positions)
        #读取网页中sequences部分的序列信息
        if linker_positions:
            #如果该蛋白中，domains之间存在linker的位置，利用上面提取到的linker position，从整条序列中提取linker序列
            sequence_to_be_crawlered=requests.get("https://www.uniprot.org/uniprot/"+ID+".fasta", headers={"Connection":"close"}, verify=False)
            sequence_info_text=sequence_to_be_crawlered.text
            sequence=''.join(sequence_info_text.split("\n")[1:])
            total_linker_number=len(linker_positions)
            current_linker_number=0
            while current_linker_number < total_linker_number:
                linker_left=linker_positions[current_linker_number][0]
                linker_right=linker_positions[current_linker_number][1]
                linker_sequence=sequence[linker_left:linker_right+1]
                linker_region=str(linker_left)+"-"+str(linker_right)
                linker_length=len(linker_sequence)
                protein_source=ID
                if len(linker_sequence) < 50 & len(linker_sequence) > 1:
                    linker_database.append(linker_sequence)
                    sheet.write(n,0,n)
                    sheet.write(n,1,linker_sequence)
                    sheet.write(n,2,protein_source)
                    sheet.write(n,3,linker_region)
                    sheet.write(n,4,linker_length)
                    #把linker都存成FASTA格式便于后期处理
                    with open("Linkers_fasta.txt", "a+") as f:
                        f.write(linker_sequence)
                        f.write("\n")
                        f.write("\n")
                    f.close()                        
                    n+=1
                else:
                    pass
                current_linker_number+=1
                print("linker_database:", linker_database)
                #print(protein_source)
                #print(linker_region)                
        else:
            pass
#print(linker_database)





#################################################################################
##得到linker基础序列信息后，调取常用多肽/蛋白计算api，计算linker其他结构信息并储存
#################################################################################

#利用GOR4二级结构预测工具计算获得的linker的二级结构
with open("Linkers_fasta.txt", "r") as f:
    for line in f:
        n=0
        if len(line) > 1:
            data="title=&notice="+line+"&ali_width=70"
            headers={"Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9",
                     "Accept-Encoding": "gzip, deflate, br",
                     "Accept-Language": "zh-CN,zh;q=0.9",
                     "Cache-Control": "max-age=0",
                     "Connection": "keep-alive",
                     "Content-Length": "96",
                     "Content-Type": "application/x-www-form-urlencoded",
                     "Host": "npsa-prabi.ibcp.fr",
                     "Origin": "https://npsa-prabi.ibcp.fr",
                     "Referer": "https://npsa-prabi.ibcp.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_gor4.html",
                     "sec-ch-ua": '"Chromium";v="92", " Not A;Brand";v="99", "Google Chrome";v="92"',
                     "sec-ch-ua-mobile": "?0",
                     "Sec-Fetch-Dest": "document",
                     "Sec-Fetch-Mode": "navigate",
                     "Sec-Fetch-Site": "same-origin",
                     "Sec-Fetch-User": "?1",
                     "Upgrade-Insecure-Requests": "1",
                     "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36"
                     }
            url="https://npsa-prabi.ibcp.fr/cgi-bin/secpred_gor4.pl"
            res=requests.post(url, data=data, headers=headers)
            res_text=res.text
            bs=BeautifulSoup(res_text,'lxml')
            secondary_info=bs.select("code")[0].text.split("\n")[4]
            print("secondary_info:", secondary_info)
            sheet.write(n,5,secondary_info)
            c_count=secondary_info.count("c")
            seq_length=len(secondary_info)
            c_rate=float(c_count/seq_length)
            sheet.write(n,6,secondary_info)
            sheet.write(n,7,c_rate)
            if c_rate <= 0.3:
                sheet.write(n,8,"True")
            else:
                sheet.write(n,8,"False")
            n+=1
f.close()
work_book.save('Linkers_from_secreted.xls')



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    