import MySQLdb
import shelve

#CREATE TABLE organisms (organism_id INTEGER(5) PRIMARY KEY AUTO_INCREMENT, organism_name VARCHAR(300), transcriptome_list VARCHAR(500),gwips_databasename VARCHAR(300),gwips_clade VARCHAR(300),gwips_organism VARCHAR(300), gwips_database VARCHAR(300), default_transcript VARCHAR(300));
#CREATE TABLE studies (study_id INTEGER(6) PRIMARY KEY AUTO_INCREMENT, organism_id INTEGER(5), study_name VARCHAR(500), paper_authors VARCHAR(300), srp_nos VARCHAR(300), paper_year VARCHAR(300), paper_pmid VARCHAR(300), paper_link VARCHAR(300), gse_nos VARCHAR(300), adapters VARCHAR(300), paper_title VARCHAR(300), description VARCHAR(2000), private TINYINT(1));
#CREATE TABLE files (file_id INTEGER(6) PRIMARY KEY AUTO_INCREMENT, organism_id INTEGER(5), study_id INTEGER(6), file_name VARCHAR(300), file_description VARCHAR(300), file_type VARCHAR(300));
#CREATE TABLE users (user_id INTEGER(6) PRIMARY KEY AUTO_INCREMENT, username VARCHAR(300), password VARCHAR(300), study_access VARCHAR(300));
#CREATE TABLE urls (url_id INTEGER(6) PRIMARY KEY AUTO_INCREMENT, url VARCHAR(2000));
#CREATE TABLE user_settings (user_id INTEGER(6) PRIMARY KEY AUTO_INCREMENT, background_col VARCHAR(300), readlength_col VARCHAR(300), metagene_fiveprime_col VARCHAR(300),metagene_threeprime_col VARCHAR(300));

#INSERT INTO user_settings VALUES(1,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(2,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(3,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(4,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(5,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(6,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");
#INSERT INTO user_settings VALUES(7,"#F2F2F7","#FF5F5B","#FF5F5B","#9ACAFF");

connection = MySQLdb.connect(host="localhost",
                             user="stephenk",
                             passwd="ribosome",
                             db="trips")

cursor = connection.cursor()
"""

cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(1,"saccharomyces_cerevisiae","Gencode_v24","yeast","yeast","S.+cerevisiae","sacCer3","YPR122W"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(2,"mycoplasma_hyorhinis","Ensemble_release37","mycoplasma_hyorhinis","mycoplasma_hyorhinis","M.+hyorinis","mh","AFX74565"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(3,"mus_musculus","Gencode_M14","mouse","mammal","Mouse","mm10","ENSMUST00000037796"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(4,"homo_sapiens","Gencode_v25","homo_sapiens","mammal","Human","hg38","ENST00000437161"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(5,"homo_sapiens_polio","Ensemble_2011","homo_sapiens","mammal","Human","hg38","polio"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(6,"rattus_norvegicus","Ensembl_r6","rattus_norvegicus","mammal","Rat","rn6","ENSRNOT00000047550"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(7,"drosophila_melanogaster","Flybase","drosophila_melanogaster","insect","Fruitfly","dm3","NM_134601"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(8,"escherichia_coli","Ensemble_k_12_ASM584v2","eschColi_K12","bacteria","Escherichia+coli+K12","eschColi_K12","gapA"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(9,"oyster","forward_transcriptome","oyster","mammal","Oyster","oyster","CHOYP_MATN1.3.5"))
cursor.execute("INSERT INTO organisms VALUES ({},'{}','{}','{}','{}','{}','{}','{}')".format(10,"caenorhabditis_elegans","Ensemble_2011","ce10","unknown","unknown","ce10","CED-4"))
connection.commit()





studyshelve = shelve.open("/home/DATA/GWIPS_viz/Trips_shelves/study_info.shelf")
allfiles_shelve = shelve.open("/home/DATA/GWIPS_viz/Trips_shelves/all_files.shelf")


study_id = 1


private_list = ["kolupaeva17","kolupaeva_resequenced_v2_feb16","Kolupaeva_resequence_apr19",
                "kolupaeva_27_sep17_rat","sonenberg17",
                "sonenberg17_human_only","Andreev_I_series","kolupaeva_27_sep17_human",
                "Andreev_K_series","Andreev_HCT116","Andreev_I_series",
                "Andreev_interferon","Andreev_GH","Andreev_new_ef","Baclaocos17_oyster"]




for org in allfiles_shelve:
    cursor.execute("SELECT organism_id from organisms WHERE organism_name = '{}';".format(org))
    org_id = cursor.fetchone()
    if org_id == None:
        continue
    org_id = org_id[0]

    org_id = int(org_id)
    for studyname in allfiles_shelve[org]:
        
        paper_authors = None
        srp_nos = None
        paper_year = None
        paper_pmid = None
        paper_link =None
        gse_nos = None
        adapters = None
        paper_title = None
        desc = None
        
        
        if studyname in private_list:
            private = 1
        else:
            private = 0



        
        if studyname in studyshelve[org]:
            if "paper_authors" in studyshelve[org][studyname]:
                paper_authors = str(studyshelve[org][studyname]["paper_authors"]).replace("'","").strip("[]")
    
            if "srp_nos" in studyshelve[org][studyname]:
                srp_nos = str(studyshelve[org][studyname]["srp_nos"]).replace("'","").strip("[]")
    
            if "paper_year" in studyshelve[org][studyname]:
                paper_year = studyshelve[org][studyname]["paper_year"]    

            if "paper_pmid" in studyshelve[org][studyname]:
                paper_pmid = studyshelve[org][studyname]["paper_pmid"]    

            if "paper_link" in studyshelve[org][studyname]:
                paper_link = (studyshelve[org][studyname]["paper_link"])
                if paper_link != None:
                    paper_link = paper_link.strip(" ")

            if "gse_nos" in studyshelve[org][studyname]:
                gse_nos = studyshelve[org][studyname]["gse_nos"]
                if gse_nos != None:
                    gse_nos = gse_nos.replace("None","")
            if "adapters" in studyshelve[org][studyname]:
                adapters = str((studyshelve[org][studyname]["adapters"])).replace("'","").strip("[]").strip(" ")

            if "paper_title" in studyshelve[org][studyname]:
                paper_title = studyshelve[org][studyname]["paper_title"]

            if "desc" in studyshelve[org][studyname]:
                
                desc = (studyshelve[org][studyname]["desc"])
                if desc != None:
                    desc = desc.strip(" ")
    
        #print "INSERT INTO studies VALUES ({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}')".format(study_id, org_id,studyname,paper_authors,srp_nos,paper_year, paper_pmid, paper_link, gse_nos, adapters, paper_title, desc)
        cursor.execute("INSERT INTO studies VALUES ({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{})".format(study_id, org_id,studyname,paper_authors,srp_nos,paper_year, paper_pmid, paper_link, gse_nos, adapters, paper_title, desc,private))
        study_id += 1
        
connection.commit()        
"""

studyshelve = shelve.open(
    "/home/DATA/GWIPS_viz/Trips_shelves/study_info.shelf")
allfiles_shelve = shelve.open(
    "/home/DATA/GWIPS_viz/Trips_shelves/all_files.shelf")

##CREATE TABLE files (file_id INTEGER(6) PRIMARY KEY, organism_id INTEGER(5), study_id INTEGER(6), file_name VARCHAR(300), file_description VARCHAR(300), file_type VARCHAR(300));

file_id = 1

for org in allfiles_shelve:
    print(org)
    cursor.execute(
        "SELECT organism_id from organisms WHERE organism_name = '{}';".format(
            org))
    org_id = cursor.fetchone()
    if org_id == None:
        continue
    org_id = org_id[0]

    org_id = int(org_id)

    #print allfiles_shelve[org]

    for studyname in allfiles_shelve[org]:
        riboseq_filenames = []
        rnaseq_filenames = []
        cursor.execute(
            "SELECT study_id from studies WHERE study_name = '{}' and organism_id = {};"
            .format(studyname, org_id))
        study_id = cursor.fetchone()
        if study_id == None:
            continue
        study_id = int(study_id[0])
        for filename in allfiles_shelve[org][studyname]["riboseq"]:
            if filename not in riboseq_filenames:
                cursor.execute(
                    "INSERT INTO files VALUES({},{},{},'{}','{}','{}')".format(
                        file_id, org_id, study_id, filename, "None",
                        "riboseq"))
                file_id += 1
                riboseq_filenames.append(filename)
        for filename in allfiles_shelve[org][studyname]["rnaseq"]:
            if filename not in rnaseq_filenames:
                cursor.execute(
                    "INSERT INTO files VALUES({},{},{},'{}','{}','{}')".format(
                        file_id, org_id, study_id, filename, "None", "rnaseq"))
                file_id += 1
                rnaseq_filenames.append(filename)

connection.commit()
