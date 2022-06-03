
# openpyxl needs to be installed
# pandas needs to be installed

from curses import meta
from os import PathLike
from textwrap import fill
import pandas
import psycopg2
import numpy as np
from Bio import Entrez

def readfile(document, inputsmaller, inputgreater):
    df = pandas.read_excel(document, index_col=None, na_values=['NA'])
    negativetop = []
    negativetopmetabolite = []
    positivetop = []
    positivetopmetabolite = []
    listpatient = []
    negativedisease = []
    positivedisease = []
    for line in df:
        if line.startswith("P") and line.endswith("Zscore"):
            sort = df.sort_values(line)
            sort = sort.loc[
                (sort[line] >= inputgreater) | (sort[line] <= -inputsmaller)]
            negativetop.append(sort[line].values[0:5])
            negativetopmetabolite.append(sort["HMDB_name"].values[0:20])
            positivetop.append(sort[line].values[-5:])
            positivetopmetabolite.append(sort["HMDB_name"].values[-20:])
            listpatient.append(line)
            sort["disease"] = sort["disease"].replace(np.nan, "None")
            negativedisease.append(sort["disease"].values[0:20])
            positivedisease.append(sort["disease"].values[-20:])
        else:
            continue

    return negativetop, negativetopmetabolite, positivetop, \
           positivetopmetabolite, listpatient, negativedisease, positivedisease


def dictionary(nt, ntm, pt, ptm, lp):
    patientlist = []
    list = []
    for p, i, i2, i3, i4 in zip(lp, nt, pt, ntm, ptm):
        patientlist.append(p)
        p = {}
        list.append(p)
        for l, m, in zip(i2, i4):
            p.update({m: l})
        for l, m in zip(i, i3):
            p.update({m: l})

    return patientlist, list


def zscore_no_duplicate(patient_values):
    value_no_duplicate = []
    for i in patient_values:
        for key, value in i.items():
            if value not in value_no_duplicate:
                value_no_duplicate.append(value)

    return value_no_duplicate


def disease_no_duplicates(patient_values):
    value_no_duplicate = []
    for i in patient_values:
        for key, value in i.items():
            value = value.replace("'", "")
            if ";" in value:
                value = value.split(";")[1:]
                for i in value:
                    i = i.lstrip()
                    if i not in value_no_duplicate:
                        value_no_duplicate.append(i)
            if "None" in value:
                if value not in value_no_duplicate:
                    value_no_duplicate.append(value)

    return value_no_duplicate


def connection_database():
    connection = psycopg2.connect(host="biocentre.nl",
                                  database="bio_jaar_2_pg_7",
                                  user="BI2_PG7",
                                  password="Blaat1234")
    cursor = connection.cursor()
    return connection, cursor


def patient_table(connection, cursor, patient_list):
    count = 0
    for i in patient_list:
        cursor.execute(
            "insert into patient (id, patient_code) values"
            "('" + str(count) + "','" + str(i) + "')")
        connection.commit()
        count += 1


def z_score_table(connection, cursor, value_no_duplicate):
    count = 0
    for i in value_no_duplicate:
        cursor.execute(
            "insert into z_score (id, z_score) values"
            "('" + str(count) + "','" + str(i) + "')")
        connection.commit()
        count += 1


def disease_table(connection, cursor, value_no_duplicate):
    count = 0
    for i in value_no_duplicate:
        cursor.execute(
            "insert into diseases (id, disease) values"
            "('" + str(count) + "','" + str(i) + "')")
        connection.commit()
        count += 1


def get_values(patient_dis_and_meta, connection, cursor, art_count):
    """ Deze functie haalt values uit de database en zet deze door om
    te text minen en de database vullen met pubmedIDs.
    """
    count = 0

    for i in patient_dis_and_meta:

        # print(count)
        for key, value in i.items():
            if ";" in value:
                value = value.split(";")[1:]
                for i in value:
                    i = i.lstrip()
                    art_count = compare(key, i, connection, cursor, art_count)
                count += 1
            else:
                continue

    return art_count

def compare(metaboliet, disease, connection, cursor, art_count):
    """ Deze funtie pakt het metaboliet en zoekt samen met disease naar
    pubmed ID's en geeft pubmed ID door om het in database te zetten.
    :param metaboliet: -str- metaboliet
    """
    Entrez.email = "A.N.Other@example.com"
    searchwords = metaboliet + " AND " + disease
    finds = Entrez.esearch(db="pubmed",
                           term=searchwords)
    result = Entrez.read(finds)

    # print(metaboliet)
    # print(disease)
    # print(result["IdList"])
    # return result["IdList"]
    metabolite_id_query = "select id from metabolite "\
        "where name_metabolite like %s"
    cursor.execute(metabolite_id_query, (metaboliet,))
    metaboliet_id = cursor.fetchall()[0][0]

    art_count = article_table(result["IdList"], connection, cursor, art_count, metaboliet_id)
    # add_data(metaboliet, disease, result["IdList"], cursor, z_score)
    return art_count

def article_table(pubmed_id, connection, cursor, art_count, metaboliet_id):
    cursor.execute("select pubmed_id from articale_metabolite")
    list_articles = []
    try:
        for value in cursor.fetchall():
            list_articles.append(value[0])
    except IndexError:
        print("list empty just yet")

    
    for i in pubmed_id:
        if i in list_articles:
            pubmed_id_query = "select id from articale_metabolite "\
                "where pubmed_id like %s"
            pub_id = str(i)
            cursor.execute(pubmed_id_query, (pub_id,))
            pubmed_id = cursor.fetchall()[0][0]
            fill_table(metaboliet_id, pubmed_id, connection, cursor)
            
        else:
            cursor.execute(
                "insert into articale_metabolite (id, pubmed_id) values"
                "('" + str(art_count) + "','" + str(i) + "')")
            connection.commit()
            fill_table(metaboliet_id, art_count, connection, cursor)
            art_count += 1
    return art_count


def fill_table(metabolite_id, pubmed_id,connection, cursor):
    check_query = "select * from articale_metabolite_metabolite where articale_metabolite_id = %s and metabolite_id = %s"
    cursor.execute(check_query,(pubmed_id, metabolite_id,))
    if not cursor.fetchall():
        tussen_query ="insert into articale_metabolite_metabolite(articale_metabolite_id, metabolite_id) "\
                            " values(%s, %s)"
        cursor.execute(tussen_query, (pubmed_id, metabolite_id,))
        connection.commit()


def main():
    document = "Output_untargeted_metabolomics.xlsx"
    inputgreater = float(1.2)
    inputsmaller = float(1.2)

    # This function is to get z_scores, metabolites, diseases and patients
    negativetop, negativetopmetabolite, positivetop, positivetopmetabolite, \
    listpatient, negativedisease, positivedisease = \
        readfile(document, inputsmaller, inputgreater)

    # This fuction is to get the patientlist filtered and get a z_score with
    # metabolite dictionary
    patient_list, patient_z_and_meta = \
        dictionary(negativetop, negativetopmetabolite,
                   positivetop, positivetopmetabolite, listpatient)

    # This fuction is to get the patientlist filtered and get a disease with
    # metabolite dictionary
    patient_list, patient_dis_and_meta = \
        dictionary(negativedisease, negativetopmetabolite,
                   positivedisease, positivetopmetabolite, listpatient)

    # This fuction is for a no duplicate list of z_scores
    noduplicate_zscore = zscore_no_duplicate(patient_z_and_meta)
    # This fuction is for a no duplicate list of diseases
    noduplicate_disease = disease_no_duplicates(patient_dis_and_meta)
    # print(noduplicate_disease)
    # connection, cursor = connection_database()
    # cursor = ""

    # print(ham)

    art_count = 0
    "DO NOT COMMENT THESE THINGS OUT. ONLY WHEN YOU WANT TO ADD THINGS TO THE DATABASE"
    connection, cursor = connection_database()
    # patient_table(connection, cursor, patient_list)
    # z_score_table(connection, cursor, noduplicate_zscore)
    # disease_table(connection, cursor, noduplicate_disease)
    get_values(patient_dis_and_meta, connection, cursor, art_count)




if __name__ == "__main__":
    main()