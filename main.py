# openpyxl needs to be installed
# pandas needs to be installed

import pandas
import psycopg2
import numpy as np

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
            sort = sort.loc[(sort[line] >= inputgreater) | (sort[line] <= -inputsmaller)]
            negativetop.append(sort[line].values[0:5])
            negativetopmetabolite.append(sort["HMDB_name"].values[0:5])
            positivetop.append(sort[line].values[-5:])
            positivetopmetabolite.append(sort["HMDB_name"].values[-5:])
            listpatient.append(line)
            sort["disease"] = sort["disease"].replace(np.nan, "None")
            negativedisease.append(sort["disease"].values[0:5])
            positivedisease.append(sort["disease"].values[-5:])
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
        "insert into patient (id,`test`) values"
            "('"+ str(count) + "','" + str(i) + "')")
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

    "DO NOT COMMENT THESE THINGS OUT. ONLY WHEN YOU WANT TO ADD THINGS TO THE DATABASE"
    # connection, cursor = connection_database()
    # patient_table(connection, cursor, patient_list)
    # z_score_table(connection, cursor, noduplicate_zscore)
    # disease_table(connection, cursor, noduplicate_disease)

main()