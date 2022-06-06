# openpyxl, pandas, Bio moet geïnstalleerd zijn.

import pandas
import psycopg2
import numpy as np
from Bio import Entrez

def leesbestand(document, inputsmaller, inputgreater, inputaantal):
    """Alle waardes die verder gebruikt worden, oftewel de HMDB_name,
    disease, patiënten worden uit het bestand gehaald en opgeslagen in lijsten.


    :param document: Een xlsx bestand
    :param inputsmaller: input van een negatieve Zscore
    :param inputgreater: input van een positieve Zscore
    :param inputaantal: aantal resultaten dat uit het xlsx bestand gehaalt moet worden.
    :return: zeven lijsten met waardes worden geretouneerd
    """
    # Bestand inlezen met behulp van pandas.
    df = pandas.read_excel(document, index_col=None, na_values=['NA'])
    # Zeven lijsten voor de waardes die opgeslagen moeten worden
    negatievetop = []
    negatievetopmetaboliet = []
    positievetop = []
    positievetopmetaboliet = []
    lijstpatienten = []
    negatieveziekte = []
    positieveziekte = []
    # for loop om over de kolom namen te gaan
    for kolom in df:
        # Als kolomnaam begint met een P en eindigt op Zscore, wordt de data opgeslagen
        if kolom.startswith("P") and kolom.endswith("Zscore"):
            sort = df.sort_values(kolom)
            sort = sort.loc[(sort[kolom] >= inputgreater) | (sort[kolom] <= -inputsmaller)]
            negatievetop.append(sort[kolom].values[0:inputaantal])
            negatievetopmetaboliet.append(sort["HMDB_name"].values[0:inputaantal])
            positievetop.append(sort[kolom].values[-inputaantal:])
            positievetopmetaboliet.append(sort["HMDB_name"].values[-inputaantal:])
            lijstpatienten.append(kolom)
            sort["disease"] = sort["disease"].replace(np.nan, "None")
            negatieveziekte.append(sort["disease"].values[0:inputaantal])
            positieveziekte.append(sort["disease"].values[-inputaantal:])
        # Als kolomnaam niet begint met een P en/of eindigt op Zscore, gebeurt er niks
        else:
            continue

    return negatievetop, negatievetopmetaboliet, positievetop, \
           positievetopmetaboliet, lijstpatienten, negatieveziekte, positieveziekte


def dictionary(nt, ntm, pt, ptm, lp):
    """ De negatievetopmetaboliet, positievetopmetaboliet en patientenlijst moeten hetzelfde
    zijn als de functie vaker wordt aangeroepen. Er wordt door alle vijf hierboven genoemde
    lijsten heen gelooped. Hierbij worden alle patienten in een lijst gezet. Daarna wordt
    de patientennaam de naam van de dictionary. Tot slot worden metabolieten gekoppeld als keys
    aan de andere meegegeven waarden als values.

    :param nt: een willekeurige negatievetop
    :param ntm: negatievetopmetaboliet
    :param pt: een willekeurige positievetop
    :param ptm: positievetopmetaboliet
    :param lp: patiëntenlijst
    :return: patientenlijst & lijst_met_dictionarys
    """
    patientenlijst = []
    lijst_met_dictionarys = []
    # Loopt tegelijk door 5 verschillende lijsten heen
    for p, i, i2, i3, i4 in zip(lp, nt, pt, ntm, ptm):
        patientenlijst.append(p)
        p = {}
        lijst_met_dictionarys.append(p)
        # Voegt positievetopmetabolieten toe als key met een andere positievetop als value
        for l, m, in zip(i2, i4):
            p.update({m: l})
        # Voegt negatievetopmetabolieten toe als key met een andere negatievetop als value
        for l, m in zip(i, i3):
            p.update({m: l})

    return patientenlijst, lijst_met_dictionarys


def disease_no_duplicates(patienten_waardes):
    """Er wordt met behulp van for loops gekeken of bepaalde waardes niet meervoudig voorkomen.
    Dit wordt gedaan om redunanten waardes niet op te slaan. Ook worden de ids die later in database
    staan gekoppelt aan de ziektes die niet dubbel zijn.

    :param patienten_waardes: lijst met dictionary van metabolieten als key en ziektes als value
    :return:value_geen_duplicatie & id_dict
    """
    value_geen_duplicatie = []
    count = -1
    id_dict = {}
    # Loop door de lijst hierbij is i iedere keer een nieuwe dictionary
    for i in patienten_waardes:
        # Per dictionary worden de keys (metabolieten) en values (ziektes) opgehaald
        for key, value in i.items():
            # ' moet uit de values gehaald worden om foutmeldingen te verkomen
            value = value.replace("'", "")
            # Er wordt gekeken of ; in values staat
            if ";" in value:
                # Er wordt gesplits op ; en de eerste waarde wordt weggehaald aangezien dit
                # een spatie is geworden afgesloten met een komma
                value = value.split(";")[1:]
                # Er wordt gekeken naar iedere ziekte die in de value staat
                for i in value:
                    # Alle spaties worden aan de linker kant weggehaald
                    i = i.lstrip()
                    # Er wordt gekeken of de ziekte niet al een keer voorkomt in de lijst
                    if i not in value_geen_duplicatie:
                        count += 1
                        id_dict[i] = count
                        value_geen_duplicatie.append(i)
            # Er wordt gekeken of None in values staat
            if "None" in value:
                # Er wordt gekeken of de ziekte niet al een keer voorkomt in de lijst
                if value not in value_geen_duplicatie:
                    count += 1
                    id_dict[value] = count
                    value_geen_duplicatie.append(value)

    return value_geen_duplicatie, id_dict


def metabolite_no_duplicates(patient_waardes):
    """Er wordt met behulp van for loops gekeken of een bepaald metaboliet al voorkomt
    in een lijst. Als hij nog niet voorkomt dan wordt deze metaboliet toegevoegd aan de lijst.

    :param patient_waardes: Een willekeurige lijst met dictionarys waar metaboliet
    de key is
    :return: value_geen_duplicatie & id_dict
    """
    count = -1
    id_dict = {}
    value_geen_duplicatie = []
    # Loopt door alle dictionarys heen die aanwezig zijn in de lijst
    for i in patient_waardes:
        # Per dictionary worden de keys (metabolieten) en random values opgehaald
        for key, value in i.items():
            # Kijkt of de metaboliet al voorkomt in de lijst
            if key not in value_geen_duplicatie:
                count += 1
                id_dict[key] = count
                value_geen_duplicatie.append(key)

    return value_geen_duplicatie, id_dict


def connection_database():
    """Er wordt een connectie gemaakt met de database

    :return: connection, cursor
    """
    # Connectie maken met de database
    connection = psycopg2.connect(host="biocentre.nl",
                                  database="bio_jaar_2_pg_7",
                                  user="BI2_PG7",
                                  password="Blaat1234",
                                  port="5900")
    # Een cursor maken
    cursor = connection.cursor()
    return connection, cursor


def patient_table(connection, cursor, patient_list):
    """De tafel met patiënten wordt gevuld

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param patient_list: lijst met alle patiënten
    :return: database patient tafel gevuld
    """
    count = 0
    # loop door de lijst met patiënten om ze een voor een toe te voegen aan de tafel
    for i in patient_list:
        cursor.execute(
        "insert into patient (id, patient_code) values"
            "('"+ str(count) + "','" + str(i) + "')")
        connection.commit()
        count += 1


def disease_table(connection, cursor, ziekte_geen_duplicatie):
    """De tafel met ziektes wordt gevuld

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param ziekte_geen_duplicatie: lijst met ziektes zonder duplicaties
    :return: database disease tafel gevuld
    """
    count = 0
    # loop door de lijst met ziektes zonder duplicates om ze een voor een toe te voegen
    # aan de tafel.
    for i in ziekte_geen_duplicatie:
        cursor.execute(
            "insert into diseases (id, disease) values"
            "('" + str(count) + "','" + str(i) + "')")
        connection.commit()
        count += 1


def metabolite_table(connection, cursor, meta_geen_duplicatie):
    """De tafel met metabolieten wordt gevuld

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param meta_geen_duplicatie: lijst met metabolieten zonder duplicaties
    :return: database metabolite tafel gevuld
    """
    count = 0
    for i in meta_geen_duplicatie:
        cursor.execute(
            "insert into metabolite (id, metabolite) values"
            "('" + str(count) + "','" + str(i.replace("'", "")) + "')")
        connection.commit()
        count += 1


def patient_metabolite_table(connection, cursor, list_patient_id, list_metabolite_id, list_z_scores):
    """De tussentafel van patiënten en metabolieten worden gevuld en gekoppeld aan elkaar. Hierbij
    worden ook de z_scores toegevoegd.

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param list_patient_id: ids van de patienten voor de foreignkeys
    :param list_metabolite_id: ids van de metabolieten voor de foreignkeys
    :param list_z_scores: lijst met z_scores
    :return: database patient_metabolite tussentafel gevuld
    """
    # loopt met behulp van een zip door alle drie de lijsten tegelijk heen
    # De waardes worden een voor een toegevoegd aan de database
    for i, i2, i3 in zip(list_z_scores, list_patient_id, list_metabolite_id):
        cursor.execute(
            "insert into metabolite_patient (z_score, patient_id, metabolite_id) values"
            "('" + str(i) + "','" + str(i2) + "','" + str(i3) + "')")
        connection.commit()


def metabolite_disease_table(connection, cursor, metaboliet_id, ziekte_id, resultaat):
    """De tussentafel van metabolieten en ziektes worden gevuld en gekoppeld aan elkaar.
    Hierbij worden ook de pubmid artikel ids toegevoegd.

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param metaboliet_id: ids van de metabolieten voor de foreignkeys
    :param ziekte_id: ids van de zietktes voor de foreignkeys
    :param resultaat: lijst gevuld met pmids
    :return: database metabolite_disease tussentafel gevuld
    """
    # loop door de lijst met pmids om ze een voor een toe te voegen aan de tafel
    for i in resultaat:
        cursor.execute(
            "insert into metabolite_diseases (article, metabolite_id, diseases_id) values"
            "('" + str(i) + "','" + str(metaboliet_id) + "','" + str(ziekte_id) + "')")
        connection.commit()


def metaboliet_ziekte_geen_duplicatie(patient_dis_and_meta):
    """Kijkt per key en value of de value (ziekte) al gekoppeld staat aan de key (metaboliet).
    Als dit niet het geval is wordt er een nieuwe dictionary gemaakt en de redunantie eruit
    gehaald. Als de key er al instaat wordt de loop overgeslagen.

    :param patient_dis_and_meta:
    :return:
    """
    geen_duplicatie_meta_ziekte = {}
    # Loopt door alle dictionarys heen die aanwezig zijn in de lijst
    for i in patient_dis_and_meta:
        # Per dictionary worden de keys (metabolieten) en de values (ziektes) opgehaald
        for key, value in i.items():
            # Er wordt gekeken of ; in values staat
            if ";" in value:
                # Als de metaboliet al voorkomt in de lijst wordt de volgende for loop begonnen
                if key in geen_duplicatie_meta_ziekte:
                    continue
                # Er wordt gesplitst op ; en de eerste waarde wordt weggehaald aangezien dit
                # een spatie is geworden afgesloten met een komma
                value = value.split(";")[1:]
                for i in value:
                    # Alle spaties worden aan de linkse kant weggehaald
                    i = i.lstrip()
                    # Kijkt of de metaboliet al voorkomt in de lijst
                    if key in geen_duplicatie_meta_ziekte:
                        geen_duplicatie_meta_ziekte[key] += [i]
                    # Als de metaboliet nog niet voorkomt wordt er een nieuwe key en value gemaakt
                    else:
                        geen_duplicatie_meta_ziekte[key] = [i]
            # Als de metaboliet nog niet voorkomt wordt er een nieuwe key en value gemaakt
            else:
                geen_duplicatie_meta_ziekte[key] = [value]
                continue

    return geen_duplicatie_meta_ziekte


def ids_koppelen(metaboliet_ziekte_geen_duplicatie, id_ziektes, id_metabolieten):
    """

    :param metaboliet_ziekte_geen_duplicatie: dictionary met metabolieten en ziektes zonder duplicaties
    :param id_ziektes: dictionary met ziektes als key en indexen als value
    :param id_metabolieten: dictionary met metabolieten als key en indexen als value
    :return: metabolieten_ziektes_ids
    """
    metabolieten_ziektes_ids = {}
    # Loopt door de dictionary met non redunanten data van metabolieten en ziektes
    for key, value in metaboliet_ziekte_geen_duplicatie.items():
        # Loopt door de dictionary met metabolieten gekoppelt met indexen
        for key2, value2 in id_metabolieten.items():
            # Kijkt of de metaboliet in de id dictionary gelijk staat aan de metaboliet
            # van de non redunanten data dictionary
            if key2 == key:
                # Loopt door de dictionary met ziektes gekoppelt met indexen
                for key3, value3 in id_ziektes.items():
                    # Loop om alle verschillende ziektes uit de lijst te krijgen bij de non
                    # redunanten data
                    for i in value:
                        # Kijkt of de ziekte in de id dictionary gelijk staat aan de ziekte
                        # van de non redunanten data dictionary
                        if key3 == i:
                            # Kijkt of de index al als key wordt gebruikt in de metabolieten_ziektes_ids
                            # dictionary
                            if value2 in metabolieten_ziektes_ids:
                                metabolieten_ziektes_ids[value2] += [value3]
                                # Kijkt of de index nog niet als key wordt gebruikt in de metabolieten_ziektes_ids
                                # dictionary
                            elif value2 not in metabolieten_ziektes_ids:
                                metabolieten_ziektes_ids[value2] = [value3]
                        else:
                            continue

    return metabolieten_ziektes_ids


def filter_z_patient_meta(id_metabolieten, patient_waardes):
    """ Er worden drie lijsten aangemaakt en er wordt gekeken of de metabolieten uit de
    twee dictionary hetzelfde zijn. Als dit het geval is worden alle ids toegevoegd aan een lijst
    en worden de z_scores toegevoegd aan de lijst_z_scores.

    :param id_metabolieten: dictionary met metabolieten en ids
    :param patient_waardes: Lijst met dictionarys met metabolieten en z_scores
    :return: lijst_patient_id, lijst_metaboliet_id & lijst_z_score
    """
    count = -1
    lijst_patient_id = []
    lijst_metaboliet_id = []
    lijst_z_scores = []
    # Loop door de lijst hierbij is i iedere keer een nieuwe dictionary
    for i in patient_waardes:
        count += 1
        # Per dictionary worden de keys (metabolieten) en values (z_scores) opgehaald
        for key2, value2 in i.items():
            # In de dictionary worden de keys (metabolieten) en values (ids) opgehaald
            for key, value in id_metabolieten.items():
                # Kijkt of de metaboliet uit de patient_waardes hetzelfde is als in de
                # id_metabolieten dictionary
                if key == key2:
                    lijst_patient_id.append(count)
                    lijst_metaboliet_id.append(value)
                    lijst_z_scores.append(value2)

    return lijst_patient_id, lijst_metaboliet_id, lijst_z_scores


def meta_and_dis_no_duplicate(patient_dis_meta):
    """Hier wordt een lijst aangemaakt met geen duplicaties van metabolieten om
    redunantie te voorkomen.

    :param patient_dis_meta: lijst met dictionarys met metabolieten en ziektes
    :return: meta_dis_geen_duplicatie
    """
    print("===", patient_dis_meta)
    meta_dis_geen_duplicatie = {}
    # loopt door de lijst met dictionaries
    for i in patient_dis_meta:
        # loopt door de dictionarys heen met metaboliet als key en ziekte als value
        for key, value in i.items():
            # Kijkt of een metaboliet nog niet voorkomt in een lijst
            if key not in meta_dis_geen_duplicatie:
                meta_dis_geen_duplicatie[key] = value

    return meta_dis_geen_duplicatie


def get_values(connection, cursor, meta_dis_geen_duplicatie, meta_ziekte_ids):
    """

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param meta_dis_geen_duplicatie: dictionary met metabolieten en ziektes zonder duplicaties
    :param meta_ziekte_ids: dictionary met metabolieten met bijbehorende ids
    """
    # Loopt tegelijk door twee dictionarys heen
    for (key, value), (key2, value2) in zip(meta_dis_geen_duplicatie.items(), meta_ziekte_ids.items()):
        # Er wordt gekeken of ; in values staat
        if ";" in value:
            # Er wordt gesplits op ; en de eerste waarde wordt weggehaald aangezien dit
            # een spatie is geworden afgesloten met een komma
            value = value.split(";")[1:]
            # loopt tegelijk door de values van allebei de dictionaries heen
            for i, i2 in zip(value, value2):
                # Gaat naar een volgende functie toe
                compare(connection, cursor, key, i, key2, i2)
        else:
            continue


def compare(connection, cursor, metaboliet, ziekte, metaboliet_id, ziekte_id):
    """

    :param connection: connectie met de database
    :param cursor: Voert een query uit
    :param metaboliet: een metaboliet als string
    :param ziekte: een ziekte als string
    :param metabolite_id: een metaboliet id als int
    :param disease_id: een ziekte id als int
    :return:
    """
    Entrez.email = "A.N.Other@example.com"
    # Zoekwoorden om te textminen
    zoekwoorden = metaboliet + " AND " + ziekte
    # Query om te textminen
    finds = Entrez.esearch(db="pubmed",
                           term=zoekwoorden)
    resultaat = Entrez.read(finds)
    # Nieuwe functie om de database te vullen
    metabolite_disease_table(connection, cursor, metaboliet_id, ziekte_id, resultaat["IdList"])




def main():
    document = "Output_untargeted_metabolomics.xlsx"
    inputgroter = float(1.2)
    inputkleiner = float(1.2)
    inputaantalresultaten = 20

    negativetop, negativetopmetabolite, positivetop, positivetopmetabolite, \
    listpatient, negativedisease, positivedisease = \
    leesbestand(document, inputkleiner, inputgroter, inputaantalresultaten)

    patient_list, patient_z_and_meta = \
    dictionary(negativetop, negativetopmetabolite,
    positivetop, positivetopmetabolite, listpatient)

    patient_list, patient_dis_and_meta = \
    dictionary(negativedisease, negativetopmetabolite,
    positivedisease, positivetopmetabolite, listpatient)

    noduplicate_disease, id_diseases = disease_no_duplicates(patient_dis_and_meta)

    noduplicate_metabolite, id_metabolites = metabolite_no_duplicates(patient_dis_and_meta)

    splitted_dictionary = metaboliet_ziekte_geen_duplicatie(patient_dis_and_meta)

    meta_disease_ids = ids_koppelen(splitted_dictionary, id_diseases, id_metabolites)

    list_patient_id, list_metabolite_id, list_z_id = filter_z_patient_meta(id_metabolites, patient_z_and_meta)

    meta_and_dis_no_duplicates = meta_and_dis_no_duplicate(patient_dis_and_meta)

    "ALLEEN NIET MEER ERUITGECOMMEND HEBBEN ALS JE DE DATABASE WILT VULLEN!"
    # connection, cursor = connection_database()
    # patient_table(connection, cursor, patient_list)
    # disease_table(connection, cursor, noduplicate_disease)
    # metabolite_table(connection, cursor, noduplicate_metabolite)
    # patient_metabolite_table(connection, cursor, list_patient_id, list_metabolite_id, list_z_id)
    # get_values(connection, cursor, meta_and_dis_no_duplicates, meta_disease_ids)


main()
