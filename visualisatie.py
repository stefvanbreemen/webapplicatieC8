# libraries
from pyvis.network import Network
import psycopg2


def connectie():
    """ Deze functie maakt connectie met de database en returned de
    connectie
    :return: connectie met de database
    """
    # Maakt connectie met de database
    connectie = psycopg2.connect(
        host="biocentre.nl",
        database="bio_jaar_2_pg_7",
        user="BI2_PG7",
        password="Blaat1234",
        port="5900")

    cursor = connectie.cursor()

    # Returnen van connectie en cursor
    return connectie, cursor


def data_from_database(connectie, cursor):
    """ Deze functie haalt data uit de database en zet de data om van
    tuples naar 2D list
    :param connectie: connectie hebben met database
    :param cursor: connectie met databas
    :return: 2D - list - metabolieten gekoppelt aan ziektes
    """
    # Data ophalen uit de database
    cursor.execute(("select metabolite, disease "
                    "from diseases full join metabolite_diseases md "
                    "on diseases.id = md.diseases_id "
                    "right join metabolite m on m.id = md.metabolite_id"
                    " where disease not like 'None'"))

    # Alle opgehaalde data opslaan
    metabolen_diseases = cursor.fetchall()

    # Twee lege lijsten aanmaken
    metabolen_diseaseslist = []
    templist = []

    # Door opgehaalde data lopen en toevoegen aan de lijsten
    for x, y in metabolen_diseases:
        templist.append(x)
        templist.append(y)
        metabolen_diseaseslist.append(templist)
        templist = []

    # Returnen van metabolen_diseaseslist
    return metabolen_diseaseslist


def graph(metabolieten_diseaseslist):
    """ Deze functie maakt een interactieve graph met metabolieten en
    bijbehorende ziektes
    """
    # Maakt het netwerk aan
    net = Network()

    # Loopt door de lijst heen en haalt de juiste gegevens eruit
    for i in range(len(metabolieten_diseaseslist)):
        net.add_node(metabolieten_diseaseslist[i][0])
        net.add_node(metabolieten_diseaseslist[i][1])

        net.add_edge(metabolieten_diseaseslist[i][0],
                     metabolieten_diseaseslist[i][1])

    # Zet de breedte en lengte op auto
    net.width = "auto"
    net.height = "auto"

    # Slaat de graph op
    net.save_graph("templates\\my_graph.html")


if __name__ == '__main__':
    conn, cursor = connectie()
    lijst_data = data_from_database(conn, cursor)
    graph(lijst_data)
