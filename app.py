import psycopg2
from flask import Flask, render_template, request
from pyvis.network import Network

app = Flask(__name__)


@app.route('/', methods=["GET"])
def website():
    """De homepagina van de website

    :return: HTML pagina met een overzicht van de website
    """
    # Returnen van HTML pagina
    return render_template("home.html")


@app.route('/zoek', methods=["POST", "GET"])
def web_home():
    """De zoekpagina van de website met een connectie aan de database

    :return: HTML pagina met een zoekfunctie op de website en resultaat uit de database
    """
    if request.method == "POST":
        # Maakt connectie met de database
        mydb = psycopg2.connect(
            host="biocentre.nl",
            database="bio_jaar_2_pg_7",
            user="BI2_PG7",
            password="Blaat1234",
            port="5900")

        cursor = mydb.cursor()

        # Haalt p_id uit form op en linkt knop aan de database
        pid = request.form.get("p_id", "")

        # Data ophalen uit de database
        cursor.execute("select * from diseases"
                       " right outer join metabolite_diseases md on diseases.id = md.diseases_id"
                       " right outer join metabolite m on md.metabolite_id = m.id"
                       " right outer join metabolite_patient mp on m.id = mp.metabolite_id"
                       " right outer join patient p on mp.patient_id = p.id"
                       " where p.patient_code like '%" + pid + "%'")

        # Alle opgehaalde data opslaan
        resultaat = cursor.fetchall()

        # Cursor en database closen
        cursor.close()
        mydb.close()

        # Returnen HTML pagina en resultaat
        return render_template("zoek.html",
                               len=len(resultaat), invoer=resultaat, pid=pid)
    else:
        # Returnen HTML pagina
        return render_template("zoek.html")


@app.route('/visualisatie', methods=["GET"])
def web_data():
    """De database van de website

    :return: HTML pagina met inhoud van de database op de website
    """
    # Visualisatie.py openen
    file = open(r'visualisatie.py', 'r').read()
    # Returnen HTML pagina en file runnen op de pagina
    return render_template("database.html", file=exec(file))


if __name__ == '__main__':
    app.run()
