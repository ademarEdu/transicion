import sqlite3
import json

conn = sqlite3.connect('molecules.db')

c = conn.cursor()

c.execute("""
    CREATE TABLE molecules (
        name TEXT,
        origin_atom_i_sp3 INTEGER,
        central_atoms_i_sp3 TEXT,
        origin_atom_i_sp2 INTEGER,
        central_atoms_i_sp2 TEXT,
        left_cation_i TEXT,
        right_cation_i TEXT,
        left_bridge_i TEXT,
        right_bridge_i TEXT
        )
    """)

# LI-7
c_atoms_i_sp3 = json.dumps([6, 15])
c_atoms_i_sp2 = json.dumps([6, 15])
l_cation_i = json.dumps([28, 29, 30, 31, 32, 33, 34, 36, 41, 42, 43, 44])
r_cation_i = json.dumps([21, 22, 23, 24, 25, 26, 27, 35, 37, 38, 39, 40])
l_bridge_i = json.dumps([11, 12, 13, 14, 19, 20])
r_bridge_i = json.dumps([7, 8, 9, 10, 17, 18])
c.execute("""INSERT INTO molecules VALUES ('LI-7', 5, ?, 5, ?, ?, ?, ?, ?)""", (c_atoms_i_sp3, c_atoms_i_sp2, l_cation_i, r_cation_i, l_bridge_i, r_bridge_i))

# LI-7-B-D-F
c_atoms_i_sp3 = json.dumps([6, 15])
c_atoms_i_sp2 = json.dumps([6, 15])
l_cation_i = json.dumps([21, 22, 23, 24, 25, 26, 27, 28, 29, 39, 41, 51, 52, 53, 54, 55, 56, 57, 58])
r_cation_i = json.dumps([30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50])
l_bridge_i = json.dumps([11, 12, 13, 14, 19, 20])
r_bridge_i = json.dumps([7, 8, 9, 10, 17, 18])
c.execute("""INSERT INTO molecules VALUES ('LI-7-B-D-F', 5, ?, 5, ?, ?, ?, ?, ?)""", (c_atoms_i_sp3, c_atoms_i_sp2, l_cation_i, r_cation_i, l_bridge_i, r_bridge_i))

# LI-7-B-D-F1
c_atoms_i_sp3 = json.dumps([6, 15])
c_atoms_i_sp2 = json.dumps([6, 15])
l_cation_i = json.dumps([21, 22, 23, 24, 25, 26, 27, 28, 29, 39, 41, 51, 52, 53, 54, 55, 56, 57, 58])
r_cation_i = json.dumps([30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50])
l_bridge_i = json.dumps([11, 12, 13, 14, 19, 20])
r_bridge_i = json.dumps([7, 8, 9, 10, 17, 18])
c.execute("""INSERT INTO molecules VALUES ('LI-7-B-D-F1', 5, ?, 5, ?, ?, ?, ?, ?)""", (c_atoms_i_sp3, c_atoms_i_sp2, l_cation_i, r_cation_i, l_bridge_i, r_bridge_i))

c.execute("SELECT * FROM molecules")
rows = c.fetchall()
for row in rows:
    print(f"{row}\n")

conn.commit()
conn.close()