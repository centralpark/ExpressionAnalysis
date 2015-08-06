import sys; print('%s %s' % (sys.executable or sys.platform, sys.version))
import MySQLdb
import sys; print('%s %s' % (sys.executable or sys.platform, sys.version))
import mysql.connector
cnx = mysql.connector.connect(user='root',database='ROCHE')
cursor = cnx.cursor()
query = 'SELECT name,name2,chrom FROM refGene WHERE name2=%s')
query = 'SELECT name,name2,chrom FROM refGene WHERE name2=%s'
cursor.execute(query, 'CD74')
cursor.execute(query, ('CD74'))
cursor.execute(query, ('CD74',))
for (name,name2,chrom) in cursor:
    print('%s\t%s\t%s\n' % (name,name2,chrom))
    
cursor
cursor.rowcount
cursor.execute(query, ('HSH',))
cursor.rowcount
cursor.execute(query, ('SRF',))
cursor.fetchall()
cursor.execute(query, ('SRF',))
cursor.rowcount
cursor.fetchall()
cursor.rowcount
cursor.arraysize
cursor._rowcount
cursor.execute(query, ('SRF',))
cursor._rowcount
len(cursor)
cursor.fetchall()
cursor.fetchone()
cursor.execute(query, ('HSH',))
cursor.fetchall()
cursor.execute(query, ('HSH',))
result = cursor.fetchall()
len(result)
cursor.execute(query, ('SRF',))
cursor.fetchone()
(name,name2,chr) = cursor.fetchone()
name
name2
a = 5
b = 10
a,b=b,a
a
b
import sys; print('%s %s' % (sys.executable or sys.platform, sys.version))
cnt
import sys; print('%s %s' % (sys.executable or sys.platform, sys.version))
import sys; print('%s %s' % (sys.executable or sys.platform, sys.version))
import os; os.chdir('/Users/HSH/Desktop')
execfile('/Users/HSH/Roche/workspace/GeneFusion/publish/20140612.py')
len(lines)
execfile('/Users/HSH/Roche/workspace/GeneFusion/publish/20140612.py')
cnt.most_common(5)
