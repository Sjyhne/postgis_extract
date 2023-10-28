import psycopg2 as psy
import yaml

from pathlib import Path

def connect():
    conn = None
    
    try:
        config = yaml.load(Path("postgisutils/connect_params.yml").open("r"), yaml.Loader)
        
        conn = psy.connect(**config)
        
        cur = conn.cursor()
        cur.execute("SELECT Version()")
        version = cur.fetchone()
        
        if version != None:
            return cur
                
    except Exception as e:
        print(e)

def close(cur):
    cur.close()
    

def get_cursor():
    return connect()

def close_cursor(cur):
    close(cur)