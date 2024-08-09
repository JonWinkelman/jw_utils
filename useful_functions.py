import random
from time import sleep
import pandas as pd

def move_cursor(mins):
    "move cursor in random dir once every min for given minutes"
    
    import pyautogui
    for i in range(mins):
        pyautogui.moveTo(random.randint(1,500), random.randint(1,500), duration = 1)
        sleep(60)

def save_df_to_doc(df):
    fp = '/Users/jonathanwinkelman/Documents/temp.csv'
    df.to_csv(fp)
    return fp