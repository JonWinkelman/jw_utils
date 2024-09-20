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


def shortest_distance(gene1, gene2):
    """
    Calculate the shortest distance between two genes.
    
    Parameters:
    gene1 (tuple): A tuple representing the start and end position of the first gene (start1, end1).
    gene2 (tuple): A tuple representing the start and end position of the second gene (start2, end2).
    
    Returns:
    int: The shortest distance between the two genes. Returns 0 if they overlap.
    """
    start1, end1 = gene1
    start2, end2 = gene2
    
    # Check if the genes overlap
    if end1 >= start2 and end2 >= start1:
        return 0
    else:
        # Calculate the shortest distance
        return min(abs(start2 - end1), abs(start1 - end2))

