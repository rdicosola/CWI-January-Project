'''
GUI for VESinverse
Created Winter 2021

@author: Rebecca DiCosola
'''

from tkinter import *

class VESinverseInputs:
    def __init__(self, window):
        window.title("VES Inverse Data Input")
        # print("ndat")
        print("First Name: %s\nLast Name: %s" )


if __name__ == '__main__':
    root = Tk()
    app = VESinverseInputs(root)
    root.mainloop()
