# Author of this program -- Narendra N. Das (JPL)

#This python program is to perform ingestion 
# of Cultivar for rice or maize

#Import libraries
from tkinter import *
#from tkinter import ttk
# import filedialog module 
from tkinter import filedialog 
# import utils from rheas.dssat
from rheas.dssat import utils as dutils

#####################################################################
window = Tk()
#Set window background color 
window.config(background = "yellow")
#Set window title
window.title("RHEAS Cultivar Tool")
#Set window size 
window.geometry('900x800')

lbl_PD = Label(window, text="Insert Cultivar into the RHEAS Database", font=('Helvetica', 12, 'bold'), background='yellow') 
lbl_PD.grid(column=0, row=0)

#####################################################################

def only_numbers(char):
    #return char.isdigit()
    return char.isdecimal()

validation = window.register(only_numbers)

#####################################################################
#Select the crop type
# Create a Tkinter variable
tkvar_crp = StringVar(window)
# Dictionary with options
choices = { 'maize', 'rice' }
tkvar_crp.set('maize') # set the default option
crpname = tkvar_crp.get()

popupMenu = OptionMenu(window, tkvar_crp, *choices)
Label(window, text="Select Crop Type", background='yellow').grid(row = 9, column = 0)
popupMenu.grid(row = 9, column =1)

def change_dropdown(*args): # on change dropdown value
    #print( tkvar_crp.get() )
    crpname = tkvar_crp.get()
    #print(crpname)

tkvar_crp.trace('w', change_dropdown) # link function to change dropdown

#####################################################################

label_blank1 = Label(window,  text = " ", width = 20, height = 1,  fg = "white", background='yellow')
label_blank1.grid(column = 0, row = 12)  

#####################################################################

p1_m = DoubleVar()
p2_m = DoubleVar()
p5_m = DoubleVar()
g2_m = DoubleVar()
g3_m = DoubleVar()
phint_m = DoubleVar()

p1_r = DoubleVar()
p2o_r = DoubleVar()
p2r_r = DoubleVar()
p5_r = DoubleVar()
g1_r = DoubleVar()
g2_r = DoubleVar()
g3_r = DoubleVar()
g4_r = DoubleVar()

def EnableCultivar():
    crpname = tkvar_crp.get()
    #print(crpname)
    if (crpname == 'maize'):
        p1_str_m.configure(state='normal')
        p2_str_m.configure(state='normal')
        p5_str_m.configure(state='normal')
        g2_str_m.configure(state='normal')
        g3_str_m.configure(state='normal')
        phint_str_m.configure(state='normal')
        Culname_str_m.configure(state='normal')
    elif (crpname == 'rice'):
        p1_str_r.configure(state='normal')
        p2r_str_r.configure(state='normal')
        p5_str_r.configure(state='normal')
        p2o_str_r.configure(state='normal')
        g1_str_r.configure(state='normal')
        g2_str_r.configure(state='normal')
        g3_str_r.configure(state='normal')
        g4_str_r.configure(state='normal')
        Culname_str_r.configure(state='normal')
    else:
        print('No crop selected.')

def DisableCultivar():
    crpname = tkvar_crp.get()
    #print(crpname)
    if (crpname == 'maize'):
        p1_str_m.configure(state='disabled')
        p2_str_m.configure(state='disabled')
        p5_str_m.configure(state='disabled')
        g2_str_m.configure(state='disabled')
        g3_str_m.configure(state='disabled')
        phint_str_m.configure(state='disabled')
        Culname_str_m.configure(state='disabled')
    elif (crpname == 'rice'):
        p1_str_r.configure(state='disabled')
        p2r_str_r.configure(state='disabled')
        p5_str_r.configure(state='disabled')
        p2o_str_r.configure(state='disabled')
        g1_str_r.configure(state='disabled')
        g2_str_r.configure(state='disabled')
        g3_str_r.configure(state='disabled')
        g4_str_r.configure(state='disabled')
        Culname_str_r.configure(state='disabled')
    else:
        print('No crop selected.')

Label(window, text="Maize Cultivar", background='yellow').grid(row=13, column=0)
Label(window, text="p1", background='yellow').grid(row=14, column=0)
Label(window, text="p2", background='yellow').grid(row=15, column=0)
Label(window, text="p5", background='yellow').grid(row=16, column=0)
Label(window, text="g2", background='yellow').grid(row=17, column=0)
Label(window, text="g3", background='yellow').grid(row=18, column=0)
Label(window, text="phint", background='yellow').grid(row=19, column=0)
Label(window, text="Maize Cultivar Name", background='yellow').grid(row=20, column=0)
p1_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p1_m)
p1_str_m.grid(column=1, row=14)
p1_str_m.focus()
p1_str_m.insert(0, '0')
p2_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p2_m)
p2_str_m.grid(column=1, row=15)
p2_str_m.focus()
p2_str_m.insert(0, '0')
p5_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p5_m)
p5_str_m.grid(column=1, row=16)
p5_str_m.focus()
p5_str_m.insert(0, '0')
g2_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g2_m)
g2_str_m.grid(column=1, row=17)
g2_str_m.focus()
g2_str_m.insert(0, '0')
g3_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g3_m)
g3_str_m.grid(column=1, row=18)
g3_str_m.focus()
g3_str_m.insert(0, '0')
phint_str_m = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = phint_m)
phint_str_m.grid(column=1, row=19)
phint_str_m.focus()
phint_str_m.insert(0, '0')
Culname_str_m = Entry(window, width=10, state='disabled')
Culname_str_m.grid(column=1, row=20)
Culname_str_m.focus()
Culname_str_m.insert(0, '0')


Label(window, text="Rice Cultivar", background='yellow').grid(row=13, column=2)
Label(window, text="p1", background='yellow').grid(row=14, column=2)
Label(window, text="p2r", background='yellow').grid(row=15, column=2)
Label(window, text="p5", background='yellow').grid(row=16, column=2)
Label(window, text="p2o", background='yellow').grid(row=17, column=2)
Label(window, text="g1", background='yellow').grid(row=18, column=2)
Label(window, text="g2", background='yellow').grid(row=19, column=2)
Label(window, text="g3", background='yellow').grid(row=20, column=2)
Label(window, text="g4", background='yellow').grid(row=21, column=2)
Label(window, text="Rice Cultivar Name", background='yellow').grid(row=22, column=2)
p1_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p1_r)
p1_str_r.grid(column=3, row=14)
p1_str_r.focus()
p1_str_r.insert(0, '0')
p2r_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p2r_r)
p2r_str_r.grid(column=3, row=15)
p2r_str_r.focus()
p2r_str_r.insert(0, '0')
p5_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p5_r)
p5_str_r.grid(column=3, row=16)
p5_str_r.focus()
p5_str_r.insert(0, '0')
p2o_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = p2o_r)
p2o_str_r.grid(column=3, row=17)
p2o_str_r.focus()
p2o_str_r.insert(0, '0')
g1_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g1_r)
g1_str_r.grid(column=3, row=18)
g1_str_r.focus()
g1_str_r.insert(0, '0')
g2_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g2_r)
g2_str_r.grid(column=3, row=19)
g2_str_r.focus()
g2_str_r.insert(0, '0')
g3_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g3_r)
g3_str_r.grid(column=3, row=20)
g3_str_r.focus()
g3_str_r.insert(0, '0')
g4_str_r = Entry(window, width=10, validate="key", validatecommand=(validation, '%S'), state='disabled', textvariable = g4_r)
g4_str_r.grid(column=3, row=21)
g4_str_r.focus()
g4_str_r.insert(0, '0')
Culname_str_r = Entry(window, width=10, state='disabled')
Culname_str_r.grid(column=3, row=22)
Culname_str_r.focus()
Culname_str_r.insert(0, '0')

#####################################################################

label_blank5 = Label(window,  text = " ", width = 20, height = 1,  fg = "white", background='yellow')
label_blank5.grid(column = 0, row = 24)  

#####################################################################

button_Establish_cultivar = Button(window,  text = "Enable Cultivar Entry", command = EnableCultivar)
button_Establish_cultivar.grid(column = 0,row = 25) 

button_Establish_cultivar = Button(window,  text = "Disable Cultivar Entry", command = DisableCultivar)
button_Establish_cultivar.grid(column = 1,row = 25) 


#####################################################################

label_blank2 = Label(window,  text = " ", width = 10, height = 1,  fg = "white", background='yellow')
label_blank2.grid(column = 0, row = 31)  

#####################################################################
def browseFiles(): 
    Shpfilename = filedialog.askopenfilename(initialdir = "/home/parallels/RHEAS/", title = "Select a File", 
                                          filetypes = (("Text files", "*.shp*"), ("all files", "*.*"))) 
    # Change label contents 
    #label_file_explorer.configure(text="File Opened: "+filename) 
    label_file_explorer.configure(text=Shpfilename)
    #print(Shpfilename)

label_file_explorer = Label(window,  text = "File Explorer to select Shapefile", width = 20, height = 2,  fg = "blue", background='yellow') 

button_explore = Button(window,  text = "Browse Files", command = browseFiles)  
   
label_file_explorer.grid(column = 0, row = 32) 
 
button_explore.grid(column = 0, row = 33) 
   
#####################################################################

label_blank3 = Label(window,  text = " ", width = 10, height = 1,  fg = "white", background='yellow')
label_blank3.grid(column = 0, row = 34)  

#####################################################################

Label(window, text="Database Name", background='yellow').grid(row=35, column=0)
db_str = Entry(window, width=10)
db_str.grid(column=1, row=35)
db_str.focus()
db_str.insert(0, 'rheas')

#####################################################################

label_blank3 = Label(window,  text = " ", width = 10, height = 1,  fg = "white", background='yellow')
label_blank3.grid(column = 0, row = 36)  
#####################################################################

ensem_no = IntVar()

Label(window, text="Ensemble Number", background='yellow').grid(row=37, column=0)
ensem_str = Entry(window, width=10, textvariable = ensem_no)
ensem_str.grid(column=1, row=37)
ensem_str.focus()
ensem_str.insert(0, '1')

#####################################################################

label_blank3 = Label(window,  text = " ", width = 10, height = 1,  fg = "white", background='yellow')
label_blank3.grid(column = 0, row = 38)  
#####################################################################


def create_rice_cultivar(p1_r, p2r_r, p5_r, p2o_r, g1_r, g2_r, g3_r, g4_r):
    return {'p1': p1_r, 'p2r': p2r_r, 'p5': p5_r, 'p2o': p2o_r, 'g1': g1_r, 'g2': g2_r, 'g3': g3_r, 'g4': g4_r}

def create_maize_cultivar(p1_m, p2_m, p5_m, g2_m, g3_m, phint_m):
    return {'p1': p1_m, 'p2': p2_m, 'p5': p5_m, 'g2': g2_m, 'g3': g3_m, 'phint': phint_m}

def PopulateCultivar():
    crpname = tkvar_crp.get()
    db_name = db_str.get()
    ensem_no = int(ensem_str.get())
    if crpname=='maize':
        p1_m = float(p1_str_m.get())
        p2_m = float(p2_str_m.get())
        p5_m = float(p5_str_m.get())
        g2_m = float(g2_str_m.get())
        g3_m = float(g3_str_m.get())
        phint_m = float(phint_str_m.get())
        maize_name = Culname_str_m.get()
        maize_culti = create_maize_cultivar(p1_m, p2_m, p5_m, g2_m, g3_m, phint_m)
        maize_culti['name'] = "'"+ maize_name + "'"
        params = [maize_culti]
        shp_filename = label_file_explorer.cget('text')
        print(shp_filename)
        print("'" + maize_name + "'")
        print('params')
        print(params)
        print('dbname: ', db_name)
        print('Ensem Number: ', ensem_no)
        #dutils.addCultivar("rheas", shp_filename, params, crop=crpname)
        dutils.addCultivar(db_name, shp_filename, params, ensem_no, crop=crpname)
    elif crpname=='rice':
        p1_r = float(p1_str_r.get())
        p2o_r = float(p2o_str_r.get())
        p2r_r = float(p2r_str_r.get())
        p5_r = float(p5_str_r.get())
        g1_r = float(g1_str_r.get())
        g2_r = float(g2_str_r.get())
        g3_r = float(g3_str_r.get())
        g4_r = float(g4_str_r.get())
        rice_name = Culname_str_r.get()
        rice_culti = create_rice_cultivar(p1_r, p2r_r, p5_r, p2o_r, g1_r, g2_r, g3_r, g4_r)
        rice_culti['name'] = "'" + rice_name + "'"
        params = [rice_culti]
        shp_filename = label_file_explorer.cget('text')
        print(shp_filename)
        print("'" + rice_name + "'")
        print(params)
        print('dbname: ', db_name)
        print('Ensem Number: ', ensem_no)
        #dutils.addCultivar("rheas", shp_filename, params, crop=crpname)
        dutils.addCultivar(db_name, shp_filename, params, ensem_no, crop=crpname)
    else:
        print('No crop selected.')

Insert_Cultivar = Button(window,  text = "Insert_Cultivar", command = PopulateCultivar) 
Insert_Cultivar.grid(column = 0, row = 39) 
#####################################################################

label_blank4 = Label(window,  text = " ", width = 10, height = 1,  fg = "white", background='yellow')
label_blank4.grid(column = 0, row = 41)  

#####################################################################

button_exit = Button(window,  text = "Exit", command = exit) 
button_exit.grid(column = 0,row = 42) 

window.mainloop()



