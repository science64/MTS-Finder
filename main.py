from threading import Thread
from tkinter import *
from tkinter.font import Font
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog, messagebox
from functions import *
from tkinter import END,Text,Tk,Frame,Entry,Label,Button,StringVar

class MyWindow():
    def __init__(self, parent):

        self.frame = Frame(parent, width=840, height=560)
        self.font = Font(family="Times New Roman", size=16)
        self.request_timeout = 30
        global root
        root = parent
        # Label(self.frame, text="Select a PSMs file!", fg='#ffe7b6', bg='#b97455', font=self.font).place(x=100, y=15)

        self.fontRadio = Font(family="Times New Roman", size=13)
        self.var = StringVar()

        ####### Browse label and button ######
        self.browseLabel = Label(self.frame, font=Font(family="Times New Roman", size=12),
                                 text="Please, choose a Peptides:")
        self.browseLabel.place(x=80, y=30)
        self.browseButton = Button(self.frame, text="Browse", justify=LEFT,
                                   font=Font(family="Times New Roman", size=12, weight='bold'), command=self.browse)
        self.browseButton.place(x=255, y=25)

        ######## Normalization radio Button #######
        self.varDecision = StringVar()
        self.norm = Radiobutton(root, font=self.fontRadio, text="Normalized values", value=True,
                               variable=self.varDecision)
        self.norm.select()
        self.norm.place(x=400, y=30)

        self.nonNorm = Radiobutton(root, font=self.fontRadio, text="Non-Normalized values", value=False,
                                variable=self.varDecision)
        self.nonNorm.place(x=580, y=30)

        Label(self.frame, text="Output Name:", font=Font(family="Times New Roman", size=12)).place(x=80, y=80)
        self.outputNamebox = ScrolledText(self.frame, font=Font(family="Times New Roman", size=12), bd=2)
        self.outputNamebox.place(x=180, y=75, width=320, height=40)

        Label(self.frame, text="Conditions:", font=Font(family="Times New Roman", size=12)).place(x=100, y=130)
        self.conditionbox = ScrolledText(self.frame, font=Font(family="Times New Roman", size=12), bd=2)
        self.conditionbox.place(x=180, y=125, width=500, height=60)

        condtionsFromText = open('condtions.txt').read()
        self.conditionbox.insert(END, condtionsFromText)

        self.browseButtonCondition = Button(self.frame, text="Browse", justify=LEFT,
                                            font=Font(family="Times New Roman", size=12, weight='bold'),
                                            command=self.browse_condition)
        self.browseButtonCondition.place(x=700, y=132)

        Label(self.frame, text="Pairs:", font=Font(family="Times New Roman", size=12)).place(x=130, y=210)
        self.pairsbox = ScrolledText(self.frame, font=Font(family="Times New Roman", size=12), bd=2)
        self.pairsbox.place(x=180, y=200, width=320, height=60)

        pairsFromText = open('pairs.txt').read()
        self.pairsbox.insert(END, pairsFromText)

        self.browseButtonPairs = Button(self.frame, text="Browse", justify=LEFT,
                                        font=Font(family="Times New Roman", size=12, weight='bold'),
                                        command=self.browse_pairs)
        self.browseButtonPairs.place(x=520, y=210)

        self.statusbar = ScrolledText(self.frame, state='disabled')
        self.statusbar.place(x=100, y=280, width=650, height=180)

        self.runbutton = Button(self.frame, text='RUN', fg='black', bg='#b4e67e',
                                font=Font(family="Times New Roman", size=18, weight='bold'),
                                command=self.runbutton_click)
        self.runbutton.place(x=250, y=480, width=150, height=50)

        self.openbutton = Button(self.frame, text='Open', fg='black', bg='#FF5733',
                                 font=Font(family="Times New Roman", size=18, weight='bold'), command=self.open_click)
        self.openbutton.place(x=450, y=480, width=150, height=50)

        self.openbutton.configure(state='disabled')

        self.frame.pack()

        self.update_status_box('\n\t >> :: MTS Finder Bot Started! :: <<\n')
        self.update_status_box('\n------------------------------------------------------------------------------\n')

    def Message(self, title, message):
        messagebox.showinfo(title=title, message=message)

    def update_status_box(self, text):
        self.statusbar.configure(state='normal')
        self.statusbar.insert(END, text)
        self.statusbar.see(END)
        self.statusbar.configure(state='disabled')

    def clear_status_box(self):
        self.statusbar.configure(state='normal')
        self.statusbar.delete(1.0, END)
        self.statusbar.see(END)
        self.statusbar.configure(state='disabled')

    def check_main_thread(self):
        root.update()
        if self.myThread.is_alive():
            root.after(1000, self.check_main_thread)
        else:
            self.x = True

    def open_click(self):
        import os
        os.startfile(f'{self.outputLocationPath}/{self.outputLocation.strip()}.xlsx')

    def browse(self):

        self.filename = filedialog.askopenfile(parent=self.frame, mode='rb', title='Choose a peptide file')
        self.filenamePretify = str(self.filename).split('/')[-1].split("'>")[0]
        if self.filenamePretify == "None":
            self.Message('Error!', 'Please choose a file!')
            return 0
        self.update_status_box(f'\n "{self.filenamePretify}" file is chosen! \n')

        self.outputLocationPath =  str(self.filename).split("'")[1].replace(str(self.filename).split("'")[1].split("/")[-1],'')

        # self.outputLocationPretify = str(self.outputLocation).split('/')[-1].split("'>")[0]
    def browse_condition(self):
        self.filename_condition = filedialog.askopenfile(parent=self.frame, mode='rb', title='Please, choose a condition text file.')
        self.filenamePretify_condition = str(self.filename_condition).split('/')[-1].split("'>")[0]
        if self.filenamePretify_condition == "None":
            self.Message('Error!', 'Please choose a file!')
            return 0
        self.outputLocationPath_condition =  str(self.filename_condition).split("'")[1].replace(str(self.filename_condition).split("'")[1].split("/")[-1],'')
        condtionsFromText = str(open(self.outputLocationPath_condition+self.filenamePretify_condition).read()).strip()
        self.conditionbox.delete('1.0', END)
        self.conditionbox.insert(END, condtionsFromText)

    def browse_pairs(self):
        self.filename_pairs = filedialog.askopenfile(parent=self.frame, mode='rb', title='Please, choose a pairs text file.')
        self.filenamePretify_pairs = str(self.filename_pairs).split('/')[-1].split("'>")[0]
        if self.filenamePretify_pairs == "None":
            self.Message('Error!', 'Please choose a file!')
            return 0
        self.outputLocationPath_pairs =  str(self.filename_pairs).split("'")[1].replace(str(self.filename_pairs).split("'")[1].split("/")[-1],'')
        pairsFromText = str(open(self.outputLocationPath_pairs+self.filenamePretify_pairs).read()).strip()
        self.pairsbox.delete('1.0', END)
        self.pairsbox.insert(END, pairsFromText)

    def runbutton_click(self):
        #self.update_status_box(f'\n "{self.varDecision.get()}" is selected! \n')
        self.runbutton.configure(state='disabled')
        self.myThread = Thread(target=self.engine)
        self.myThread.daemon = True
        self.myThread.start()
        root.after(1000, self.check_main_thread)

    def engine(self):
        try:
            self.update_status_box(f'\n The file is reading..! \n')
            if '.xlsx' in self.filenamePretify:
                self.fileRead = pd.read_excel(self.outputLocationPath+self.filenamePretify)
            elif '.txt' in self.filenamePretify:
                self.fileRead = pd.read_csv(self.outputLocationPath+self.filenamePretify,sep='\t',header=0)
        except Exception as e:
            self.update_status_box(f'\n Error is "{e}"! \n')
            self.Message('Error!', 'An Error Occured, please choose a file before run!')
            return 0

        self.update_status_box(f'\n The file is read! \n')

        # conditions = ['Light', '0DMSO', '0DMSO', '0DMSO', 'Rotenone', 'Rotenone', 'Rotenone', 'Antimycin', 'Antimycin', 'Antimycin', 'Boost']
        self.conditions = self.conditionbox.get("1.0", END)

        condtionsFile = open(f'{self.outputLocationPath}/condtions.txt', 'w')
        condtionsFile.write(self.conditions)
        condtionsFile.close()

        conditionsFinal = self.conditions.split(',')
        conditionsFinal[-1] = conditionsFinal[-1].strip()

        # pairs = [['0DMSO', 'Rotenone'], ['0DMSO', 'Antimycin']]
        self.pairs = self.pairsbox.get("1.0", END)

        pairsFile = open(f'{self.outputLocationPath}/pairs.txt', 'w')
        pairsFile.write(self.pairs)
        pairsFile.close()

        pairsFinal = self.pairs.split(';')
        pairsFinal[-1] = pairsFinal[-1].strip()
        pairsFinal = [i.strip() for i in pairsFinal]
        pairsFinal = [i.lstrip() for i in pairsFinal]
        # [print(x) for x in thislist]
        pairsFinal = [pairs.split('/') for pairs in pairsFinal]

        self.update_status_box(f'\n Conditions: {self.conditions.strip()} \n')

        self.update_status_box(f'\n Pairs: {self.pairs.strip()} \n')

        self.update_status_box(f'\n Running..! \n')

        normalization = self.varDecision.get()

        self.data = mtsFinderEngine(self.fileRead, conditionsFinal, pairsFinal, normalization)
        # self.data['Accession'] = self.data.index
        # self.data['Gene Symbol'] = ''

        self.update_status_box(f'\n Completed..! \n')
        self.update_status_box(f'\n Data is saving..! \n')

        self.outputLocation = self.outputNamebox.get("1.0", END)

        try:
            self.data.to_excel(f'{self.outputLocationPath}/{self.outputLocation.strip()}.xlsx', index=False, engine="openpyxl")
            # self.datamePROD = pd.read_excel(f'{self.outputLocationPath}/{self.outputLocation.strip()}.xlsx')
            # self.data = GeneNameEngine(self.datamePROD)
            # self.data = mito_human(self.data)
            # self.data.to_excel(f'{self.outputLocationPath}/{self.outputLocation.strip()}.xlsx', index=False,
            #                    engine="openpyxl")
            print(f'{self.outputLocation.strip()}!')
            self.update_status_box(f'\n Saved as {self.outputLocation.strip()}! \n')
            self.Message('Finished!', 'Application Completed!')
            self.openbutton.configure(state='normal')
            self.runbutton.configure(state='normal')
        except Exception as e:
            self.update_status_box(f'\n Error is "{e}"! \n')
            self.Message('Error!', 'An Error Occured, please fix it and rerun!')

if __name__ == '__main__':

    root = Tk()

    root.title("MTS Finder App v3.5 by S. Bozkurt @2023")
    root.geometry("840x560")
    root.resizable(0, 0)

    # get screen width and height
    ws = root.winfo_screenwidth()  # width of the screen
    hs = root.winfo_screenheight()  # height of the screen

    # calculate x and y coordinates for the Tk root window
    x = (ws / 2) - (840 / 2)
    y = (hs / 2) - (560 / 2)

    # set the dimensions of the screen
    # and where it is placed
    root.geometry('%dx%d+%d+%d' % (840, 560, x, y))

    # root.iconphoto(False, tkinter.PhotoImage(file='icon.ico'))
    MyWindow(root)
    root.wm_iconbitmap('files//icon.ico')
    root.mainloop()