import Tkinter as tk
import tkFileDialog
import tkMessageBox
import os

from tkValidatingEntry import MaxLengthEntry


class MainPage(tk.Frame):
	"""
	main frame of program
	"""

	def __init__(self, parent, controller):
		tk.Frame.__init__(self, parent)
		self.controller = controller
		self.bind_class("Text", "<Control-a>", self.select_all)

		Step1Label = tk.Label(self, text="PAIRWISE SEQUENCE ALIGNMENT", font=('System', 12, 'bold'), background="#202b3d", foreground="#ffffff")
		Step1Label.place(x=450, y=15, anchor="center")

		Step2Label = tk.Label(self, text="Scoring using BLOSUM62 with linear gap score", background="#202b3d", foreground="#ffffff")
		Step2Label.place(x=450, y=35, anchor="center")

		S1Label = tk.Label(self, text="Enter 1st protein sequence:", background="#202b3d", foreground="#ffffff")
		S1Label.place(x=20, y=40)

		S2Label = tk.Label(self, text="Enter 2nd protein sequence:", background="#202b3d", foreground="#ffffff")
		S2Label.place(x=20, y=230)

		S1TextField = tk.Text(self, height=8, width=107, background="#2e3d55", foreground="#ffffff")
		S1TextField.place(x=20, y=60)
		S1TextField.bind("<Tab>", self.focus_next_window)

		S1btn = tk.Button(self, text="Browse", command=lambda: self.upload(S1TextField),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff")
		S1btn.place(x=130, y=200)
		S1btn.bind('<Return>', lambda e: self.upload(S1TextField))
		S1btn.bind("<Tab>", self.focus_next_window)

		S2TextField = tk.Text(self, height=8, width=107, background="#2e3d55", foreground="#ffffff")
		S2TextField.place(x=20, y=250)
		S2TextField.bind("<Tab>", self.focus_next_window)

		S2btn = tk.Button(self, text="Browse", command=lambda: self.upload(S2TextField),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff")
		S2btn.place(x=130, y=390)
		S2btn.bind('<Return>', lambda e: self.upload(S2TextField))
		S2btn.bind("<Tab>", self.focus_next_window)

		UploadS1Label = tk.Label(self, text="or Upload a text file", background="#202b3d", foreground="#ffffff")
		UploadS1Label.place(x=20, y=202)

		UploadS2Label = tk.Label(self, text="or Upload a text file", background="#202b3d", foreground="#ffffff")
		UploadS2Label.place(x=20, y=392)

		# gap option
		gapLabel = tk.Label(self, text="Set gap score:", background="#202b3d", foreground="#ffffff")
		gapLabel.place(x=20, y=425)

		gapScore = MaxLengthEntry(self, "", 4, width=10, background="#2e3d55", foreground="#ffffff")
		gapScore.place(x=100, y=427)
		gapScore.bind("<Tab>", self.focus_next_window)

		# SubmitBtn
		SubmitBtn = tk.Button(self, text="Submit", font=('', 9,'bold'), command=lambda: self.retrieveInputs(S1TextField, S2TextField, gapScore),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff")
		SubmitBtn.place(x=20, y=460)
		SubmitBtn.bind('<Return>', lambda e: self.retrieveInputs(S1TextField, S2TextField, gapScore))

	def upload(self, textfield):
		"""
		This method open text file specified by user in current directory,
		and insert everything within into the textbox.
		"""
		cwd = os.getcwd()
		filename = tkFileDialog.askopenfilename(initialdir=cwd, title="Select file",
												filetypes=[("text files", "*.txt")])
		if filename != "":
			f = open(filename, 'r')
			textfield.delete('1.0', 'end-1c')
			textfield.insert(1.0, f.read())
			f.close()

	def retrieveInputs(self, S1textfield, S2textfield, gapScore):
		"""
		Utility method to retrieve and validate inputs from textbox.
		Display ResultPage on success, return on failure
		"""

		# retrieve input, remove whitespace and newline characters
		seq1 = S1textfield.get("1.0", 'end-1c')
		seq1 = ''.join(seq1.split())

		seq2 = S2textfield.get("1.0", 'end-1c')
		seq2 = ''.join(seq2.split())

		gap = gapScore.get()

		# validate inputs, and call showResults on success
		if not seq1 or not seq2 or not gap:
			tkMessageBox.showinfo("Error!", "Please enter the sequences to align, as well as the gap score.")
			return

		try:
			gap = int(gap)
		except ValueError:
			tkMessageBox.showinfo("Error!", "Please enter only integers into gap score.")
			return

		if not seq1.isalpha() or not seq2.isalpha():
			tkMessageBox.showinfo("Error!", "Please enter only alphabets into both sequences.")
			return

		if not self.run(seq1, seq2, gap):
			tkMessageBox.showinfo("Error!", "An unidentified sequence character has caused an error in sequence "
											"alignment. Please verify and try again.")
			return
		else:
			self.controller.frames["ResultPage"].display()
			self.controller.show_frame("ResultPage")
			S1textfield.delete('1.0', 'end-1c')
			S2textfield.delete('1.0', 'end-1c')
			gapScore.delete(0, 'end')

	def run(self, seq1, seq2, gap):
		"""
		This method updates the needleman object in main and invokes methods of Needleman class,
		to print the results to console. Return True on success and False on failure
		"""

		self.controller.needleman.gapScore = gap
		self.controller.needleman.seq1 = seq1
		self.controller.needleman.seq2 = seq2
		if self.controller.needleman.alignSequences():
			return True
		else:
			return False

	def select_all(self, e):
		"""Control-a keyboard shortcut"""
		e.widget.tag_add("sel", "1.0", "end")

	def focus_next_window(self, e):
		"""tab keyboard shortcut"""
		e.widget.tk_focusNext().focus()
		return "break"
