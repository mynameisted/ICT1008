import Tkinter as tk
import ttk
import tkMessageBox
import tkFileDialog
import os

class ResultPage(tk.Frame):

	def __init__(self, parent, controller):
		tk.Frame.__init__(self, parent)
		ttkStyle = ttk.Style()
		ttkStyle.theme_use('clam')
		ttkStyle.configure("TSeparator", background="#202b3d")
		self.controller = controller

	def display(self):
		"""
		This method displays all widgets in the resultPage frame
		"""

		heading = tk.Label(self, text="PAIRWISE SEQUENCE ALIGNMENT", font=("System", 12, "bold"), background="#202b3d", foreground="#ffffff")
		heading.place(x=450, y=30, anchor="center")

		subheading = tk.Label(self, text="Scoring using BLOSUM62 with linear gap score", background="#202b3d", foreground="#ffffff")
		subheading.place(x=450, y=50, anchor="center")

		# Optimal Alignment
		label = tk.Label(self, text="Optimal Alignment:", font=('', 9), background="#202b3d", foreground="#ffffff")
		label.place(x=20, y=80)


		alignedSeqFrame = tk.Frame(self, height=80, width=860)
		alignedSeqFrame.place(x=20,y=100)

		seqTextFrame = tk.Frame(alignedSeqFrame, width=580, height=350)
		seqTextFrame.pack(fill="both", expand=True)
		seqTextFrame.grid_propagate(False)
		seqTextFrame.grid_rowconfigure(0, weight=1)
		seqTextFrame.grid_columnconfigure(0, weight=1)

		seqText = tk.Text(seqTextFrame, borderwidth=3, background="#2e3d55", foreground="#ffffff")
		seqText.config(font=("consolas", 13))
		seqText.grid(row=0, column=0, sticky="nsew",ipadx=10,ipady=10)

		arr = self.splitSeq(50)
		for i in range(len(arr[0])-1,-1, -1):
			seqText.insert('1.0',str(arr[0][i]) + '\n')
			seqText.insert('2.0',str(arr[1][i]) + '\n')
			seqText.insert('3.0',str(arr[2][i]) + '\n\n')
		# create a Scrollbar and associate it with txt
		scrollbar = tk.Scrollbar(seqTextFrame, command=seqText.yview)
		scrollbar.grid(row=0, column=1, sticky='nsew')
		seqText['yscrollcommand'] = scrollbar.set
		seqText['state'] = 'disabled'

		# Results
		tk.Label(self, text="Alignment Results with Gap Penalty of " +
							str(self.controller.needleman.gapScore) + ":", font=('', 9), background="#202b3d", foreground="#ffffff").place(x=630, y=80)

		tk.Label(self, text="    Identity: " + str(self.controller.needleman.alignmentResults["Identity"])
							+ "/" + str(self.controller.needleman.alignmentResults["Length"]) + ' (' +
							str(round(self.controller.needleman.alignmentResults["IdentityPercent"],
									  2)) + "%)", background="#202b3d", foreground="#ffffff").place(x=630, y=110)

		tk.Label(self, text="    Similarity: " + str(self.controller.needleman.alignmentResults["Similarity"])
							+ "/" + str(self.controller.needleman.alignmentResults["Length"]) + ' (' +
							str(round(self.controller.needleman.alignmentResults["SimilarityPercent"],
									  2)) + "%)", background="#202b3d", foreground="#ffffff").place(x=630, y=130)

		tk.Label(self, text="    Match Count: " + str(self.controller.needleman.alignmentResults["Match"]) + ' (' +
							str(round(self.controller.needleman.alignmentResults["MatchPercent"], 2))
							+ "%)", background="#202b3d", foreground="#ffffff").place(x=630, y=160)

		tk.Label(self, text="    Mismatch Count: " + str(self.controller.needleman.alignmentResults["Mismatch"]) + ' (' +
							str(round(self.controller.needleman.alignmentResults["MismatchPercent"], 2))
							+ "%)", background="#202b3d", foreground="#ffffff").place(x=630, y=180)

		tk.Label(self, text="    Gap Count: " + str(self.controller.needleman.alignmentResults["Gap"]) + ' (' +
							str(round(self.controller.needleman.alignmentResults["GapPercent"], 2))
							+ "%)", background="#202b3d", foreground="#ffffff").place(x=630, y=200)

		tk.Label(self, text="    Alignment Path Score: " +
							str(self.controller.needleman.alignmentResults["PathScore"]), background="#202b3d", foreground="#ffffff").place(x=630, y=230)

		ttk.Separator(self, orient="horizontal").place(x=635, y=260, height=5, width=160)
		tk.Label(self, text="    Sequence 1 Length: " +
							str(len(self.controller.needleman.seq1)), background="#202b3d", foreground="#ffffff").place(x=630, y=270)

		tk.Label(self, text="    Sequence 2 Length: " +
							str(len(self.controller.needleman.seq2)), background="#202b3d", foreground="#ffffff").place(x=630, y=290)

		tk.Label(self, text="    Aligned Sequence Length: " +
							str(len(self.controller.needleman.tracebackSequences[0])), background="#202b3d", foreground="#ffffff").place(x=630, y=320)

		tk.Label(self, text="Alignment runtime:" + str("{:10.8f}".format(self.controller.needleman.runTime)) +
							"s", font=('', 8), background="#202b3d", foreground="#ffffff").place(x=20, y=460)

		tk.Label(self, text="| Memory usage:" + str(self.controller.needleman.memSpace), font=('', 8), background="#202b3d", foreground="#ffffff").place(x=175, y=460)

		tk.Button(self, text="Back", font=('', 9,'bold'), command=lambda: self.remove(),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff").place(x=650, y=460)
		tk.Button(self, text="Save Matrices", font=('', 9, 'bold'), command=lambda: self.save("matrix"),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff").place(x=700, y=460)
		tk.Button(self, text="Save Results", font=('', 9, 'bold'), command=lambda: self.save("results"),background="#3598dc", foreground="#ffffff",activebackground="#58606b",activeforeground="#ffffff").place(x=800, y=460)

	def remove(self):
		"""
		This method removes all widgets in the resultPage frame
		"""
		for widget in self.controller.get_frame("ResultPage").winfo_children():
			widget.destroy()

		self.controller.show_frame("MainPage")

	def save(self, savetype):

		filename = tkFileDialog.asksaveasfilename(initialdir="", title="Select file",
												  filetypes=[("text files", "*.txt")])
		if filename:
			if filename[-4:] != ".txt":
				filename += ".txt"

			f = open(filename, "w")
			if f:
				if savetype is "matrix":
							"""
							This method saves all relevant matrices into user's system
							"""
							f.write("===============\n")
							f.write("SCORING MATRIX\n")
							f.write("===============\n\n")
							f.write('\n'.join([''.join(['{:>6}'.format(item) for item in row]) for row in self.controller.needleman.scoringMatrix]))
							f.write("\n\n==============================================\n")
							f.write("TRACEBACK MATRIX (1- Top, 2-Left, 3-Diagonal)\n")
							f.write("==============================================\n\n")
							f.write('\n'.join([''.join(['{:>6}'.format(item) for item in row]) for row in self.controller.needleman.tracebackMatrix]))
				elif savetype is "results":
						"""
						This method saves all relevant results into user's system
						"""
						arr = self.splitSeq()
						f.write("==========================================================\n"
							"Alignment results using BLOSUM62 with a gap score of %d" % self.controller.needleman.gapScore +
							"\n==========================================================\n\n"
							"  Optimal Alignment:\n\n")

						for i in range(0, len(arr[0])):
							f.write('\t' + str(arr[0][i]) + '\n')
							f.write('\t' + str(arr[1][i]) + '\n')
							f.write('\t' + str(arr[2][i]) + '\n\n')

						f.write("\n  Optimal Sequence Length: " + str(len(self.controller.needleman.tracebackSequences[0]))
								+ "\n\n  Identity: " + str(self.controller.needleman.alignmentResults["Identity"])
								+ "/" + str(self.controller.needleman.alignmentResults["Length"]) + ' (' +
								str(round(self.controller.needleman.alignmentResults["IdentityPercent"], 2)) + "%)"
								+ "\n  Similarity: " + str(self.controller.needleman.alignmentResults["Similarity"])
								+ "/" + str(self.controller.needleman.alignmentResults["Length"]) + ' (' +
								str(round(self.controller.needleman.alignmentResults["SimilarityPercent"], 2)) + "%)"
								+ "\n\n  Match Count: " + str(self.controller.needleman.alignmentResults["Match"])
								+ ' (' + str(round(self.controller.needleman.alignmentResults["MatchPercent"], 2)) + "%)"
								+ "\n  Mismatch Count: " + str(self.controller.needleman.alignmentResults["Mismatch"])
								+ ' (' + str(round(self.controller.needleman.alignmentResults["MismatchPercent"], 2)) + "%)"
								+ "\n  Gap Count: " + str(self.controller.needleman.alignmentResults["Gap"])
								+ ' (' + str(round(self.controller.needleman.alignmentResults["GapPercent"], 2)) + "%)"
								+ "\n\n  Alignment Path Score: " + str(self.controller.needleman.alignmentResults["PathScore"]))
				f.close()
		return



	def splitSeq(self,bp=50):
		"""
		Utility method to split sequence strings into substrings (default of 50) for writing to file
		returns an array of arrays containing the substrings on success
		"""
		str1 = ''.join(self.controller.needleman.tracebackSequences[0])
		str2 = ''.join(self.controller.needleman.tracebackSequences[1])
		str3 = ''.join(self.controller.needleman.tracebackSequences[2])

		arr1 = []
		arr2 = []
		arr3 = []
		for i in range(0, len(str1) + bp, bp):
			if i != 0:
				arr1.append(str1[i - bp:i])
				arr2.append(str2[i - bp:i])
				arr3.append(str3[i - bp:i])

		arr = [arr1, arr2, arr3]
		return arr