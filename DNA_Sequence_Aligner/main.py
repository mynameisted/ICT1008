import Tkinter as tk

from Needleman import Needleman
from MainPage import MainPage
from ResultPage import ResultPage

class App(tk.Tk):

	def __init__(self, *args, **kwargs):
		tk.Tk.__init__(self, *args, **kwargs)
		self.title("Sequence Alignment for Similarity in Biological Data")
		self.geometry("900x500")
		self.resizable(False, False)
		self.iconbitmap('favicon.ico')
		self.frames = {}
		self.needleman = Needleman()

		container = tk.Frame(self)

		container.pack(side="top", fill="both", expand=True)
		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)

		for F in (MainPage, ResultPage):
			page_name = F.__name__
			frame = F(parent=container, controller=self)
			frame.config(background="#202b3d")
			self.frames[page_name] = frame
			frame.grid(row=0, column=0, sticky="nsew")

		self.show_frame("MainPage")

	def show_frame(self, page_name):
		"""
		Switches frame to the frame specified in argument: page_name
		"""
		frame = self.frames[page_name]
		frame.tkraise()

	def get_frame(self, page_name):
		"""
		Getter method for frame of page specified in argument: page_name
		"""
		return self.frames[page_name]


if __name__ == "__main__":
	app = App()
	app.mainloop()