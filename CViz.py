import tkinter as tk
from tkinter import Tk, Button, Checkbutton, BooleanVar, Canvas, Toplevel, filedialog
from tkinter.ttk import Combobox, Notebook, Frame
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fasta_data_dict = {}
snp_positions = {}
canvas = None
dataset_dropdowns = []


def load_fasta_files(file_paths):
    global fasta_data_dict
    fasta_data_dict.update(
        {record.id: record.seq for file_path in file_paths for record in SeqIO.parse(file_path, "fasta")})
    print("Loaded FASTA files:", file_paths)


def refill_dataset_dropdowns():
    for dataset_dropdown in dataset_dropdowns:
        dataset_dropdown["values"] = list(fasta_data_dict.keys())


def open_file_dialog():
    file_paths = filedialog.askopenfilenames(
        filetypes=[("Fasta Files", ("*.fasta", "*.fa", "*.fna", "*.ffn", "*.faa", "*.frn"))])
    if file_paths:
        load_fasta_files(file_paths)
        refill_dataset_dropdowns()
        create_visualization_window()


def display_chromosome_3d(chromosome_id, highlight_snps=False):
    sequence = fasta_data_dict.get(chromosome_id)
    if not sequence:
        print(f"Error: Chromosome {chromosome_id} not found in FASTA data")
        return

    t = np.linspace(0, 2 * np.pi, 100)
    x, y, z = 10 * (1 + np.sin(t)), 10 * (t + 2 * np.pi * np.cos(t)), 10 * t

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z)
    ax.set_title(f'Chromosome {chromosome_id} Visualization')
    ax.set_xlabel('X-axis'), ax.set_ylabel('Y-axis'), ax.set_zlabel('Z-axis')
    plt.show()


def on_dropdown_select(event, dropdown_var, highlight_snps):
    display_chromosome_3d(dropdown_var.get(), highlight_snps)


def update_visualization(highlight_snps):
    for dataset_dropdown in dataset_dropdowns:
        on_dropdown_select(None, dataset_dropdown, highlight_snps)


def create_visualization_window():
    global canvas, dataset_dropdowns
    vis_window = Toplevel()
    vis_window.title("Chromosome Visualization")
    highlight_snps = BooleanVar()
    snp_checkbox = Checkbutton(vis_window, text="Highlight SNPs", variable=highlight_snps,
                               command=lambda: update_visualization(highlight_snps.get()))
    snp_checkbox.pack()

    dataset_dropdowns = []
    for _ in range(4):
        dataset_var = tk.StringVar(vis_window)
        dataset_dropdown = Combobox(vis_window, textvariable=dataset_var)
        dataset_dropdown.bind("<<ComboboxSelected>>",
                              lambda event, var=dataset_var: on_dropdown_select(event, dataset_var, highlight_snps))
        dataset_dropdown.pack(side=tk.LEFT, padx=5)
        dataset_dropdowns.append(dataset_dropdown)

    chromosome_dropdown = Combobox(vis_window, textvariable=tk.StringVar(vis_window))
    chromosome_dropdown.place(relx=0.5, rely=0.1, anchor=tk.CENTER)
    fasta_button = Button(vis_window, text="Upload FASTA File", command=open_file_dialog)
    fasta_button.pack()
    canvas = Canvas(vis_window, bg="white", width=500, height=500)
    canvas.pack()
    tab_control = Notebook(vis_window)
    extra_tab = Frame(tab_control)
    tab_control.add(extra_tab, text="Extra Tab")
    tab_control.pack()
    print("Visualization window created")


if __name__ == "__main__":
    start_screen = Tk()
    open_file_button = Button(start_screen, text="Open FASTA Files", command=open_file_dialog)
    open_file_button.pack()
    start_screen.mainloop()
