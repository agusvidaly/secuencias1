import random
import requests
from collections import Counter
import plotly.graph_objects as go
import numpy as np
import pandas as pd
#from scipy.stats import chi2_contingency
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt
import inspect

nucleotides = ['A', 'C', 'T', 'G']
aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
codones = {codon: 0 for codon in ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']}
genetic_code = {
  "TTT": "F", "TTC": "F",  # Fenilalanina
  "TTA": "L", "TTG": "L",  # Leucina
  "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",  # Leucina
  "ATT": "I", "ATC": "I", "ATA": "I",  # Isoleucina
  "ATG": "M",  # Metionina (START)
  "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",  # Valina
  "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",  # Serina
  "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",  # Prolina
  "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",  # Treonina
  "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",  # Alanina
  "TAT": "Y", "TAC": "Y",  # Tirosina
  "TAA":  "", "TAG": "",  # STOP
  "CAT": "H", "CAC": "H",  # Histidina
  "CAA": "Q", "CAG": "Q",  # Glutamina
  "AAT": "N", "AAC": "N",  # Asparagina
  "AAA": "K", "AAG": "K",  # Lisina
  "GAT": "D", "GAC": "D",  # Ácido aspártico
  "GAA": "E", "GAG": "E",  # Ácido glutámico
  "TGT": "C", "TGC": "C",  # Cisteína
  "TGA":  "",  # STOP
  "TGG": "W",  # Triptófano
  "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",  # Arginina
  "AGT": "S", "AGC": "S",  # Serina
  "AGA": "R", "AGG": "R",  # Arginina
  "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"  # Glicina
          }
aminoacidos = {aa: 0 for aa in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']}

# 1a)
def generate_random_DNA_sequence(length): # Genera una secuencia de ntd al azar con iguales probabilidades para cada ntd
  sequence = ''
  for i in range(length):
    sequence += random.choice(nucleotides)
  return sequence

# 1a)
def generate_random_DNA_weighted(): # Genera una secuencia  (string) en la que el contenido GC lo define el usuario.
  # input= La longitud de la secuencia y la frecuencia de GC
  # output= string de la secuencia
  length=int(input('Seleccione el largo de la secuencia de nt: '))
  CG=float(input("Seleccione la frecuencia GC (entre 0 y 1): ")) #frecuencia de CG (entre 0 y 1)
  AT=1-CG #frecuencia de AT (se calcula sola)
  sequence=""
  frec=[AT/2,CG/2,AT/2,CG/2] # frecuencias de cada par de ntds
  for i in range(length):
    sequence += np.random.choice(nucleotides,p=frec) # agrega nucleotidos con la probabilidad definida
#  print("la secuencia ponderada es:",sequence)
  return(sequence)

# 1b)
def generate_random_protein_sequence(length):
  sequence_protein = ''
  for i in range(length):
    sequence_protein += random.choice(aminoacids)
  return sequence_protein

# 1b)
probabilidades=[0.00,0.1,0.00,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]

def generate_protein_sequence_weighted(): # Genera un str de secuencia aminoacídica donde todas los aá tienen != probabilidad. Las frecuencias se definen arriba
  lista_num_azar=[]
  numeros=list(range(0,20))
  sec_prot=""
  len_sec_prot=int(input("Seleccione el largo de la secuencia de aminoácidos: "))
  for i in range(len_sec_prot):
    lista_num_azar.append(np.random.choice(numeros,p=probabilidades))
  for i in lista_num_azar:
    sec_prot+=aminoacids[int(i)-1]
  return(sec_prot)

# 1c) i)
def frecBases(sequence):
  lista_frec=[]
  cod_gen_unicos=set(sequence)
  lista_frec.append(sequence.count("A")*100/len(sequence))
  lista_frec.append(sequence.count("C")*100/len(sequence))
  lista_frec.append(sequence.count("T")*100/len(sequence))
  lista_frec.append(sequence.count("G")*100/len(sequence))
  return(lista_frec)

def graficoFrecBases(sequence,name=None): # Gráfico de frecuencias de cada base
  frec_bases=frecBases(sequence)
  if name is None:
    frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
    name = [name for name, value in frame_info.f_locals.items() if value is sequence][0]
  x_values = nucleotides # Valores para el eje X
  y_values = frec_bases # Valores para el eje Y
  fig = go.Figure(data=go.Bar(x=x_values, y=y_values,marker_color='#576EF4'))
  fig.update_layout(
    title_text=f'Frecuencias de las bases en la secuencia: {name}',
    xaxis_title='Bases',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 100]),  # Para que el eje y vaya de 0 a 100%
    template='plotly')
  fig.show()

# 1c) ii)
def frecCod(sequence,name=None): # Separa los codones y los cuenta. Devuelve
  for cd in range(0,len(sequence)-2,3):
    if sequence[cd:cd+3] in codones :
      codones[sequence[cd:cd+3]] +=1 #Cuenta las veces que aparece cada codon.
    codonesTot =sum(list(codones.values()))
    porcentajesCod = [(codon*100)/codonesTot for codon in list(codones.values())] #Valores para el eje Y.
    y_values_cod = porcentajesCod
  return(y_values_cod)

def graficoFredCod(sequence,name=None): # Gráfico de frecuencias de cada codon
  if name is None:
    frame_info = inspect.currentframe().f_back
    name = [name for name, value in frame_info.f_locals.items() if value is sequence][0]
  y_values_cod=frecCod(sequence)
  x_values_cod = list(codones.keys())
  fig = go.Figure(data=go.Bar(x=x_values_cod, y=y_values_cod,marker_color='#576EF4'))
  fig.update_layout(
    title_text=f"Porcentaje de Cada Codon en la Secuencia: {name}",
    xaxis_title='Codones',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 10]),  # Para que el eje y vaya de 0 a 100%
    template='plotly')
  fig.show()

# 1c) iii)
def translate(sequence): # Traduce una secuencia de ntd a una de aa segun el codigo establecido arriba
  global genetic_code
  amino_acid_sequence_t = []
  for i in range(0, len(sequence) - 2, 3):
    codon = sequence[i:i+3]
    amino_acid = genetic_code.get(codon)
    amino_acid_sequence_t.append(amino_acid) # es la lista con aa de la secuencia dada
    aa_seq_t=''.join(amino_acid_sequence_t)
  return(aa_seq_t)

def graficoFrecAa(aa_sequence,name=None):
  # Calcula la frecuencia de aa en una secuencia dada.
  # Toma en cuenta si la secuencia ingresada es de ntd (en cuyo caso la traduce) o si son aa.
  if name is None:
    frame_info = inspect.currentframe().f_back
    name = [name for name, value in frame_info.f_locals.items() if value is aa_sequence][0]
  if all(char in nucleotides for char in aa_sequence):
    print('perá que hay que traducir esta secuencia de nucleótidos')
    aa_sequence_translated=translate(aa_sequence)
    graficoFrecAa(aa_sequence_translated)
  else:
    for a in range(len(aa_sequence)):
      if aa_sequence[a]=='X' or aa_sequence[a]=='Z':
          continue
      else:
        aminoacidos[aa_sequence[a]] += 1 #Cuenta las veces que aparece cada aa
    x_values_aa = list(aminoacidos.keys())
    aminoacidosTot =sum(list(aminoacidos.values()))
    porcentajesAA = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
    # Ordenar los datos en función del porcentaje
    sorted_data = sorted(zip(x_values_aa, porcentajesAA), key=lambda pair: pair[1])
    # Separar las listas ordenadas
    x_values_sorted, y_values_sorted = zip(*sorted_data)
    y_values_aa_sorted = porcentajesAA
    fig = go.Figure(data=go.Bar(x=x_values_sorted, y=y_values_sorted,marker_color='#576EF4'))
    fig.update_layout(
      title_text=f'Porcentaje de Cada Aminoácido en la Secuencia: {name}',
      xaxis_title='Aminoácidos',
      yaxis_title='Porcentaje (%)',
      yaxis=dict(range=[0, 10]),  # Ajusta el rango del eje y si es necesario
      template='plotly')
    fig.show()

def process_sequence(sequence,name=None):
  if all(char in nucleotides for char in sequence):
    print("Processing nucleotide sequence...")
    graficoFrecBases(sequence,name)
    graficoFredCod(sequence,name)
    graficoFrecAa(sequence,name)
  else:
    print("Processing amino acid sequence...")
    graficoFrecAa(sequence,name)

# 1d)
def count_orf_lengths(secuencia,name=None): # output: una lista con las longitudes de los orfs
  stop = ['TAA', 'TAG', 'TGA']
  for frame in range(3):
    contador=  0
    largoORF=[]
    for cod in range(frame,len(secuencia)-2,3):
      if secuencia[cod:cod+3] in stop:
        largoORF.append(contador)
        contador = 0 # si encuentra un codon de stop resetea el contador a 0
      else:
        contador +=1
  return(largoORF)

def graficoOrfs(secuencia,name=None):
  if name is None:
    frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
    name = [name for name, value in frame_info.f_locals.items() if value is secuencia][0]
  y_values=count_orf_lengths(secuencia,name=None) # largo de los distintos ORF
  media=np.mean(y_values)
  fig=go.Figure()
  fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  fig.update_layout(title_text = f'Secuencia: {name}')
  fig.update_xaxes(title_text="Longitud ORF (codones)")
  fig.update_yaxes(title_text="Frecuencia")
  fig.show()

def graficoOrfs3ML(seq,name=None): #cuenta los orfs y hace un histograma para cada marco de lectura
    if name is None:
      frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
      name = [name for name, value in frame_info.f_locals.items() if value is seq][0]
    stop = ['TAA', 'TAG', 'TGA']
    for frame in range(3):
      contador=  0
      largoORF=[]
      for cod in range(frame,len(seq)-2,3):
          if seq[cod:cod+3] in stop:
              largoORF.append(contador)
              contador = 0 # si encuentra un codon de stop resetea el contador a 0
          else:
              contador +=1 # si no es un codon de stop suma de a 1 al contador.
      y_values=largoORF # largo de los distintos ORF
      media=np.mean(y_values)
      fig=go.Figure()
      fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
      fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
      fig.update_layout(title_text = f'Marco de lecura - {frame+1}, Secuencia: {name}')
      fig.update_xaxes(title_text="Longitud ORF (codones)")
      fig.update_yaxes(title_text="Frecuencia")
      fig.show()
    return # muestra los tres 3 histogramas

def graficoOrfs6ML(seq,name=None):
  if name is None:
      frame_info = inspect.currentframe().f_back
      name = [name for name, value in frame_info.f_locals.items() if value is seq][0]
  graficoOrfs3ML(seq,name)
  reverse_complement=Seq(seq).reverse_complement()
  graficoOrfs3ML(str(reverse_complement),name)


# 2a)
def graficoFrecBases_comparacion(seq1,seq2):
  serie1=frecBases(seq1)
  serie2=frecBases(seq2)
  fig = go.Figure()
  fig.add_trace(go.Bar(
    x=nucleotides,
    y=serie1,
    name='secuencia',
    marker_color='blue'))
  fig.add_trace(go.Bar(
    x=nucleotides,
    y=serie2,
    name='reversa complementaria',
    marker_color='lightskyblue'))
  fig.update_layout(
    title='Frecuencias de las bases en la secuencia y la reversa complementaria',
    xaxis_title='Bases',
    yaxis_title='Frecuencias (%)',
    barmode='group')

  fig.show()

def graficoFredCod_comparacion(seq1,seq2):
  y_values_cod1=frecCod(seq1)
  y_values_cod2=frecCod(seq2)
  x_values_cod = list(codones.keys())
  fig = go.Figure()
  fig.add_trace(go.Bar(
    x=x_values_cod,
    y=y_values_cod1,
    name='secuencia',
    marker_color='blue'))
  fig.add_trace(go.Bar(
    x=x_values_cod,
    y=y_values_cod2,
    name='reversa complementaria',
    marker_color='lightskyblue'))
  fig.update_layout( 
    title='Frecuencias de codones de la secuencia y la reversa complementaria',
    xaxis_title='Codones',
    yaxis_title='Frecuencias',
    barmode='group')
  
  fig.show()

def graficoFrecAa_comparacion(seq1,seq2):
  # prot al azar
  for a in range(len(seq1)):
    if seq1[a]=='X' or seq1[a]=='Z':
          continue
    else:
      aminoacidos[seq1[a]] += 1 #Cuenta las veces que aparece cada aa
  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_azar = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
  # prot real
  for a in range(len(seq2)):
    if seq2[a]=='X' or seq2[a]=='Z':
      continue
    else:
        aminoacidos[seq2[a]] += 1 #Cuenta las veces que aparece cada aa
  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_real = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.
  fig = go.Figure()
  # Agregar barras para los aminoácidos de la secuencia al azar
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_azar,
    name='Secuencias',
    marker_color='blue'))
  # Agregar barras para los aminoácidos de las secuencias reales
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_real,
    name='Reversa complementaria',
    marker_color='lightskyblue'))
  # Configurar el layout del gráfico
  fig.update_layout(
    title='Comparación de Porcentajes de Aminoácidos en Secuencia de Proteína y Secuencia al Azar',
    xaxis_title='Aminoácidos',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 10]),  # Ajusta el rango del eje Y si es necesario
    barmode='group',  # Agrupar barras lado a lado
    template='plotly')

  fig.show()

def comparacionEjercicio2a(lenght):
  secuencia=generate_random_DNA_sequence(lenght)
  reversa_complementaria=Seq(secuencia)
  reversa_complementaria=reversa_complementaria.reverse_complement()
  secuencia=Seq(secuencia)
  secuencia_traducida=translate(secuencia)
  reversa_complementaria_traducida=translate(reversa_complementaria)
  graficoFrecBases_comparacion(secuencia,reversa_complementaria)
  graficoFredCod_comparacion(secuencia,reversa_complementaria)
  graficoFrecAa_comparacion(secuencia_traducida,reversa_complementaria_traducida)

# 2b)
def graficoOrfs_comparacion(length):
  secuencia=generate_random_DNA_sequence(length)
  reversa_complementaria=Seq(secuencia)
  reversa_complementaria=reversa_complementaria.reverse_complement() 
  # seq 1
  serie1=count_orf_lengths(secuencia)
  media=np.mean(serie1)
  fig = go.Figure()
  fig.add_trace(go.Histogram(
    x=serie1,
    nbinsx=300,
    name='Secuencia',
    marker_color='blue',
    opacity=0.55))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # seq 2 
  serie2=count_orf_lengths(reversa_complementaria)
  media=np.mean(serie2)
  fig.add_trace(go.Histogram(
    x=serie2,
    nbinsx=300,
    name='Reverse Complement ORF Lengths',
    marker_color='lightskyblue',
    opacity=0.65))
  fig.add_vline(x=media, line_dash="dash", line_color="violet", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # juntos
  fig.update_layout(
    title='Distribution of ORF Lengths Between Stop Codons',
    xaxis_title='ORF Length',
    yaxis_title='Frequency',
    barmode='overlay')
  fig.show()


# 2b)
# URL del archivo GenBank
archivo_gb = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/sequence.gb'

def atributos(url):
  # Descarga el archivo desde la URL
  response = requests.get(url)
  response.raise_for_status()  # Lanza un error si la descarga falla

  # Guarda el contenido del archivo en un objeto StringIO
  from io import StringIO
  file_content = StringIO(response.text)

  # Lee el archivo GenBank desde el contenido descargado
  record = SeqIO.read(file_content, "genbank")
  
  print(f"ID: {record.id}")
  print(f"Name: {record.name}")
  print(f"Description: {record.description}")
  print(f"Sequence: {record.seq[:100]}...")  # Muestra los primeros 100 nucleótidos

# 2c)
def abrir_archivo(url): #abre el url, genera una lista con cada secuencia y devuelve un str con todas las secuencias seguidas.
  responseNt = requests.get(url)
  secuencia = responseNt.text #lo transforma en archivo de texto
  archivo = secuencia.splitlines()
  lista_sec = []
  secuencia_parcial = ''
  for line in range(len(archivo)):
    if archivo[line].startswith('>'):
      nro_linea = line+1
      while not archivo[nro_linea].startswith('>') and nro_linea+1 != len(archivo):
        seq = archivo[nro_linea].strip()
        secuencia_parcial +=seq
        nro_linea += 1
      lista_sec.append(secuencia_parcial)
  print('se leyeron las',len(lista_sec),'secuencias')
  return(lista_sec)

# 3a) i)
def histogramasOrf(secuencia,name=None):
  if name is None:
    frame_info = inspect.currentframe().f_back # para obtener el nombre de la variable y pasarla al titulo del histograma
    name = [name for name, value in frame_info.f_locals.items() if value is secuencia][0]
  longitud_ORF=count_orf_lengths(secuencia)
  y_values = longitud_ORF
  media =np.mean(y_values)
  fig = go.Figure()
  fig.add_trace(go.Histogram(x=y_values,nbinsx=300))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  fig.update_xaxes(title_text="Longitud ORF (codones)")
  fig.update_yaxes(title_text="Frecuencia")
  fig.update_layout(title=f'Longitudes de ORFs en la secuencia {name}')
  fig.show() #histgrama de ORFs

def histogramasJuntos(sec_cRNA,sec_RNA_azar):
  serie_azar=count_orf_lengths(sec_RNA_azar)
  serie_real=count_orf_lengths(sec_cRNA)
  # DNA al azar
  y_values=serie_azar
  media =np.mean(y_values)
  fig=go.Figure()
  fig.add_trace(go.Histogram(
    x=serie_azar,
    name='Secuencia de nucleótidos al azar',
    opacity=0.5,
    marker_color='blue'))
  fig.add_vline(x=media, line_dash="dash", line_color="violet", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # cRNA
  y_values=serie_real
  media =np.mean(y_values)
  fig.add_trace(go.Histogram(
    x=serie_real,
    name='cRNA',
    opacity=0.5,
    marker_color='green'))
  fig.add_vline(x=media, line_dash="dash", line_color="red", annotation_text=f"Media = {media:.2f}", annotation_position="top right")
  # Superponer los histogramas en lugar de apilarlos
  fig.update_layout(
    barmode='overlay',  # Superponer las barras
    title_text='Histograma con dos series',  # Título del gráfico
    xaxis_title_text='Valor',  # Etiqueta del eje X
    yaxis_title_text='Frecuencia')  # Etiqueta del eje Y

  fig.show()

# 3a) ii)
def FrecAaJuntos(sec_aa_real,sec_aa_azar):
  # prot al azar
  for a in range(len(sec_aa_azar)):
    if sec_aa_azar[a]=='X' or sec_aa_azar[a]=='Z' or sec_aa_azar[a]=='B':
          continue
    else:
      aminoacidos[sec_aa_azar[a]] += 1 #Cuenta las veces que aparece cada aa

  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_azar = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.

  # prot real
  for a in range(len(sec_aa_real)):
    if sec_aa_real[a]=='X' or sec_aa_real[a]=='Z' or sec_aa_real[a]=='B':
      continue
    else:
        aminoacidos[sec_aa_real[a]] += 1 #Cuenta las veces que aparece cada aa

  x_values_aa = list(aminoacidos.keys())
  aminoacidosTot =sum(list(aminoacidos.values()))
  porcentajesAA_real = [(aa*100)/aminoacidosTot for aa in list(aminoacidos.values())] #Valores para el eje Y.

  fig = go.Figure()
  # Agregar barras para los aminoácidos de la secuencia al azar
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_azar,
    name='Secuencias de Proteínas Reales',
    marker_color='blue'))

  # Agregar barras para los aminoácidos de las secuencias reales
  fig.add_trace(go.Bar(
    x=x_values_aa,
    y=porcentajesAA_real,
    name='Secuencia al Azar',
    marker_color='lightskyblue'))

  # Configurar el layout del gráfico
  fig.update_layout(
    title='Comparación de Porcentajes de Aminoácidos en Secuencia de Proteína y Secuencia al Azar',
    xaxis_title='Aminoácidos',
    yaxis_title='Porcentaje (%)',
    yaxis=dict(range=[0, 10]),  # Ajusta el rango del eje Y si es necesario
    barmode='group',  
    template='plotly')

  fig.show()



# 1) histograma de las longitudes de los orfs
#sec_cRNA=abrir_archivo(archivo_cRNA)
#sec_cRNA=''.join(sec_cRNA)
#longitud_azar=len(sec_cRNA)
#sec_RNA_azar = generate_random_DNA_sequence(longitud_azar)
#histogramasJuntos(sec_cRNA,sec_RNA_azar) #toma la secuencia del urlRNAc y una hecha al azar
#histogramasOrf(sec_RNA_azar)
#histogramasOrf(sec_cRNA)

# 2) graficos con las frecuencias de los codones
#graficoFredCod(sec_RNA_azar)
#graficoFredCod(sec_cRNA)

# 3) graficos con las frecuencias de los aa
#sec_aa_real=abrir_archivo(archivo_prot)
#sec_aa_real=''.join(sec_aa_real)
#longitud_azar_prot=len(sec_aa_real)
#sec_aa_azar = generate_random_protein_sequence(longitud_azar_prot) 
#graficoFrecAa(sec_aa_azar,'secuencias de aa generadas al azar') 
#graficoFrecAa(sec_aa_real,'35 proteínas reales elegidas al azar')
#FrecAaJuntos(sec_aa_real,sec_aa_azar)

# 4) graficos con las frecuencias de los aa de secuencias traducidas
#sec_RNA=abrir_archivo(archivo_RNA)
#sec_RNA=''.join(sec_RNA)
#longitud_azar=len(sec_RNA)
#sec_RNA_azar = generate_random_DNA_sequence(longitud_azar)
#graficoFrecAa(sec_RNA_azar) #toma la secuencia de url_para_traducir y hace una de nts al azar
#graficoFrecAa(sec_RNA)
#FrecAaJuntos(sec_RNA_azar,sec_RNA)

def largo_proteinas_reales(secuencias):
  print(f'Se usaron {len(secuencias)} secuencias de proteinas reales.')
  lista_longitudes=[]
  for sec in secuencias:
    lista_longitudes.append(len(sec))
  y_values = lista_longitudes # largo de las distintas proteinas

  # histograma
  median = np.median(y_values) 
  fig = go.Figure()
  histogram = go.Histogram(x=y_values, nbinsx=300)
  fig.add_trace(histogram)
  fig.add_vline(x=median, line_dash="dash", line_color="red", annotation_text=f"Median = {median:.2f}", annotation_position="top right")
  fig.update_xaxes(title_text="Longitud en aminoácidos")
  fig.update_yaxes(title_text="Frecuencia")
  return fig.show()


# Example ussage

# URL de la secuencia FASTA con 10500 proteínas
urlProteinas = 'https://raw.githubusercontent.com/Nehuenpg/Fasta-Sequences/main/SecuenciaProteinas'

responseProt = requests.get(urlProteinas)
fasta_data_prot = responseProt.text #lo transforma en archivo de texto
archivoProt = fasta_data_prot.splitlines()


#lista_prot=abrir_archivo(archivoProt)
#print(largo_proteinas_reales(lista_prot))

def generar_proteinas_azar(n_proteinas=10500, media_log=5.8, std_log=0.6): 
  """
  Genera secuencias de proteínas de longitud aleatoria basada en una distribución log-normal.
  Parámetros:
  - n_proteinas: Número de proteínas a generar.
  - media_log: Media logarítmica de la longitud (log-normal).
  - std_log: Desviación estándar logarítmica de la longitud.
  Retorna:
  - Una lista con las secuencias de proteínas generadas.
  - Una lista con las longitudes de las proteínas generadas.
  """
    
  longitudes = np.random.lognormal(mean=media_log, sigma=std_log, size=n_proteinas).astype(int) # Generar longitudes a partir de una distribución log-normal
  secuencias_azar = []
  #amino_acidos = 'ACDEFGHIKLMNPQRSTVWY'
  for longitud in longitudes:
    secuencia = ''.join(random.choices(aminoacids, k=longitud)) #cambié amino_acidos por la lista aminoacids definida al ppio
    secuencias_azar.append(secuencia)

  print(f'Se usaron {len(longitudes)} secuencias de proteinas generadas al azar.')
  y_values = longitudes # largo de las distintas proteinas
  median = np.median(y_values)

  # Crear el histograma
  fig = go.Figure()
  histogram = go.Histogram(x=y_values, nbinsx=300)
  fig.add_trace(histogram)

  # Añadir una línea vertical para la media
  fig.add_vline(x=median, line_dash="dash", line_color="red", annotation_text=f"Median = {median:.2f}", annotation_position="top right")
  # Actualizar los ejes
  fig.update_xaxes(title_text="Longitud en aminoácidos")
  fig.update_yaxes(title_text="Frecuencia")

  # Mostrar el gráfico
  return fig.show()

#print(generar_proteinas_azar(n_proteinas=10500, media_log=5.8, std_log=0.6))


#FrecAaJuntos(secuencia_traducida,reversa_complementaria_traducida)