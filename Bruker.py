class Acqus :
	''' Clase Acqus. Contine informacion de los parametros usados en la medicion.'''
	def __init__(self) :
		self.dic = None # dictionary
		pass
	def metodo (self) :
		# returns 
		pass
class PulseProg :
	'''Clase PulseProg. Contine informacion sobre el programa de pulsos utilizado en la medicion'''
	def __init__(self) :
		self.atributo = None # 
		pass
	def metodo (self) :
		# returns 
		pass
class Crudos (Datos) :
	def __init__(self) :
		self.fid = None # Fid
		pass
	def metodo (self) :
		# returns 
		pass
class Datos (Medicion) :
	'''Un comentario sobre la clase'''
	def __init__(self) :
		self.acqus = None # Acqus
		self.procs = None # Procs
		self.pulseprog = None # PulseProg
		self.atributo = valor # tipo
		pass
	def setters (self) :
		# returns 
		pass
	def metodo (self, parametro_in, patametro_out) :
		# returns tipo
		pass
class Medicion :
	def __init__(self) :
		self.directorio = None # string
		pass
class Procesados (Datos) :
	def __init__(self) :
		self.espectro = None # Espectro
		pass
	def metodo (self) :
		# returns 
		pass
class Procs :
	'''Clase Procs. Contine informacion de los parametros usados en el procesamiento de TopSpin, se encuentran el el archivo "procs"'''
	def __init__(self) :
		self.atributo = None # 
		pass
	def metodo (self) :
		# returns 
		pass
