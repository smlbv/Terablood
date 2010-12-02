# -----------------------------------------------------------------------------#
# Proyecto: Terablood	        Universidad de Los Andes
# Autor: Oscar Castillo			ol.castillo28@uniandes.edu.co
# Sconstruct
# Script utilizado para construir el proyecto
# -----------------------------------------------------------------------------#
env = Environment(CXXFLAGS=['-lm','-fopenmp'])

# -----------------------------------------------------------------------------#
# Construir todas las librerias utilizadas
# NOTA: $PATH../c++/lib se debe agregar a la variable LD_LIBRARY_PATH
# -----------------------------------------------------------------------------#

env.SharedLibrary('lib/ibm', ['src/ibm.cpp'])
env.SharedLibrary('lib/lbm', ['src/lbm.cpp'])
env.SharedLibrary('lib/fem', ['src/fem.cpp'])

env.SharedLibrary('lib/mesh', ['src/mesh.cpp'])
env.SharedLibrary('lib/elemento', ['src/elemento.cpp'])

env.SharedLibrary('lib/fluid', ['src/fluid.cpp'])
env.SharedLibrary('lib/fronteras', ['src/fronteras.cpp'])

librerias = ['ibm', 'lbm', 'fem', 'mesh', 'elemento', 'fluid', 'fronteras', 'libgomp']

# -----------------------------------------------------------------------------#
# Compilacion y generacion de ejecutables
# -----------------------------------------------------------------------------#
#Program('bin/esfera', 'src/esfera.cpp', LIBS=librerias, LIBPATH='lib')
#Program('bin/poiseuille', 'src/poiseuille.cpp', LIBS=librerias, LIBPATH='lib')
#Program('bin/celula', 'src/celula.cpp', LIBS=librerias, LIBPATH='lib')
#Program('bin/cortante', 'src/cortante.cpp', LIBS=librerias, LIBPATH='lib')
env.Program('bin/treading', 'src/tank-treading.cpp', LIBS=librerias, LIBPATH='lib', )
