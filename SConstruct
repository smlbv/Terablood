# -----------------------------------------------------------------------------#
# Proyecto: Terablood	        Universidad de Los Andes
# Autor: Oscar Castillo			ol.castillo28@uniandes.edu.co
# Sconstruct
# Script utilizado para construir el proyecto
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# Construir todas las librerias utilizadas
# NOTA: $PATH../c++/lib se debe agregar a la variable LD_LIBRARY_PATH
# -----------------------------------------------------------------------------#
SharedLibrary('lib/ibm', ['src/ibm.cpp'])
SharedLibrary('lib/mesh', ['src/mesh.cpp'])
SharedLibrary('lib/elemento', ['src/elemento.cpp'])


# -----------------------------------------------------------------------------#
# Compilacion y generacion de ejecutables
# -----------------------------------------------------------------------------#
Program('bin/esfera', 'src/esfera.cpp', LIBS=['ibm','mesh', 'elemento'], LIBPATH='lib')
