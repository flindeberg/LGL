######################################

All files distributed with LGL fall under the terms of the
GNU General Public License, and are copyright (c) 2002,2003 Alex Adai.

Updated 2014 Barrett Lyon and The Opte Project Team 

LGL on the web at: http://bioinformatics.icmb.utexas.edu/lgl

Much thanks to the Marcotte lab for testing.

######################################
# Table of contents

  0 Before compiling!!
  I Setup and Installation

######################################
# 0

You must have javac and jar in your path to compile JAVA programs. If you don't the
compiles will fail. (This fork removed the old build.sh, which pointed at JDK 1.5;
use the Makefile instead.)

######################################
# I

Change in to the same directory that this README is in and run `make jarfiles`
(or `make` in the repository root). Classes are compiled with `javac --release 11`.

This should produce 3 .jar files

LGLLib.jar - The full viewer library that can produce images or do the GUI 
viewer. 


For backwards compatibility the older files are also produced:


LGLView.jar - The 2D viewer to look at layouts

ImageMaker.jar - Generates high resolution images of layouts. Run without 
arguments to get the usage. This program is very crude, and will probably 
be expanded as soon as time permits. 

The source code is there so feel free to dive in!

Numerous warnings will be produced, but it should compile.

######################################
