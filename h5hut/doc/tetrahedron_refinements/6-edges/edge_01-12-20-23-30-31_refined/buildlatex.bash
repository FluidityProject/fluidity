#!/bin/bash

function write_file() {
	local	-r	fname_texfile=$1
	local	-r	fname_image=$2
	local	-r	annotation=$3

	cat > "${fname_texfile}" <<-EOF
\documentclass[a4paper,12pt]{article}

\usepackage{graphicx}
\thispagestyle{empty}
\textwidth 190mm
\usepackage[left=10mm,top=10mm,right=10mm,nohead,nofoot]{geometry}
\begin{document}
\includegraphics[width=180mm]{${fname_image}}

${annotation}
\end{document}
EOF
}

function compile_texfile() {
	local	-r	fname_texfile=$1

	pdflatex "${fname_texfile}"
}

fname_image="tets.png"
fname_texfile=${fname_image%.*}.tex
annotation=${PWD##*/}
annotation=${annotation//[a-z_]}
annotation="Refined edges: ${annotation//-/ }"

write_file "${fname_texfile}" "${fname_image}" "${annotation}"
compile_texfile "${fname_texfile}"
