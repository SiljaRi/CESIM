library(shiny)
library(fishplot) #graphic output
library(rhandsontable) #writable tables
library(shinyjs) #show and hide
library(openxlsx) #data export excel

source("util/Phylogeny.R")#code for creating a phylogeny
source("util/CheckTable.R")#code for checking phylogeny table
source("util/Mutations.R")#code for mutation creation
source("util/CNVchecks.R")#code for checking cnv table

#for colored phylogeny table
color_renderer = "
    function(instance, td, row, col, prop, value, cellProperties) {
        Handsontable.renderers.TextRenderer.apply(this, arguments);
        if (instance.params) {
            clr = instance.params.ColorCode
            clr = clr instanceof Array ? clr : [clr]
            td.style.background =  clr[row]
        }
    }"
