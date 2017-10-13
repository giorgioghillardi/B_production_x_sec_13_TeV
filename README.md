# B_production_x_sec_13_TeV
Files for studies on B particle (see the list in interface/channels.h)

MADE IN LIP OF LISBOA 


The analysis start from the output of Bfinder the other project in my repository (https://github.com/giorgioghillardi/Bfinder)

In "bin" there are the principal macroes that does all the steps to rech cross sections and branching fractios

In "interface" there are headers files with the tools for the macros in "bin"




- myloop_new.cc

is the main code that do the cuts needed for the anlysis.
The command work likes
myloop_new --channel 1 --mc 0 --truth 0 --cuts 1 --debug 0 --output somewhere

select the channel form interface/channels.h
select if are data or MonteCarlo (--mc 0 or 1)
select in case of MC if you want the Gen check
select if you want cout debugging
select where you want the outpu


- myloop_gen


-calculate_bin_syst.cc

--channel 1 --syst cbpdf  --ptmin 0 --ptmax 150 --ymin 0.00 --ymax 2.25

indicate the systematics error (signal_pdf or mass_window or combined_syst)
