# XeBRA MC Processor v1.0 June 2019
This is the general framework for XeBRA processor. The input of this processor are MC simulations (so far only ER).

To install this in your machine, download the repository and in your terminal type:

cmake .; make

An example of how to run this processor :
./XebraMCProcessor -d "/path/to/your/working/folder/" -i "your_root_file_without_the_extension"

Where -d is followed by your working directory that contains your input file.

This processor will generate a file called "your_root_file_without_the_extension_Sorted.root"

The principal branches generated contain:

Xp, Yp, Zp : vector - E weighted position of deposited energy per precluster
Xp_RMS, Yp_RMS, Zp_RMS : vector - E weighted position RMS per precluster 

Ed : vector - energy per precluster
nScat : number of preclusters (=scatters) (for single scatters set a cut for this at ns_scatters ==1 in your output root file)


MC Truth: 
type_primary : the mother nucleus of that decay
xpri, ypri, zpri, e_pri: Position and energy of the primary particle

