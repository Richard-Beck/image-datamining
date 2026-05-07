input = "/share/andor_lab/Jackson/FUCCI/NCI-N87/Overlays/MaxProjections/2410_FUCCI_FoF11_241016_fucci.nucleus_Processed001/2410_FUCCI_FoF11_241016_fucci.nucleus_Processed001.tif";
output = "/home/4473331/projects/image-datamining/data_inventory/andor_lab/fucci_ncin87_fof11_overlay_max_projection.png";

open(input);
run("RGB Color");
saveAs("PNG", output);
run("Close All");
