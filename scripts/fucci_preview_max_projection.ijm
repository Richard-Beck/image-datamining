input_dir = "/share/andor_lab/Jackson/FUCCI/NCI-N87/MaxProjections/2410_FUCCI_FoF11_241016_fucci.nucleus_Processed001/";
prefix = "2410_FUCCI_FoF11_241016_fucci.nucleus_Processed001";
output = "/home/4473331/projects/image-datamining/data_inventory/andor_lab/fucci_ncin87_fof11_max_projection_rgb.png";

open(input_dir + prefix + "_ch00.tif");
rename("ch00_red");
run("Enhance Contrast", "saturated=0.35 normalize");

open(input_dir + prefix + "_ch01.tif");
rename("ch01_gray");
run("Enhance Contrast", "saturated=0.35 normalize");

open(input_dir + prefix + "_ch02.tif");
rename("ch02_green");
run("Enhance Contrast", "saturated=0.35 normalize");

run("Merge Channels...", "c1=ch00_red c2=ch01_gray c3=ch02_green create keep");
rename("fucci_composite");

Stack.setChannel(1);
run("Red");
Stack.setChannel(2);
run("Grays");
Stack.setChannel(3);
run("Green");

run("RGB Color");
saveAs("PNG", output);
run("Close All");
