input_dir = "/share/andor_lab/Jackson/FUCCI/NCI-N87/fucci_ncin87/FoF11_241016_fucci.nucleus/";
output = "/home/4473331/projects/image-datamining/data/fucci_multichannel_stacks/NCI-N87_FoF11_241016_fucci_nucleus_3ch_70z.tif";

run("Close All");
for (z = 0; z < 70; z++) {
    zstr = "" + z;
    if (z < 10)
        zstr = "0" + zstr;
    for (c = 0; c < 3; c++) {
        cstr = "" + c;
        if (c < 10)
            cstr = "0" + cstr;
        open(input_dir + "FoF11_241016_fucci.nucleus_z" + zstr + "_ch" + cstr + ".tif");
        rename("z" + zstr + "_ch" + cstr);
    }
}

run("Images to Stack", "name=FoF11_sequence title=[] use");

run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=70 frames=1 display=Composite");
Stack.setChannel(1);
run("Red");
Stack.setChannel(2);
run("Grays");
Stack.setChannel(3);
run("Green");

saveAs("Tiff", output);
run("Close All");
