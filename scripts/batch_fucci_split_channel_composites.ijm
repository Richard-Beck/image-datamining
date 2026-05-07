output_dir = "/home/4473331/projects/image-datamining/data/fucci_max_projection_composites/";

processRoot("/share/andor_lab/Jackson/FUCCI/NCI-N87/MaxProjections/", "NCI-N87");
processRoot("/share/andor_lab/Jackson/FUCCI/MDAMB/MaxProjections/", "MDAMB");

function processRoot(root, label) {
    dirs = getFileList(root);
    for (i = 0; i < dirs.length; i++) {
        name = dirs[i];
        if (!endsWith(name, "/"))
            continue;
        name = substring(name, 0, lengthOf(name) - 1);
        if (name == "Overlays" || name == "MetaData")
            continue;

        dir = root + name + "/";
        ch00 = dir + name + "_ch00.tif";
        ch01 = dir + name + "_ch01.tif";
        ch02 = dir + name + "_ch02.tif";
        if (!File.exists(ch00) || !File.exists(ch01) || !File.exists(ch02)) {
            print("Skipping incomplete set: " + dir);
            continue;
        }

        makeComposite(ch00, ch01, ch02, output_dir + label + "_" + name + ".png");
    }
}

function makeComposite(ch00, ch01, ch02, output) {
    run("Close All");

    open(ch00);
    rename("ch00_red");
    run("Enhance Contrast", "saturated=0.35 normalize");

    open(ch01);
    rename("ch01_gray");
    run("Enhance Contrast", "saturated=0.35 normalize");

    open(ch02);
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
    print("Saved " + output);
    run("Close All");
}
