input = "/home/4473331/projects/image-datamining/data/fucci_multichannel_stacks/NCI-N87_FoF11_241016_fucci_nucleus_3ch_70z.tif";

open(input);
getDimensions(width, height, channels, slices, frames);
print("width=" + width);
print("height=" + height);
print("channels=" + channels);
print("slices=" + slices);
print("frames=" + frames);
print("bitDepth=" + bitDepth());
run("Close All");
