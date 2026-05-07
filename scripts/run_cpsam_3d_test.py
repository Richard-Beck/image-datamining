import os
import tifffile as tf
from cellpose import models

inp = os.path.expanduser("~/projects/image-datamining/data/fucci_multichannel_stacks/NCI-N87_FoF11_241016_fucci_nucleus_3ch_70z.tif")
outdir = os.path.expanduser("~/projects/image-datamining/data/cellposeSAM_test")
os.makedirs(outdir, exist_ok=True)

img = tf.imread(inp)
print("input:", img.shape, img.dtype, img.min(), img.max())

model = models.CellposeModel(gpu=True, pretrained_model="cpsam")

masks, flows, styles = model.eval(
    img,
    do_3D=True,
    z_axis=0,
    channel_axis=1,
    anisotropy=1.44
)

print("masks:", masks.shape, masks.dtype, "objects:", int(masks.max()))

mask_path = os.path.join(outdir, "NCI-N87_FoF11_cpsam_3d_masks.tif")
tf.imwrite(mask_path, masks.astype("uint32"))
print("wrote:", mask_path)
