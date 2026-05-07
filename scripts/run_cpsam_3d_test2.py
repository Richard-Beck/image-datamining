import os
import argparse
import numpy as np
import tifffile as tf
from scipy import ndimage as ndi
from cellpose import models


def pnorm(x, lo=1, hi=99.8, eps=1e-6):
    x = x.astype(np.float32)
    a, b = np.percentile(x, (lo, hi))
    return np.clip((x - a) / (b - a + eps), 0, 1)


def summarize_masks(masks):
    ids, counts = np.unique(masks[masks > 0], return_counts=True)
    print("objects:", len(ids))
    if len(ids):
        qs = np.quantile(counts, [0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1])
        print("voxel count quantiles:", np.round(qs, 1))
        for t in [100, 500, 1000, 2000, 5000, 10000]:
            print(f"objects < {t} voxels:", int(np.sum(counts < t)))


def save_preview_stack(img2, out_path):
    # img2 is Z,C,Y,X float 0-1. Save as uint16 for easy inspection.
    tf.imwrite(out_path, (np.clip(img2, 0, 1) * 65535).astype(np.uint16), imagej=True, metadata={"axes": "ZCYX"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=os.path.expanduser("~/projects/image-datamining/data/fucci_multichannel_stacks/NCI-N87_FoF11_241016_fucci_nucleus_3ch_70z.tif"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/projects/image-datamining/data/cellposeSAM_test/fucci_fused"))
    ap.add_argument("--anisotropy", type=float, default=1.44)
    ap.add_argument("--min_size", type=int, default=2000)
    ap.add_argument("--bf_channel", type=int, default=0)
    ap.add_argument("--fucci_red_channel", type=int, default=1)
    ap.add_argument("--fucci_green_channel", type=int, default=2)
    ap.add_argument("--invert_bf", action="store_true")
    ap.add_argument("--crop", type=int, default=0, help="Optional centered XY crop size, e.g. 512. Use 0 for full image.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    img = tf.imread(os.path.expanduser(args.input))
    print("input:", img.shape, img.dtype, img.min(), img.max())

    if img.ndim != 4:
        raise ValueError(f"Expected Z,C,Y,X image, got shape {img.shape}")

    if args.crop and args.crop < img.shape[2]:
        z, c, y, x = img.shape
        cy, cx = y // 2, x // 2
        h = args.crop // 2
        img = img[:, :, cy-h:cy+h, cx-h:cx+h]
        print("cropped:", img.shape)

    bf_raw = img[:, args.bf_channel, :, :]
    red_raw = img[:, args.fucci_red_channel, :, :]
    green_raw = img[:, args.fucci_green_channel, :, :]

    bf = pnorm(bf_raw)
    if args.invert_bf:
        bf = 1 - bf

    red = pnorm(red_raw)
    green = pnorm(green_raw)

    # FUCCI guide channel for segmentation only.
    # max() avoids cancellation and prevents one channel being invisible in mixed cell-cycle states.
    fucci_fused = np.maximum(red, green)

    # Optional mild smoothing of FUCCI guide to make it less speckly.
    fucci_fused = ndi.gaussian_filter(fucci_fused, sigma=(0.5, 0.8, 0.8))

    # Cellpose input: Z,C,Y,X
    img2 = np.stack([bf, fucci_fused], axis=1).astype(np.float32)
    print("cellpose input:", img2.shape, img2.dtype, img2.min(), img2.max())

    preview_path = os.path.join(args.outdir, "bf_plus_fucci_fused_input_ZCYX.tif")
    save_preview_stack(img2, preview_path)
    print("wrote preview input:", preview_path)

    model = models.CellposeModel(gpu=True, pretrained_model="cpsam")

    res = model.eval(
        img2,
        do_3D=True,
        z_axis=0,
        channel_axis=1,
        anisotropy=args.anisotropy,
        min_size=args.min_size
    )

    masks = res[0] if isinstance(res, tuple) else res
    print("masks:", masks.shape, masks.dtype, masks.min(), masks.max())
    summarize_masks(masks)

    tag = f"bf_fucciFused_min{args.min_size}"
    if args.invert_bf:
        tag += "_invertBF"
    if args.crop:
        tag += f"_crop{args.crop}"

    mask_path = os.path.join(args.outdir, f"NCI-N87_FoF11_cpsam_3d_{tag}_masks.tif")
    tf.imwrite(mask_path, masks.astype(np.uint32))
    print("wrote masks:", mask_path)


if __name__ == "__main__":
    main()
